# This material is part of the SAVI-Space created by Malte Korn and Matthias Rarey
# at the Center for Bioinformatics, University of Hamburg, with support from
# Marc Nicklaus (NIH,NCI), Phil Judson, Raphael Klein (BioSolveIT GmbH) and
# Christian Lemmen (BioSolveIT GmbH).
# Philip Judson and Marc Nicklaus provided the LHASA transform rules and assisted with their knowledge about SAVI-2020.
# SAVI-Space and all its components including this file/directory is licensed under CC-BY 4.0.

from tqdm import tqdm
import pandas as pd
from pathlib import Path
import logging
import json
import datetime

from rdkit import Chem

from utils import load_reactions, load_kill_patterns, load_pattern_json, get_filtered_bbs_path, get_pattern_smarts, mask_exocyclic_doublebonds, smarts_exocyclic_doublebonds, initialize_logger

def rdkit_check(smiles, smarts, disallow_exocyclic_doublebonds=False, uniquify=False):
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return 0  # Invalid SMILES, return early

    if len(smarts) == 2:
        unique_matches = smarts[1]
        smarts = smarts[0]
    else:
        unique_matches = None

    if disallow_exocyclic_doublebonds:
        mol = mask_exocyclic_doublebonds.normalize(mol)
        if not mol:
            logging.error("Failed to normalize molecule")
            raise ValueError("Failed to normalize molecule")
        smarts = smarts_exocyclic_doublebonds(smarts)

    patt = Chem.MolFromSmarts(smarts)
    if not patt:
        logging.error(f"Invalid SMARTS pattern: {smarts}")
        raise ValueError(f"Invalid SMARTS pattern: {smarts}")

    patt_atom_mapping = [atom.GetAtomMapNum() for atom in patt.GetAtoms()]
    matches = mol.GetSubstructMatches(patt, uniquify=uniquify)
    valid_matches = []

    if matches:
        for match in matches:
            if unique_matches:
                if unique_matches == [-1]:
                    match_candidate = [-1]
                else:
                    match_candidate = [m for (m, a) in zip(match, patt_atom_mapping) if a in unique_matches]
            else:
                match_candidate = [m for (m, a) in zip(match, patt_atom_mapping) if a != 0]
            if not match_candidate:
                logging.error("No atom mapping found in match")
                raise ValueError("No atom mapping found in match")
            if match_candidate not in valid_matches:
                valid_matches.append(match_candidate)

    return len(valid_matches)


def rdkit_pattern_matching(bbs_path, outpath, reaction_index, reactions, reaction_pattern, disallow_exocyclic_doublebonds=False, unique_matches=None, uniquify=False):
    reaction_outpath = get_filtered_bbs_path(outpath, reaction_index, reaction_pattern)
    reaction_outpath.mkdir(parents=True, exist_ok=True)
    pattern_smarts = get_pattern_smarts(reactions, reaction_index, unique_matches)
    df = pd.read_csv(bbs_path, sep=" ", header=None, names=["smiles", "name"])
    for name, patt in pattern_smarts.items():
        df[name] = df["smiles"].apply(lambda x: rdkit_check(x, patt, disallow_exocyclic_doublebonds=disallow_exocyclic_doublebonds, uniquify=uniquify))
    df.to_csv(reaction_outpath / f"results.csv")


def generate_sublist_pairs(series1, series2, compatibility_dict, kill_indices):
    # Group elements of series1 and series2 by their group indices
    group_dict1 = series1.groupby(series1).apply(lambda x: x.index.tolist()).to_dict()
    group_dict2 = series2.groupby(series2).apply(lambda x: x.index.tolist()).to_dict()

    result = []
    kill_idx = []

    # Generate pairs of sublists based on compatibility
    for g1 in group_dict1:
        if g1 not in compatibility_dict:
            continue
        compatible_groups = compatibility_dict[g1]
        for g2 in compatible_groups:
            if g2 in group_dict2:
                new_pair = (group_dict1[g1], group_dict2[g2])
                # Check if new pair is a subset of any existing pair
                is_subset = False
                for pair in result:
                    if set(new_pair[0]).issubset(pair[0]) and set(new_pair[1]).issubset(pair[1]):
                        is_subset = True
                        break
                    elif set(pair[0]).issubset(new_pair[0]) and set(pair[1]).issubset(new_pair[1]):
                        result.remove(pair)
                        break
                if not is_subset:
                    result.append(new_pair)
                    kill_idx.append(kill_indices[g1])
    return result, kill_idx

def kill_statements(reaction_index, outpath, kill_patterns, reactions, reaction_pattern, disallow_exocyclic_doublebonds=False):
    reaction_outpath = get_filtered_bbs_path(outpath, reaction_index, reaction_pattern)
    df = pd.read_csv(reaction_outpath / f"results.csv", index_col=0)
    kill_pattern = kill_patterns[reaction_index]
    pattern_smarts = get_pattern_smarts(reactions, reaction_index)
    supbar = tqdm(enumerate(pattern_smarts), desc="reactant", leave=False)
    kill_df = pd.DataFrame()
    name_dict = {}
    for idx1, pattern_name in supbar:
        if pattern_name.endswith("_r2"):
            continue
        logging.info(f"Applying kill function for {pattern_name[:-3].replace('_', ' ')}")
        reaction_idx, pattern_idx, sub_pattern_idx, reactant = pattern_name.split("_")
        if f"{reaction_idx}_{pattern_idx}" not in name_dict:
            name_dict[f"{reaction_idx}_{pattern_idx}"] = {}
        name_dict[f"{reaction_idx}_{pattern_idx}"][f"{sub_pattern_idx}"] = []
        temp_name = "_".join([reaction_idx, pattern_idx, sub_pattern_idx])
        if (df[temp_name + "_r1"]==1).sum() > 0:
            kill_output1 = df[df[temp_name + "_r1"]>0]["smiles"].apply(lambda x: kill(x, pattern_smarts[temp_name + "_r1"], kill_pattern, disallow_exocyclic_doublebonds=disallow_exocyclic_doublebonds))
        else:
            kill_output1 = pd.Series()
        if (df[temp_name + "_r2"]==1).sum() > 0:
            kill_output2 = df[df[temp_name + "_r2"]>0]["smiles"].apply(lambda x: kill(x, pattern_smarts[temp_name + "_r2"], kill_pattern, disallow_exocyclic_doublebonds=disallow_exocyclic_doublebonds))
        else:
            kill_output2 = pd.Series()
        compatibility_dict = {}
        temp_kill_indices = {}
        for adata in kill_output1:
            for bdata in kill_output2:
                if str(adata) in compatibility_dict and str(bdata) in compatibility_dict[str(adata)]:
                    continue
                pairs_temp = create_pairs(adata, bdata)
                if pairs_temp:
                    if str(adata) not in compatibility_dict:
                        compatibility_dict[str(adata)] = []
                    compatibility_dict[str(adata)].append(str(bdata))
                    if str(adata) in temp_kill_indices:
                        assert temp_kill_indices[str(adata)] == pairs_temp # check if pairs are the same
                    temp_kill_indices[str(adata)] = pairs_temp
        pairs, kill_indices = generate_sublist_pairs(kill_output1.astype(str), kill_output2.astype(str), compatibility_dict, temp_kill_indices)
        if pairs:
            kill_outputs = pd.DataFrame()
            for idx2, pair in enumerate(pairs):
                temp_kill_output1 = kill_output1.copy()
                temp_kill_output2 = kill_output2.copy()
                temp_kill_output1[:] = kill_indices[idx2]
                temp_kill_output2[:] = kill_indices[idx2]
                temp_kill_output1[pair[0]] = None
                temp_kill_output2[pair[1]] = None
                # number of positions in index depending on number of pairs
                temp_kill_output1.name = temp_name + f"_r1_kill_{idx2:0{len(str(len(pairs)))}}"
                temp_kill_output2.name = temp_name + f"_r2_kill_{idx2:0{len(str(len(pairs)))}}"
                name_dict[f"{reaction_idx}_{pattern_idx}"][f"{sub_pattern_idx}"].append(f"kill_{idx2:0{len(str(len(pairs)))}}")
                kill_outputs = pd.concat([kill_outputs, pd.concat([temp_kill_output1, temp_kill_output2], axis=1)], axis=1)
                logging.info(f"added {f"kill_{idx2:0{len(str(len(pairs)))}}"}")
            logging.info(f"added {len(pairs)} kill statement(s)")
        else:
            kill_output1 = kill_output1.apply(lambda x: list(x.keys())[0] if x[list(x.keys())[0]] else None)
            kill_output2 = kill_output2.apply(lambda x: list(x.keys())[0] if x[list(x.keys())[0]] else None)
            kill_output1.name = temp_name + "_r1_kill_0"
            kill_output2.name = temp_name + "_r2_kill_0"
            name_dict[f"{reaction_idx}_{pattern_idx}"][f"{sub_pattern_idx}"].append("kill_0")
            kill_outputs = pd.concat([kill_output1, kill_output2], axis=1)
            logging.info("added kill_0")
        kill_df = pd.concat([kill_df, kill_outputs], axis=1)
    df = pd.concat([df, kill_df] , axis=1)
    # order columns alphabetically but without the first two  columns ["smiles", "name"]
    columns = list(df.columns)
    columns.remove("smiles")
    columns.remove("name")
    columns.sort()
    df = df[["smiles", "name"] + columns]
    df.to_csv(reaction_outpath / f"results.csv")
    with open(reaction_outpath / f"kill_names.json", "w") as f:
        json.dump(name_dict, f, indent=4)
    get_stats(reaction_outpath, apply_kill=True)

def validate_reactants(bbs_path, outpath, indices=None, apply_kill=False, reaction_path=None, kill_path=None, disallow_exocyclic_doublebonds=False, uniquify=True, use_unique_matches=False):
    start = datetime.datetime.now()
    bbs_path = Path(bbs_path)
    if reaction_path:
        reaction_pattern = load_pattern_json(filename=reaction_path)
        reactions, unique_matches = load_reactions(filename=reaction_path)
    else:
        reaction_pattern = load_pattern_json()
        reactions, unique_matches = load_reactions()
    if kill_path:
        kill_patterns = load_kill_patterns(filename=kill_path)
    else:
        kill_patterns = load_kill_patterns()

    if not indices:
        indices = reaction_pattern.keys()

    if not use_unique_matches:
        unique_matches = None

    for reaction_index in indices:
        logging.info(f"Start reactants validation for {reaction_pattern[reaction_index]['name']}")
        rdkit_pattern_matching(bbs_path, outpath, reaction_index, reactions, reaction_pattern, disallow_exocyclic_doublebonds, unique_matches=unique_matches, uniquify=uniquify)
        if apply_kill:
            kill_statements(reaction_index, outpath, kill_patterns, reactions, reaction_pattern, disallow_exocyclic_doublebonds)
        else:
            reaction_outpath = get_filtered_bbs_path(outpath, reaction_index, reaction_pattern)
            pattern_smarts = get_pattern_smarts(reactions, reaction_index)
            name_dict = {}
            for pattern_name in pattern_smarts.keys():
                reaction_idx, pattern_idx, sub_pattern_idx, reactant = pattern_name.split("_")
                if f"{reaction_idx}_{pattern_idx}" not in name_dict:
                    name_dict[f"{reaction_idx}_{pattern_idx}"] = {}
                name_dict[f"{reaction_idx}_{pattern_idx}"][f"{sub_pattern_idx}"] = ["no_kill"]
            with open(reaction_outpath / f"kill_names.json", "w") as f:
                json.dump(name_dict, f, indent=4)
            get_stats(reaction_outpath, apply_kill=False)
    total_time = (datetime.datetime.now() - start)
    logging.info(f"Total time: {total_time}")


def get_stats(reaction_outpath, apply_kill=False):
    df = pd.read_csv(reaction_outpath / f"results.csv", index_col=0, low_memory=False)
    columns = list(set(df.columns) - set(["smiles", "name"]))
    info_add = " and number of killed matches" if apply_kill else ""
    logging.info("number of 1 or more matches" + info_add + ":")
    stats_df = ((df[[c for c in columns if "kill" not in c]] >= 1.0)).sum(axis=0)
    if apply_kill:
        kill_stats = ((~df[[c for c in columns if "kill" in c]].isna())).sum(axis=0)
        stats_df = pd.concat([stats_df, kill_stats], axis=0)
    stats_df = stats_df.sort_index()
    stats_df.to_csv(reaction_outpath / f"stats.csv", header=False)
    lines = str(stats_df).split("\n")[:-1]
    for line in lines:
        logging.info(line)
    logging.info("#"*50)

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-b","--bbs_path", type=str, help="path to the reactants file")
    parser.add_argument("-n","--name", type=str, default="SAVI-Space", help="name of the space and so the folder created in the output. Default: SAVI-Space")
    parser.add_argument("-o","--outpath", type=str, default=None, help="path where create the output. Default: current directory")
    parser.add_argument("-i","--indices", nargs="+", type=str, default=None, help="list of rxn numbers to include")
    parser.add_argument("--apply_kill", action="store_true", help="apply kill function to rdkit output")
    parser.add_argument("--cluster", action="store_true", help="disable tqdm for runs on cluster")
    parser.add_argument("--reaction_path", type=str, default=None, help="path to the reactions file")
    parser.add_argument("--kill_path", type=str, default=None, help="path to the kill file")
    parser.add_argument("--disallow_exocyclic_doublebonds", action="store_true", help="disallow exocyclic doublebonds for aromatic systems")
    parser.add_argument("--uniquify", action="store_true", help="uniquify rdkit matches")
    parser.add_argument("--unique_matches", action="store_true", help="use unique matches")
    parser.add_argument("--overwrite", action="store_true", help="overwrite output directory")
    args = parser.parse_args()

    initialize_logger(outpath / "reactants_validation.log")

    if args.apply_kill:
        from kill_statements import kill, create_pairs

    if args.cluster:
        from utils import tqdm

    if not args.outpath:
        outpath = Path.cwd() / args.name
    else:
        outpath = Path(args.outpath) / args.name
    if outpath.exists() and not args.overwrite and not args.indices:
        logging.error(f"{outpath} already exists. Use --overwrite to overwrite the directory.")
        raise FileExistsError(f"{outpath} already exists. Use --overwrite to overwrite the directory.")
    elif outpath.exists() and args.overwrite:
        import shutil
        shutil.rmtree(outpath)

    outpath.mkdir(parents=True, exist_ok=True)

    logging.info(f"Start reactants validation for {args.name}")

    validate_reactants(args.bbs_path,
                               outpath,
                               indices=args.indices,
                               apply_kill=args.apply_kill,
                               reaction_path=args.reaction_path,
                               kill_path=args.kill_path,
                               disallow_exocyclic_doublebonds=args.disallow_exocyclic_doublebonds,
                                uniquify=args.uniquify,
                                use_unique_matches=args.unique_matches)
    
if __name__=="__main__":
    main()