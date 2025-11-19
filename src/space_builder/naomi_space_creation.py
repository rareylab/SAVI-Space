# This material is part of the SAVI-Space created by Malte Korn and Matthias Rarey
# at the Center for Bioinformatics, University of Hamburg, with support from
# Marc Nicklaus (NIH,NCI), Phil Judson, Raphael Klein (BioSolveIT GmbH) and
# Christian Lemmen (BioSolveIT GmbH).
# Philip Judson and Marc Nicklaus provided the LHASA transform rules and assisted with their knowledge about SAVI-2020.
# SAVI-Space and all its components including this file/directory is licensed under CC-BY 4.0.

from pathlib import Path
import pandas as pd
import json
import datetime

from utils import load_reactions, load_pattern_json, single_matching, get_valid, get_filtered_bbs_path


def get_bbs_for_reactant(path, reaction_idx, reaction, reactions, apply_kills=True):
    patterns_json = load_pattern_json()
    bbs_path = get_filtered_bbs_path(path, reaction_idx, patterns_json)
    output_path = path / "space_tool_input"
    output_path.mkdir(parents=True, exist_ok=True)
    bbs = pd.read_csv(bbs_path / f"results.csv", index_col=0, low_memory=False)
    with open(bbs_path / f"kill_names.json", 'r') as json_file:
        kill_names = json.load(json_file)
    bbs["index"] = bbs.index
    reaction_path = output_path / "Naomi/reactions"
    reaction_path.mkdir(parents=True, exist_ok=True)
    molecule_path = output_path / "molecules"
    molecule_path.mkdir(parents=True, exist_ok=True)
    processed = [pd.DataFrame(columns=bbs.columns), pd.DataFrame(columns=bbs.columns)]
    for sub_reaction_idx, sub in enumerate(kill_names[f"{reaction_idx}_{reaction}"]):
        sub_reaction_name = f"{reaction_idx}_{reaction}_{sub}"
        print(f"Starting sub reaction {sub_reaction_name}")
        for kill_name in kill_names[f"{reaction_idx}_{reaction}"][sub]:
            reactant_paths = []
            smiles = []
            for reactant_idx in range(1, 3):
                if apply_kills:
                    smiles.append(get_valid(bbs, reactant_idx, sub_reaction_name, kill_name))
                    if len(smiles[-1]) == 0:
                        continue
                    reactant_paths.append(molecule_path / f"{reaction_idx}/rxn{sub_reaction_name}_{reactant_idx}_{kill_name}.smi")
                else:
                    kill_name = "no_kill"
                    smiles.append(bbs[single_matching(bbs, reactant_idx, sub_reaction_name)])
                    if len(smiles[-1]) == 0:
                        continue
                    reactant_paths.append(molecule_path / f"{reaction_idx}/rxn{sub_reaction_name}_{reactant_idx}.smi")

            if set(smiles[0].index).issubset(set(processed[0].index)):
                len_before = len(smiles[1])
                smiles[1] = smiles[1].drop(processed[1].index, errors='ignore')
                print(f"Skipping {len_before-len(smiles[1])} of {len_before} reactants for reaction {sub_reaction_name}_{kill_name} due to already processed reactant pairs")
            if set(smiles[1].index).issubset(set(processed[1])):
                len_before = len(smiles[0])
                smiles[0] = smiles[0].drop(processed[0].index, errors='ignore')
                print(f"Skipping {len_before-len(smiles[0])} of {len_before} reactants for reaction {sub_reaction_name}_{kill_name} due to already processed reactant pairs")

            if len(smiles) != 2 or len(smiles[0]) == 0 or len(smiles[1]) == 0:
                print(f"Skipping reaction {kill_name} due to missing reactants")
                continue

            reactant_paths[-1].parent.mkdir(parents=True, exist_ok=True)
            for i, s in enumerate(smiles):
                if "name" in s.columns:
                    s.to_csv(reactant_paths[i], header=False, index=False, columns=["smiles", "name"], sep=" ")
                else:
                    s.to_csv(reactant_paths[i], header=False, index=False, columns=["smiles", "index"], sep=" ")

            # create spacelight json file
            json_path = reaction_path / str(reaction_idx)
            json_path.mkdir(parents=True, exist_ok=True)
            try:
                pattern = reactions[sub_reaction_idx]
            except:
                print(reaction_idx, sub_reaction_idx)
                print(reactions)
                exit()
            smirks = [{"components": [1,2], "reaction": pattern}]

            spacelight_json = {
                "topologies" : [
                    {
                        "name": f"rxn{sub_reaction_name}_{kill_name}",
                        "reactions": smirks,
                        "reagentGroups": [
                            {
                                "groupId": 1,
                                "reagents": str(reactant_paths[0])
                            },
                            {
                                "groupId": 2,
                                "reagents": str(reactant_paths[1])
                            }
                        ]
                    }
                ]
            }
            with open(json_path / f"rxn{sub_reaction_name}_{kill_name}.json", 'w') as json_file:
                json.dump(spacelight_json, json_file, indent=4)

            processed[0] = pd.concat([processed[0], smiles[0]])
            processed[1] = pd.concat([processed[1], smiles[1]])


def create_bbs(path, reaction_index, reaction_path=None, apply_kills=True):
    print(f"Creating BBs for reaction {reaction_index}")
    if reaction_path:
        reactions, unique_matches  = load_reactions(reaction_path)
    else:
        reactions, unique_matches  = load_reactions()
    print(reactions)
    for reaction in range(len(reactions[reaction_index])):
        get_bbs_for_reactant(path, reaction_index, reaction, reactions[reaction_index][reaction], apply_kills=apply_kills)


def create_spacelight_db_command(path,
                           name="SAVI-Space",
                           n_cores="4",
                           use_cactus_aromaticity=False,
                           deprotection_path=None,
                           single_space=None):
    output_path = path / "spaces/Naomi"
    output_path.mkdir(parents=True, exist_ok=True)
    input_path = path / "space_tool_input/Naomi/reactions"
    if single_space:
        input_path = input_path / str(single_space)

    return (f"$SPACELIGHT generate -i {input_path.absolute()} -f {(output_path / f'{name}.tfsdb').absolute()} " +
    f"-c {1 if use_cactus_aromaticity else 0} {f'-d {deprotection_path}' if deprotection_path else ''} -p {n_cores}")

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Run spacelight on a reaction')
    parser.add_argument("-n","--name", type=str, default="SAVI-Space", help="name of the space. Default: SAVI-Space")
    parser.add_argument("-o","--outpath", type=str, default=None, help="path where the space is located. Default: current directory")
    parser.add_argument("-i","--indices", nargs="+", type=str, default=None, help="reaction indices to run")
    parser.add_argument("-c", "--use_cactus_aromaticity", action="store_true", help="use cactus aromaticity")
    parser.add_argument("--deprotection_path", type=str, default=None, help="path to the deprotection file")
    parser.add_argument("--reaction_path", type=str, default=None, help="path to the reactions file")
    parser.add_argument("--n_cores", type=int, default=4, help="number of cores to use")
    parser.add_argument("--no_kills", action="store_true", help="skip kills")
    parser.add_argument("--overwrite", action="store_true", help="overwrite existing files")
    args = parser.parse_args()

    if args.outpath:
        path = Path(args.outpath) / "args.name"
    else:
        path = Path.cwd() / args.name

    if not (path / 'processed_bbs').exists():
        raise Exception(f"Path {path} or 'processed_bbs' does not exist. Please run 'reactants_validation.py' first.")

    if (path / 'space_tool_input/Naomi').exists() and not args.overwrite:
        raise Exception(f"Path {path / 'space_tool_input/Naomi'} already exists. Please remove it or set the '--overwrite' flag.")
    elif (path / 'space_tool_input/Naomi').exists():
        import shutil
        shutil.rmtree(path / 'space_tool_input/Naomi')

    if (path / 'spaces/Naomi').exists() and not args.overwrite:
        raise Exception(f"Path {path / 'spaces/Naomi'} already exists. Please remove it or set the '--overwrite' flag.")
    elif (path / 'spaces/Naomi').exists():
        import shutil
        shutil.rmtree(path / 'spaces/Naomi')

    print(f"Creating {'Single Space(s)' if args.indices else 'Space'}\n"+
            f"Path: {path}\n"+
            ("Using cactus aromaticity\n" if args.use_cactus_aromaticity else ""))

    if args.reaction_path:
        patterns_json = load_pattern_json(args.reaction_path)
    else:
        patterns_json = load_pattern_json()

    reaction_indices = list(patterns_json.keys())
    if args.indices:
        for idx in args.indices:
            print(f"Reaction index: {idx}")
            if idx not in reaction_indices:
                raise Exception(f"Reaction index {idx} not found")
        reaction_indices = args.indices

    for reaction_idx in reaction_indices:
        create_bbs(path, reaction_idx, reaction_path=args.reaction_path, apply_kills=not args.no_kills)


    bash_script = [
        "#!/bin/bash",
        "",
        "echo 'Running spacelight...'",
        "",
        "if [ -z \"$SPACELIGHT\" ]; then",
        "    echo 'Please set spacelight binary path as SPACELIGHT environment variable'",
        "    echo 'export SPACELIGHT=/path/to/spacelight'",
        "    exit 1",
        "fi",
        ""]

    if args.indices:
        for idx in args.indices:
            bash_script.append(create_spacelight_db_command(path, name=f"{args.name}_{idx}", use_cactus_aromaticity=args.use_cactus_aromaticity, n_cores=args.n_cores, deprotection_path=args.deprotection_path, single_space=idx))
    else:
        bash_script.append(create_spacelight_db_command(path, name=args.name, use_cactus_aromaticity=args.use_cactus_aromaticity, n_cores=args.n_cores, deprotection_path=args.deprotection_path))
    bash_script.append("")
    bash_script.append("echo 'Finished spacelight'")

    with open(path / "create_naomi_space.sh", 'w') as script_file:
        script_file.write("\n".join(bash_script))

    print("__________________________________________________________")
    print(f"For creating space, set the paths in the script '{(path / 'create_naomi_space.sh').absolute()}' and run the following command:")
    print("")
    print(f"bash {(path / 'create_naomi_space.sh').absolute()}")
    print("__________________________________________________________")