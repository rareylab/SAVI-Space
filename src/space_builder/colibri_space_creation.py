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


def get_bbs_for_reactant(path, name, reaction_idx, reaction, reactions, apply_kills=True, overwrite=False):
    patterns_json = load_pattern_json()
    bbs_path = get_filtered_bbs_path(path, reaction_idx, patterns_json)
    output_path = path / "space_tool_input"
    output_path.mkdir(parents=True, exist_ok=True)
    bbs = pd.read_csv(bbs_path / f"results.csv", index_col=0, low_memory=False)
    with open(bbs_path / f"kill_names.json", 'r') as json_file:
        kill_names = json.load(json_file)
    bbs["index"] = bbs.index
    reaction_path = output_path / "Colibri/reactions"
    reaction_path.mkdir(parents=True, exist_ok=True)
    molecule_path = output_path / "molecules"
    molecule_path.mkdir(parents=True, exist_ok=True)
    list_path = output_path / f"Colibri/{name}.list"
    for sub_reaction_idx, sub in enumerate(kill_names[f"{reaction_idx}_{reaction}"]):
        sub_reaction_name = f"{reaction_idx}_{reaction}_{sub}"
        print(f"Starting sub reaction {sub_reaction_name}")
        for kill_name in kill_names[f"{reaction_idx}_{reaction}"][sub]:
            smiles = []
            for reactant_idx in range(1, 3):
                if apply_kills:
                    smiles.append(get_valid(bbs, reactant_idx, sub_reaction_name, kill_name))
                    sub_reaction_idx_path = lambda x: molecule_path / f"{reaction_idx}/rxn{sub_reaction_name}_{x}_{kill_name}.smi"

                else:
                    kill_name = "no_kill"
                    print("Single matching", single_matching(bbs, reactant_idx, sub_reaction_name).sum())
                    smiles.append(bbs[single_matching(bbs, reactant_idx, sub_reaction_name)])
                    sub_reaction_idx_path = lambda x: molecule_path / f"{reaction_idx}/rxn{sub_reaction_name}_{x}.smi"

            if len(smiles) != 2 or len(smiles[0]) == 0 or len(smiles[1]) == 0:
                print(f"Skipping reaction {kill_name} due to missing reactants")
                continue

            reactant_paths = [sub_reaction_idx_path(1), sub_reaction_idx_path(2)]
            reactant_paths[-1].parent.mkdir(parents=True, exist_ok=True)
            for i, s in enumerate(smiles):
                if "name" in s.columns:
                    s.to_csv(reactant_paths[i], header=False, index=False, columns=["smiles", "name"], sep=" ")
                else:
                    s.to_csv(reactant_paths[i], header=False, index=False, columns=["smiles", "index"], sep=" ")

            # create BSIT list file

            try:
                pattern = reactions[sub_reaction_idx]
            except:
                print(reaction_idx, sub_reaction_idx)
                print(reactions)
                exit()

            #save reaction as smirks
            with open(reaction_path / f"rxn{sub_reaction_name}_{kill_name}.smirks", 'w') as smirks_file:
                smirks_file.write(pattern)

            with open(list_path, 'a') as list_file:
                list_file.write(f"{str(reaction_path / f'rxn{sub_reaction_name}_{kill_name}.smirks')} ; " +
                                f"{str(reactant_paths[0])} ; {str(reactant_paths[1])}\n")


def create_bbs(path, name, reaction_index, reaction_path=None, apply_kills=True):
    print(f"Creating BBs for reaction {reaction_index}")
    if reaction_path:
        reactions, unique_matches  = load_reactions(reaction_path)
    else:
        reactions, unique_matches  = load_reactions()
    for reaction in range(len(reactions[reaction_index])):
        get_bbs_for_reactant(path, name, reaction_index, reaction, reactions[reaction_index][reaction], apply_kills=apply_kills)


def create_Colibri_commands(path, name="SAVI-Space", use_scripts=False):
    output_path = path / "spaces/Colibri"
    output_path.mkdir(parents=True, exist_ok=True)
    list_path = path / f"space_tool_input/Colibri/{name}.list"

    creation_script = ["#!/bin/bash", "", "echo 'Creating Space...'", ""]
    if use_scripts:
        creation_script.append('if [ -z "$SCRIPT_GENERATE_SPACELIGHT_DB_FILE" ]; then')
        creation_script.append(f"    export SCRIPT_GENERATE_SPACELIGHT_DB_FILE=path-to-generate_spacelight_db_file.py")
        creation_script.append(f"    export SCRIPT_GENERATE_SINGLE_FRAGMENT_SPACES=path-to-generate_single_fragment_spaces.py")
        creation_script.append(f"    export SCRIPT_MERGE_SINGLE_FRAGMENT_SPACES=path-to-merge_single_fragment_spaces.py")
        creation_script.append("")
        creation_script.append(f"    export TOOL_SPACELIGHT_DB_CREATOR=path-to-spacelight_db_creator")
        creation_script.append(f"    export TOOL_EXECUTABLE_REACTION_SYNTHESIZER=path-to-reaction_synthesizer")
        creation_script.append(f"    export TOOL_FRAGSPACE_MERGER=path-to-fragspace_merger")
        creation_script.append("")
        creation_script.append(f"    export FUNCTIONAL_GROUP_LABEL=path-to-FunctionalGroupLabel.txt")
        creation_script.append(f"    export PROTECTING_GROUPS=path-to-ProtectingGroups.txt")
        creation_script.append("fi")
        creation_script.append("")
        creation_script.append(f"python $SCRIPT_GENERATE_SPACELIGHT_DB_FILE -t $TOOL_SPACELIGHT_DB_CREATOR " +
                f"-i {list_path} -o {output_path} -f {name} " +
                "--tool_executable_reaction_synthesizer $TOOL_EXECUTABLE_REACTION_SYNTHESIZER " +
                f"-s $FUNCTIONAL_GROUP_LABEL -p $PROTECTING_GROUPS --logfile {output_path / 'SPACELIGHT_DB_CREATOR.log'}")
        creation_script.append("")
        creation_script.append(f"python $SCRIPT_GENERATE_SINGLE_FRAGMENT_SPACES -t $TOOL_EXECUTABLE_REACTION_SYNTHESIZER " +
                f"-i {list_path} -o {output_path / 'singleFragSpaces'} -s $FUNCTIONAL_GROUP_LABEL -p $PROTECTING_GROUPS --logfile {output_path / 'SINGLE_FRAGMENT_SPACES.log'}")
        creation_script.append("")
        creation_script.append(f"python $SCRIPT_MERGE_SINGLE_FRAGMENT_SPACES -t $TOOL_FRAGSPACE_MERGER " +
                f"-i {list_path} -d {output_path / 'singleFragSpaces'} -o {output_path / 'mergedFragSpace'} " +
                f"-f {name} -a {output_path / f'{name}.tfsdb'} --space_name {name} --vendor_name 'AMD University of Hamburg' " +
                f"--release_date {datetime.datetime.now().strftime('%Y-%m-%d')} --vendor_type academic --vendor_webpage 'https://software.zbh.uni-hamburg.de/' " +
                f"--vendor_image_file /work/korn/SAVI-Space/SAVISpace_icon.png --logfile {output_path / 'MERGE_SINGLE_FRAGMENT_SPACES.log'}")
    else:
        raise NotImplementedError("Please set 'use_scripts' to True")
    return creation_script

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
        path = Path(args.outpath) / args.name
    else:
        path = Path.cwd() / args.name

    if not (path / 'processed_bbs').exists():
        raise Exception(f"Path {path} or 'processed_bbs' does not exist. Please run 'reactants_validation.py' first.")

    if (path / 'space_tool_input/Colibri').exists() and not args.overwrite:
        raise Exception(f"Path {path / 'space_tool_input/Colibri'} already exists. Please remove it or set the '--overwrite' flag.")
    elif (path / 'space_tool_input/Colibri').exists():
        import shutil
        shutil.rmtree(path / 'space_tool_input/Colibri')

    if (path / 'spaces/Colibri').exists() and not args.overwrite:
        raise Exception(f"Path {path / 'spaces/Colibri'} already exists. Please remove it or set the '--overwrite' flag.")
    elif (path / 'spaces/Colibri').exists():
        import shutil
        shutil.rmtree(path / 'spaces/Colibri')

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

    list_path = path / "space_tool_input/Colibri/"
    list_path.mkdir(parents=True, exist_ok=True)
    with open(list_path / f"{args.name}.list", 'w') as list_file:
        pass

    for idx, reaction_idx in enumerate(reaction_indices):
        create_bbs(path, args.name, reaction_idx, reaction_path=args.reaction_path, apply_kills=not args.no_kills)

    with open(list_path / f"{args.name}.list", 'r') as list_file:
        lines = list_file.readlines()
    # put merge commands in between
    number_of_merges = 4
    loc = 0
    while loc < len(lines):
        lines.insert(loc, "@merge\n")
        loc += (len(lines)//number_of_merges + 1)
    with open(list_path / f"{args.name}.list", 'w') as list_file:
        list_file.write("".join(lines))


    bash_script = create_Colibri_commands(path, name=args.name, use_scripts=True)
    with open(path / "create_Colibri_space.sh", 'w') as script_file:
        script_file.write("\n".join(bash_script))
    print("__________________________________________________________")
    print(f"For creating space, set the paths in the script '{(path / "create_Colibri_space.sh").absolute()}' and run the following command:")
    print("")
    print(f"bash {(path / 'create_Colibri_space.sh').absolute()}")
    print("__________________________________________________________")