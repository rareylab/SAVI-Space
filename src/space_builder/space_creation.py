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
import logging

from utils import load_reactions, load_pattern_json, single_matching, get_valid, get_filtered_bbs_path, initialize_logger

def setup_naomi(pattern, sub_reaction_name, kill_name, reactant_paths, reaction_path, reaction_idx):
    
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
    # create spacelight json file
    json_path = reaction_path / str(reaction_idx)
    json_path.mkdir(parents=True, exist_ok=True)
    with open(json_path / f"rxn{sub_reaction_name}_{kill_name}.json", 'w') as json_file:
        json.dump(spacelight_json, json_file, indent=4)


def setup_colibri(pattern, sub_reaction_name, kill_name, reactant_paths, reaction_path, list_path):
    with open(reaction_path / f"rxn{sub_reaction_name}_{kill_name}.smirks", 'w') as smirks_file:
        smirks_file.write(pattern)

    with open(list_path, 'a') as list_file:
        list_file.write(f"{str(reaction_path / f'rxn{sub_reaction_name}_{kill_name}.smirks')} ; " +
                        f"{str(reactant_paths[0])} ; {str(reactant_paths[1])}\n")



def initialize_data_for_reactant(mode, path, name, reaction_idx, reaction, reactions, reaction_path=None, apply_kills=True, singlespace=False):
    if reaction_path:
        patterns_json = load_pattern_json(reaction_path)
    else:
        patterns_json = load_pattern_json()
    bbs_path = get_filtered_bbs_path(path, reaction_idx, patterns_json)
    output_path = path / "space_tool_input"
    output_path.mkdir(parents=True, exist_ok=True)
    if singlespace:
        list_path = output_path / f"{mode}/{name}_singlespace_{reaction_idx}.list"
    else:
        list_path = output_path / f"{mode}/{name}.list"
    bbs = pd.read_csv(bbs_path / f"results.csv", index_col=0, low_memory=False)
    with open(bbs_path / f"kill_names.json", 'r') as json_file:
        kill_names = json.load(json_file)
    bbs["index"] = bbs.index
    reaction_path = output_path / f"{mode}/reactions"
    reaction_path.mkdir(parents=True, exist_ok=True)
    molecule_path = output_path / "molecules"
    molecule_path.mkdir(parents=True, exist_ok=True)
    for sub_reaction_idx, sub in enumerate(kill_names[f"{reaction_idx}_{reaction}"]):
        sub_reaction_name = f"{reaction_idx}_{reaction}_{sub}"
        logging.info(f"Process sub reaction {sub_reaction_name}")
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

            if len(smiles) != 2 or len(smiles[0]) == 0 or len(smiles[1]) == 0:
                logging.warning(f"Skipping reaction {kill_name} due to missing reactants")
                continue
            
            reactant_paths[-1].parent.mkdir(parents=True, exist_ok=True)
            for i, s in enumerate(smiles):
                if "name" in s.columns:
                    s.to_csv(reactant_paths[i], header=False, index=False, columns=["smiles", "name"], sep=" ")
                else:
                    s.to_csv(reactant_paths[i], header=False, index=False, columns=["smiles", "index"], sep=" ")

            pattern = reactions[sub_reaction_idx]
            if mode == 'colibri':
                setup_colibri(pattern, sub_reaction_name, kill_name, reactant_paths, reaction_path, list_path)
            elif mode == 'naomi':
                setup_naomi(pattern, sub_reaction_name, kill_name, reactant_paths, reaction_path, reaction_idx)


def initialize_data(mode, path, name, reaction_index, reaction_path=None, apply_kills=True, singlespace=False):
    logging.info(f"Filter BBs for reaction {reaction_index}")
    if reaction_path:
        reactions, unique_matches  = load_reactions(reaction_path)
    else:
        reactions, unique_matches  = load_reactions()
    for reaction in range(len(reactions[reaction_index])):
        initialize_data_for_reactant(mode, path, name, reaction_index, reaction, reactions[reaction_index][reaction], reaction_path=reaction_path, apply_kills=apply_kills, singlespace=singlespace)


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
    
def create_Naomi_commands(path, args, name="SAVI-Space", indices=None):
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

    if indices:
        for idx in indices:
            bash_script.append(create_spacelight_db_command(path, name=f"{name}_{idx}", use_cactus_aromaticity=args.use_cactus_aromaticity, n_cores=args.n_cores, deprotection_path=args.deprotection_path, single_space=idx))
    else:
        bash_script.append(create_spacelight_db_command(path, name=name, use_cactus_aromaticity=args.use_cactus_aromaticity, n_cores=args.n_cores, deprotection_path=args.deprotection_path))
    bash_script.append("")
    bash_script.append("echo 'Finished spacelight'")

    with open(path / "create_naomi_space.sh", 'w') as script_file:
        script_file.write("\n".join(bash_script))

    print("__________________________________________________________")
    print(f"For creating space, set the paths in the script '{(path / 'create_naomi_space.sh').absolute()}' and run the following command:")
    print("")
    print(f"bash {(path / 'create_naomi_space.sh').absolute()}")
    print("__________________________________________________________")


def create_colibri_script(name, output_path, list_path, vendor_image_path):
    creation_script = []
    creation_script.append(f"python $SCRIPT_GENERATE_SPACELIGHT_DB_FILE -t $TOOL_SPACELIGHT_DB_CREATOR " +
            f"-i {list_path} -o {output_path} -f {name} " +
            "--tool_executable_reaction_synthesizer $TOOL_EXECUTABLE_REACTION_SYNTHESIZER " +
            f"-s $FUNCTIONAL_GROUP_LABEL -p $PROTECTING_GROUPS --logfile {output_path / 'SPACELIGHT_DB_CREATOR.log'}")
    # creation_script.append("")
    # creation_script.append(f"python $SCRIPT_GENERATE_SINGLE_FRAGMENT_SPACES -t $TOOL_EXECUTABLE_REACTION_SYNTHESIZER " +
    #         f"-i {list_path} -o {output_path / 'singleFragSpaces'} -s $FUNCTIONAL_GROUP_LABEL -p $PROTECTING_GROUPS --logfile {output_path / 'SINGLE_FRAGMENT_SPACES.log'}")
    # creation_script.append("")
    # creation_script.append(f"python $SCRIPT_MERGE_SINGLE_FRAGMENT_SPACES -t $TOOL_FRAGSPACE_MERGER " +
    #         f"-i {list_path} -d {output_path / 'singleFragSpaces'} -o {output_path / 'mergedFragSpace'} " +
    #         f"-f {name} -a {output_path / f'{name}.tfsdb'} --space_name {name} --vendor_name 'AMD University of Hamburg' " +
    #         f"--release_date {datetime.datetime.now().strftime('%Y-%m-%d')} --vendor_type academic --vendor_webpage 'https://software.zbh.uni-hamburg.de/' " +
    #         f"--vendor_image_file {vendor_image_path} --logfile {output_path / 'MERGE_SINGLE_FRAGMENT_SPACES.log'}")
    
    return creation_script


def create_colibri_commands(path, name="SAVI-Space", indices=None):
    output_path = path / "spaces/colibri"
    output_path.mkdir(parents=True, exist_ok=True)
    vendor_image_path = Path(__file__).parent.parent.parent / "images/SAVISpace.png"
    bashscripts_exports = ["#!/bin/bash", "", "echo 'Creating Space...'", ""]
    bashscripts_exports.append('if [ -z "$SCRIPT_GENERATE_SPACELIGHT_DB_FILE" ]; then')
    bashscripts_exports.append(f"    export SCRIPT_GENERATE_SPACELIGHT_DB_FILE=path-to-generate_spacelight_db_file.py")
    bashscripts_exports.append(f"    export SCRIPT_GENERATE_SINGLE_FRAGMENT_SPACES=path-to-generate_single_fragment_spaces.py")
    bashscripts_exports.append(f"    export SCRIPT_MERGE_SINGLE_FRAGMENT_SPACES=path-to-merge_single_fragment_spaces.py")
    bashscripts_exports.append("")
    bashscripts_exports.append(f"    export TOOL_SPACELIGHT_DB_CREATOR=path-to-spacelight_db_creator")
    bashscripts_exports.append(f"    export TOOL_EXECUTABLE_REACTION_SYNTHESIZER=path-to-reaction_synthesizer")
    bashscripts_exports.append(f"    export TOOL_FRAGSPACE_MERGER=path-to-fragspace_merger")
    bashscripts_exports.append("")
    bashscripts_exports.append(f"    export FUNCTIONAL_GROUP_LABEL=path-to-FunctionalGroupLabel.txt")
    bashscripts_exports.append(f"    export PROTECTING_GROUPS=path-to-ProtectingGroups.txt")
    bashscripts_exports.append("fi")
    bashscripts_exports.append("")

    if indices:
        bashscript_path = path / "scripts"
        bashscript_path.mkdir(parents=True, exist_ok=True)
        bashscripts = []
        for idx in indices:
            name = f"{name}_singlespace_{idx}"
            list_path = path / f"space_tool_input/colibri/{name}.list"
            creation_script = create_colibri_script(name, output_path, list_path, vendor_image_path)
    
            with open(bashscript_path / f"create_colibri_singespace_{idx}.sh", 'w') as script_file:
                script_file.write("\n".join(creation_script))
            bashscripts.append(str((bashscript_path / f"create_colibri_singespace_{idx}.sh").absolute()))
        with open(path / f"create_colibri_singespace.sh", 'w') as script_file:
            script_file.write("\n".join(bashscripts_exports))
            script_file.write("\n".join(bashscripts))
    else:
        list_path = path / f"space_tool_input/colibri/{name}.list"
        creation_script = create_colibri_script(f"{name}", output_path, list_path, vendor_image_path)
    
        with open(path / "create_colibri_space.sh", 'w') as script_file:
            script_file.write("\n".join(bashscripts_exports))
            script_file.write("\n".join(creation_script))
        
        with open(list_path, 'r') as list_file:
            lines = list_file.readlines()
        if len(lines) > 20:
            # put merge commands in between
            number_of_merges = 4
            loc = 0
            while loc < len(lines):
                lines.insert(loc, "@merge\n")
                loc += (len(lines)//number_of_merges + 1)
            with open(list_path, 'w') as list_file:
                list_file.write("".join(lines))
            
        print("__________________________________________________________")
        print(f"For creating space, set the paths in the script '{(path / "create_colibri_space.sh").absolute()}' and run the following command:")
        print("")
        print(f"bash {(path / 'create_colibri_space.sh').absolute()}")
        print("__________________________________________________________")


def main():
    import argparse

    parser = argparse.ArgumentParser(description='Run spacelight on a reaction')
    parser.add_argument('space_type', choices=['colibri', 'naomi'], help='Type of space to generate (colibri or naomi)')
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
        
    initialize_logger(path / "space_creation.log")

    if not (path / 'processed_bbs').exists():
        logging.error(f"Path {path} or 'processed_bbs' does not exist. Please run 'reactants_validation.py' first.")
        raise Exception(f"Path {path} or 'processed_bbs' does not exist. Please run 'reactants_validation.py' first.")

    if (path / f'space_tool_input/{args.space_type}').exists() and not args.overwrite:
        logging.error(f"Path {path / 'space_tool_input/{args.space_type}'} already exists. Please remove it or set the '--overwrite' flag.")
        raise Exception(f"Path {path / 'space_tool_input/{args.space_type}'} already exists. Please remove it or set the '--overwrite' flag.")
    elif (path / 'space_tool_input/{args.space_type}').exists():
        import shutil
        shutil.rmtree(path / 'space_tool_input/{args.space_type}')

    if (path / 'spaces/{args.space_type}').exists() and not args.overwrite:
        logging.info(f"Path {path / 'spaces/{args.space_type}'} already exists. Please remove it or set the '--overwrite' flag.")
        raise Exception(f"Path {path / 'spaces/{args.space_type}'} already exists. Please remove it or set the '--overwrite' flag.")
    elif (path / 'spaces/{args.space_type}').exists():
        import shutil
        shutil.rmtree(path / 'spaces/{args.space_type}')

    logging.info(f"Creating {'Single Space(s)' if args.indices else 'Space'}\n"+
            f"Path: {path}\n"+
            ("Using cactus aromaticity\n" if args.use_cactus_aromaticity else "\n"))

    if args.reaction_path:
        patterns_json = load_pattern_json(args.reaction_path)
    else:
        patterns_json = load_pattern_json()

    reaction_indices = list(patterns_json.keys())
    if args.indices:
        for idx in args.indices:
            logging.info(f"Reaction index: {idx}")
            if idx not in reaction_indices:
                logging.info(f"Reaction index {idx} not found")
                raise Exception(f"Reaction index {idx} not found")
        reaction_indices = args.indices

    for reaction_idx in reaction_indices:
        initialize_data(args.space_type, path, args.name, reaction_idx, reaction_path=args.reaction_path, apply_kills=not args.no_kills, singlespace=bool(args.indices))

    if args.space_type == 'colibri':        
        create_colibri_commands(path, args.name, args.indices)
    elif args.space_type == 'naomi':
        create_Naomi_commands(path, args, args.name, args.indices)
    
    
if __name__ == "__main__":
    main()