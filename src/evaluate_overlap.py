# This material is part of the SAVI-Space created by Malte Korn and Matthias Rarey
# at the Center for Bioinformatics, University of Hamburg, with support from
# Marc Nicklaus (NIH,NCI), Phil Judson, Raphael Klein (BioSolveIT GmbH) and
# Christian Lemmen (BioSolveIT GmbH).
# Philip Judson and Marc Nicklaus provided the LHASA transform rules and assisted with their knowledge about SAVI-2020.
# SAVI-Space and all its components including this file/directory is licensed under CC-BY-NC 4.0.

import pandas as pd
import json
from pathlib import Path
import itertools
from tqdm import tqdm
import argparse
from space_builder.utils import killing, single_matching
import random

def evaluate_overlap(space_path_1, space_path_2, output_path, examples):

    space_path_1 = Path(space_path_1)
    space_path_2 = Path(space_path_2)
    space_name_1 = space_path_1.name
    space_name_2 = space_path_2.name

    output_path = Path(output_path)

    reactants_paths1 = list((space_path_1 / "processed_bbs").glob("*"))
    reactants_paths2 = list((space_path_2 / "processed_bbs").glob("*"))

    output_path = output_path / f"overlap_{space_name_1}_{space_name_2}"
    output_path.mkdir(parents=True, exist_ok=True)

    result = {}

    pbar = tqdm(list(reactants_paths1))
    for reactants_path1 in pbar:
        reaction_index = reactants_path1.name.split("_")[0]
        reactants_path2 = [reactant_path for reactant_path in reactants_paths2 if reaction_index in str(reactant_path)]
        if len(reactants_path2) != 1:
            pbar.write(f"Reaction {reaction_index} not found in {space_name_2}")
            continue
        else:
            reactants_path2 = reactants_path2[0]
        df1 = pd.read_csv(reactants_path1 / "results.csv", index_col=0, low_memory=False)
        df2 = pd.read_csv(reactants_path2 / "results.csv", index_col=0, low_memory=False)
        with open(reactants_path1 / "kill_names.json") as json_file:
            kill_names1 = json.load(json_file)
        with open(reactants_path2 / "kill_names.json") as json_file:
            kill_names2 = json.load(json_file)
        assert len(df1) == len(df2)
        pairs1 = set()
        pairs2 = set()
        # space 1
        for sub_reaction_name in tqdm(kill_names1, leave=False):
            for subsub_name in tqdm(kill_names1[sub_reaction_name], desc=sub_reaction_name, leave=False):
                for kill_name in tqdm(kill_names1[sub_reaction_name][subsub_name], desc=subsub_name, leave=False):
                    print(sub_reaction_name, kill_name)
                    match_reactant1 = single_matching(df1, 1, sub_reaction_name + "_"+subsub_name.split("_")[0], kill_name)
                    match_reactant2 = single_matching(df1, 2, sub_reaction_name + "_" + subsub_name.split("_")[0], kill_name)
                    killed_reactant1 = killing(df1, 1, sub_reaction_name + "_"+subsub_name.split("_")[0], kill_name)
                    killed_reactant2 = killing(df1, 2, sub_reaction_name + "_"+subsub_name.split("_")[0], kill_name)
                    valid_reactant1 = ~killed_reactant1 & match_reactant1
                    valid_reactant2 = ~killed_reactant2 & match_reactant2
                    for a, b in itertools.product(df1[valid_reactant1].index, df1[valid_reactant2].index):
                        pairs1.add(int(a) * len(df1) + int(b))
        #space 2
        for sub_reaction_name in tqdm(kill_names2, leave=False):
            for subsub_name in tqdm(kill_names2[sub_reaction_name], desc=sub_reaction_name, leave=False):
                for kill_name in tqdm(kill_names2[sub_reaction_name][subsub_name], desc=subsub_name, leave=False):
                    match_reactant1 = single_matching(df2, 1, sub_reaction_name + "_"+subsub_name.split("_")[0], kill_name)
                    match_reactant2 = single_matching(df2, 2, sub_reaction_name + "_"+subsub_name.split("_")[0], kill_name)
                    killed_reactant1 = killing(df2, 1, sub_reaction_name + "_"+subsub_name.split("_")[0], kill_name)
                    killed_reactant2 = killing(df2, 2, sub_reaction_name + "_"+subsub_name.split("_")[0], kill_name)
                    valid_reactant1 = ~killed_reactant1 & match_reactant1
                    valid_reactant2 = ~killed_reactant2 & match_reactant2
                    for a, b in itertools.product(df2[valid_reactant1].index, df2[valid_reactant2].index):
                        pairs2.add(int(a) * len(df2) + int(b))

        intersection = len(pairs1.intersection(pairs2))
        pbar.write(f"Reaction {reaction_index} has {len(pairs1)} pairs in space 1 and {len(pairs2)} pairs in space 2")
        pbar.write(f"Reaction {reaction_index} has {intersection} pairs in common")
        result[reaction_index] = {}
        result[reaction_index]["size_" + space_name_1] = len(pairs1)
        result[reaction_index]["size_" + space_name_2] = len(pairs2)
        result[reaction_index]["overlap"] = intersection
        if examples:
            example_pairs = random.sample(sorted(pairs1 - pairs2), min(100, len(pairs1 - pairs2)))
            example_pairs = [(df1.loc[a // len(df1)].smiles, df1.loc[a % len(df1)].smiles) for a in example_pairs]
            with open(output_path / f"{reaction_index}_only_{space_name_1}.txt", "w") as file:
                file.write("\n".join([f"{a} {b}" for a, b in example_pairs]))
            example_pairs = random.sample(sorted(pairs2 - pairs1), min(100, len(pairs2 - pairs1)))
            example_pairs = [(df2.loc[a // len(df2)].smiles, df2.loc[a % len(df2)].smiles) for a in example_pairs]
            with open(output_path / f"{reaction_index}_only_{space_name_2}.txt", "w") as file:
                file.write("\n".join([f"{a} {b}" for a, b in example_pairs]))

    pd.DataFrame(result).T.to_csv(output_path / "result.csv")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--space_paths", type=str, nargs=2, required=True)
    parser.add_argument("--output_path", type=str, required=True)
    parser.add_argument("-e", "--examples", action="store_true")
    args = parser.parse_args()

    evaluate_overlap(args.space_paths[0], args.space_paths[1],
                      args.output_path, args.examples)