import pandas as pd
import json
from pathlib import Path
import itertools
from tqdm import tqdm
import argparse
from space_builder.utils import killing, single_matching

def count(space_path, output_path):

    output_path = Path(output_path)

    reactants_paths = list(Path(f"{space_path}/processed_bbs").glob("*"))

    products = {}
    products_wo_kill = {}
    reactants1 = {}
    reactants2 = {}
    reactants1_wo_kill = {}
    reactants2_wo_kill = {}

    pbar = tqdm(reactants_paths)
    for reactant_path in pbar:
        reaction_index = reactant_path.name.split("_")[0]
        print(f"\tCounting reactants and products in {reaction_index}")
        reactants_path = [reactant_path for reactant_path in reactants_paths if reaction_index in str(reactant_path)][0]
        df = pd.read_csv(reactants_path / "results.csv", low_memory=False)
        with open(reactants_path / "kill_names.json") as json_file:
            kill_names = json.load(json_file)
        products[reaction_index] = set()
        products_wo_kill[reaction_index] = set()
        reactants1[reaction_index] = set()
        reactants2[reaction_index] = set()
        reactants1_wo_kill[reaction_index] = set()
        reactants2_wo_kill[reaction_index] = set()
        for sub_reaction_name in kill_names:
            for subsub_name in tqdm(kill_names[sub_reaction_name], desc=sub_reaction_name, leave=False):
                for kill_name in tqdm(kill_names[sub_reaction_name][subsub_name], desc=sub_reaction_name, leave=False):
                    match_reactant1 = single_matching(df, 1, sub_reaction_name + "_"+subsub_name.split("_")[0], kill_name)
                    reactants1_wo_kill[reaction_index].update(df[match_reactant1].index)
                    match_reactant2 = single_matching(df, 2, sub_reaction_name + "_" + subsub_name.split("_")[0], kill_name)
                    reactants2_wo_kill[reaction_index].update(df[match_reactant2].index)
                    for a, b in itertools.product(df[match_reactant1].index, df[match_reactant2].index):
                        products_wo_kill[reaction_index].add(int(a) * len(df) + int(b))
                    killed_reactant1 = killing(df, 1, sub_reaction_name + "_"+subsub_name.split("_")[0], kill_name)
                    killed_reactant2 = killing(df, 2, sub_reaction_name + "_"+subsub_name.split("_")[0], kill_name)
                    valid_reactant1 = ~killed_reactant1 & match_reactant1
                    reactants1[reaction_index].update(df[valid_reactant1].index)
                    valid_reactant2 = ~killed_reactant2 & match_reactant2
                    reactants2[reaction_index].update(df[valid_reactant2].index)
                    for a, b in itertools.product(df[valid_reactant1].index, df[valid_reactant2].index):
                        products[reaction_index].add(int(a) * len(df) + int(b))

        products[reaction_index] = len(products[reaction_index])
        products_wo_kill[reaction_index] = len(products_wo_kill[reaction_index])
        reactants1[reaction_index] = len(reactants1[reaction_index])
        reactants2[reaction_index] = len(reactants2[reaction_index])
        reactants1_wo_kill[reaction_index] = len(reactants1_wo_kill[reaction_index])
        reactants2_wo_kill[reaction_index] = len(reactants2_wo_kill[reaction_index])

        pd.concat([pd.DataFrame(products, index=["products"]).T,
                pd.DataFrame(products_wo_kill, index=["not_killed_products"]).T,
                pd.DataFrame(reactants1, index=["reactants1"]).T,
                pd.DataFrame(reactants1_wo_kill, index=["not_killed_reactants1"]).T,
                pd.DataFrame(reactants2, index=["reactants2"]).T,
                pd.DataFrame(reactants2_wo_kill, index=["not_killed_reactants2"]).T], axis=1).to_csv(output_path / f"count_SAVISpace_{space}.csv")



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--space_paths", type=str, nargs="+", help="path of the space to count reactants in")
    parser.add_argument("--output_path", type=str, help="path to save the results")
    parser.add_argument("--progress", type=bool, default=False, help="Show progress bar")
    args = parser.parse_args()

    if not args.progress:
        tqdm = lambda x, *args, **kwargs: x

    for space in args.space_paths:
        print(f"Counting reactants and products in {space}")
        count(space, args.output_path)