# This material is part of the SAVI-Space created by Malte Korn and Matthias Rarey
# at the Center for Bioinformatics, University of Hamburg, with support from
# Marc Nicklaus (NIH,NCI), Phil Judson, Raphael Klein (BioSolveIT GmbH) and
# Christian Lemmen (BioSolveIT GmbH).
# Philip Judson and Marc Nicklaus provided the LHASA transform rules and assisted with their knowledge about SAVI-2020.
# SAVI-Space and all its components including this file/directory is licensed under CC-BY 4.0.

import argparse
import pandas as pd
import numpy as np
from tqdm import tqdm
from pathlib import Path

from rdkit import Chem
from rdkit import RDLogger
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.Descriptors import ExactMolWt
from rdkit.Chem.GraphDescriptors import BertzCT

from subprocess import PIPE, run

tqdm.pandas()

_normalization_transforms = """
//	Name	SMIRKS
Sulfoxid to S=O\t[S+:1]-[O-:2]>>[S+0:1]=[O;+0:2]
Sulfoxid to S=O\t[S+:1]-[O:2]>>[S+0:1]=[O;+0:2]
Protonize Oxigens\t[O-:1]>>[OH;+0:1]
"""
_normalizer_params = rdMolStandardize.CleanupParameters()
_normalizer = rdMolStandardize.NormalizerFromData(
    _normalization_transforms, _normalizer_params
)

def standardize(mol):
    if mol is None:
        return None

    if mol.GetNumAtoms() <= 1:
        return None

    # removeHs, disconnect metal atoms, normalize the molecule, reionize the molecule
    clean_mol = rdMolStandardize.Cleanup(mol)

    # if many fragments, get the "parent" (the actual mol we are interested in)
    parent_clean_mol = rdMolStandardize.FragmentParent(clean_mol)

    uncharged_parent_clean_mol = _normalizer.normalize(parent_clean_mol)

    invalid_atoms = rdMolStandardize.DisallowedAtomsValidation([Chem.Atom(atom) for atom in ["Sn","Hg","Cu","Sm","Mg", "Pd", "Ni", "W", ]])
    if invalid_atoms.validate(uncharged_parent_clean_mol):
        return None

    #remove atom maps
    [a.SetAtomMapNum(0) for a in uncharged_parent_clean_mol.GetAtoms()]

    # # remove stereochemistry
    # Chem.RemoveStereochemistry(uncharged_parent_clean_mol)

    return uncharged_parent_clean_mol

def naomi_check(file_path, naomi_check_binary):
    naomi_check_binary
    p = run([naomi_check_binary,
        "--input", file_path.absolute(),
        "--smarts", "[*]",
        "--out", file_path.parent / "temp.csv"], stdout=PIPE, stderr=PIPE)
    df = pd.read_csv(file_path.parent / "temp.csv")
    (file_path.parent / "temp.csv").unlink()
    # check if dataframe have two or one columns and name them name and ID
    if df.shape[1] == 3:
        df[["molecule","name"]].to_csv(file_path, header=None, index=False, sep=" ")
    elif df.shape[1] == 2:
        df["molecule"].to_csv(file_path, header=None, index=False, sep=" ")
    else:
        raise ValueError("Input file must have one or two columns")

def divide_disconnected_smiles(file_path):
    with open(file_path, "r") as f:
        smiles = f.readlines()
    results = dict()
    with open(f"{file_path.stem}_fragments.smi", "w") as f:
        for s in tqdm(smiles, desc="Dividing disconnected smiles"):
            s = s.split(" ", maxsplit=1)
            for frag in s[0].split("."):
                if len(s) == 2:
                    f.write(f"{frag} {s[1]}")
                else:
                    f.write(f"{frag}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=str, required=True, help="Input smilesfile to standardize")
    parser.add_argument("--divide_disconnected_smiles", action="store_true", help="Divide disconnected smiles")
    parser.add_argument("--naomi_check_binary", type=str, help="Path to naomi_check binary")
    parser.add_argument("--max_heavy_atoms", type=int, default=None, help="Maximal number of heavy atoms")
    parser.add_argument("--min_heavy_atoms", type=int, default=None, help="Minimal number of heavy atoms")
    parser.add_argument("--max_molecular_weight", type=float, default=None, help="Maximal molecular weight")
    parser.add_argument("--min_molecular_weight", type=float, default=None, help="Minimal molecular weight")
    parser.add_argument("--max_bertz_ct", type=float, default=None, help="Maximal bertz ct")
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose output")
    args = parser.parse_args()



    process_log = []

    if not args.verbose:
        RDLogger.DisableLog('rdApp.*')

    input_file = Path(args.input)
    if input_file.suffix.lower() == ".sdf":
        sdf_entrys = 0
        sdf_valid_entrys = 0
        if not Path(f"{input_file.stem}.smi").exists():
            mols = Chem.SDMolSupplier(input_file)
            output_file = Path(f"{input_file.stem}.smi")
            try:
                name_prop = None
                with open(output_file, "w") as f:
                    for mol in tqdm(mols, desc="Converting sdf to smi"):
                        sdf_entrys += 1
                        if mol is not None:
                            sdf_valid_entrys += 1
                            if not name_prop or not name_prop in mol.GetPropNames():
                                names_prop =  [s for s in mol.GetPropNames() if any(sub in s.lower() for sub in ("id", "name"))]
                                name_prop = min(names_prop, key=len) if names_prop else None
                                if (name_prop):
                                    print(f"Found property: {name_prop}. Use as molecule name")
                            if (name_prop):
                                f.write(f"{Chem.MolToSmiles(mol)} {mol.GetProp(name_prop)}\n")
                            else:
                                f.write(f"{Chem.MolToSmiles(mol)}")
            except Exception as e:
                if output_file.is_file():
                    output_file.unlink()
                raise Exception(f"Something went wrong while converting the sd file.\nError: {e}") from e

            print(f"Converted {sdf_valid_entrys}/{sdf_entrys} valid sdf entrys to smi")
        else:
            print(f"SD file already converted. Remove `{input_file.stem}.smi` if you want to convert again.")
            exit()
        process_log.append(["sdf_entrys", sdf_entrys])
        process_log.append(["sdf_valid_entrys", sdf_valid_entrys])
    elif input_file.suffix.lower() == ".smi":
        smiles_entrys = 0
        smiles_valid_entrys = 0
        with open(input_file, "r") as f:
            for line in f:
                smiles_entrys += 1
                if Chem.MolFromSmiles(line.split(" ")[0]) is not None:
                    smiles_valid_entrys += 1
        print(f"Checked {smiles_valid_entrys}/{smiles_entrys} valid smiles entrys")
        process_log.append(["smiles_entrys", smiles_entrys])
        process_log.append(["smiles_valid_entrys", smiles_valid_entrys])
    else:
        raise ValueError("Input file must be either sdf or smi")

    if args.divide_disconnected_smiles:
        output_fragments_file = Path(f"{input_file.stem}_fragments.smi")
        if not output_fragments_file.exists():
            divide_disconnected_smiles(Path(f"{input_file.stem}.smi"))
        smiles_entrys = 0
        with open(output_fragments_file, "r") as f:
            smiles_entrys = len(f.readlines())
        print(f"Diveded {smiles_entrys} smiles")
        process_log.append(["divided_smiles", smiles_entrys])

    output_fragments_file = Path(f"{input_file.stem}{'_fragments' if args.divide_disconnected_smiles else ''}.smi")
    if args.naomi_check_binary is not None:
        naomi_check(output_fragments_file, args.naomi_check_binary)
        smiles_entrys = 0
        with open(output_fragments_file, "r") as f:
            smiles_entrys = len(f.readlines())
        print(f"{smiles_entrys} valid smiles with naomi_check")
        process_log.append(["naomi_valid_smiles", smiles_entrys])

    # load or create csv with additional informations like MW
    df_name = f"{input_file.stem}{'_fragments' if args.divide_disconnected_smiles else ''}.csv"
    if not Path(df_name).exists():
        df = pd.read_csv(output_fragments_file, header=None, sep=" ", names=["smiles", "name"])
    else:
        df = pd.read_csv(df_name)
        # check if loaded df hold the same smiles
        df_compare = pd.read_csv(output_fragments_file, header=None, sep=" ", names=["smiles", "name"])["smiles"].sort_values().reset_index(drop=True)
        if not df_compare.equals(df["smiles"].sort_values().reset_index(drop=True)):
            raise ValueError(f"Existing csv file do not hold the same molecules as the input {input_file.suffix.lower()[1:]} file")

    # load mols
    df["molecule"] = df["smiles"].progress_apply(lambda x: Chem.MolFromSmiles(x) if x is not None else np.nan)

    if (args.min_heavy_atoms is not None or args.max_heavy_atoms is not None) and "heavy_atoms" not in df.columns:
        df["heavy_atoms"] = df["molecule"].progress_apply(lambda x: x.GetNumHeavyAtoms() if x is not None else np.nan)

    if (args.min_molecular_weight is not None or args.max_molecular_weight is not None) and "MW" not in df.columns:
        df["MW"] = df["molecule"].progress_apply(lambda x: ExactMolWt(x) if x is not None else np.nan)
    #save
    df.drop(columns=["molecule"]).to_csv(df_name, index=False)

    # remove molecules with more than max_heavy_atoms or less than min_heavy_atoms
    if args.min_heavy_atoms is not None:
        name = f"{input_file.stem}_min-heavy-atoms_{args.min_heavy_atoms}"
        df = df[df["heavy_atoms"] >= args.min_heavy_atoms]
        print(f"{df.shape[0]} molecules with more than {args.min_heavy_atoms} heavy atoms")
        process_log.append(["min_heavy_atoms", df.shape[0]])

    output_name = input_file.stem
    if args.max_heavy_atoms is not None:
        output_name = f"{output_name}_max-heavy-atoms_{args.max_heavy_atoms}"
        df = df[df["heavy_atoms"] <= args.max_heavy_atoms]
        print(f"{df.shape[0]} molecules with less than {args.max_heavy_atoms} heavy atoms")
        process_log.append(["max_heavy_atoms", df.shape[0]])

    if args.min_molecular_weight is not None:
        output_name = f"{output_name}_min-molecular-weight_{args.min_molecular_weight}"
        df = df[df["MW"] >= args.min_molecular_weight]
        print(f"{df.shape[0]} molecules with more than {args.min_molecular_weight} molecular weight")
        process_log.append(["min_molecular_weight", df.shape[0]])

    if args.max_molecular_weight is not None:
        output_name = f"{output_name}_max-molecular-weight_{args.max_molecular_weight}"
        df = df[df["MW"] <= args.max_molecular_weight]
        print(f"{df.shape[0]} molecules with less than {args.max_molecular_weight} molecular weight")
        process_log.append(["max_molecular_weight", df.shape[0]])

    if args.max_bertz_ct is not None:
        output_name = f"{output_name}_max-bertz-ct_{args.max_bertz_ct}"
        df = df[df["molecule"].progress_apply(lambda x: BertzCT(x) if x is not None else np.nan) <= args.max_bertz_ct]
        print(f"{df.shape[0]} molecules with less than {args.max_bertz_ct} bertz ct")
        process_log.append(["max_bertz_ct", df.shape[0]])

    print("Standardize molecules:")
    df["molecule"] = df["molecule"].progress_apply(standardize)
    df = df.dropna(subset=["molecule"])
    print(f"{df.shape[0]} valid molecules after standardization")
    process_log.append(["standardization", df.shape[0]])
    df["smiles"] = df["molecule"].progress_apply(lambda x: Chem.MolToSmiles(x) if x is not None else np.nan)

    # combine duplicates by smiles and concatenate names
    size_before = df.shape[0]
    df = df.groupby("smiles").agg({"name": lambda x: " ".join(x)}).reset_index()
    print(f" Removed {size_before - df.shape[0]} duplicates")
    process_log.append(["removed_duplicates", size_before - df.shape[0]])

    print(f"Writing {df.shape[0]} molecules to {output_name}{'_fragments' if args.divide_disconnected_smiles else ''}_standardized.smi")
    process_log.append(["final_molecules", df.shape[0]])

    df[['smiles', 'name']].to_csv(f"{output_name}{'_fragments' if args.divide_disconnected_smiles else ''}_standardized.smi", header=None, index=False, sep=" ")

    with open(f"{output_name}{'_fragments' if args.divide_disconnected_smiles else ''}_standardized_process_log.txt", "w") as f:
        for item in process_log:
            f.write(f"{item[0]}: {item[1]}\n")