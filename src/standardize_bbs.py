# This material is part of the SAVI-Space created by Malte Korn and Matthias Rarey
# at the Center for Bioinformatics, University of Hamburg, with support from
# Marc Nicklaus (NIH,NCI), Phil Judson, Raphael Klein (BioSolveIT GmbH) and
# Christian Lemmen (BioSolveIT GmbH).
# Philip Judson and Marc Nicklaus provided the LHASA transform rules and assisted with their knowledge about SAVI-2020.
# SAVI-Space and all its components including this file/directory is licensed under CC-BY 4.0.

import argparse
import pandas as pd
from pandarallel import pandarallel
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

def naomi_check(file, naomi_check_binary):
    file_path = Path(file)
    naomi_check_binary
    p = run([naomi_check_binary,
        "--input", file,
        "--smarts", "[*]",
        "--out", file_path.parent / "temp.csv"], stdout=PIPE, stderr=PIPE)
    df = pd.read_csv(file_path.parent / "temp.csv")
    (file_path.parent / "temp.csv").unlink()
    # check if dataframe have two or one columns and name them name and ID
    if df.shape[1] == 3:
        df[["molecule","name"]].to_csv(file, header=None, index=False, sep=" ")
    elif df.shape[1] == 2:
        df["molecule"].to_csv(file, header=None, index=False, sep=" ")
    else:
        raise ValueError("Input file must have one or two columns")

def divide_disconnected_smiles(file):
    with open(file, "r") as f:
        smiles = f.readlines()
    results = dict()
    with open(f"{file[:-4]}_fragments.smi", "w") as f:
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
    parser.add_argument("--n_jobs", type=int, default=1, help="Number of jobs to run")
    args = parser.parse_args()


    pandarallel.initialize(progress_bar=True, nb_workers=args.n_jobs)

    process_log = []

    if not args.verbose:
        RDLogger.DisableLog('rdApp.*')

    name = args.input[:-4]
    file = args.input
    file_type = args.input[-3:]

    if file_type == "sdf":
        sdf_entrys = 0
        sdf_valid_entrys = 0
        if not Path(f"{name}.smi").exists():
            mols = Chem.SDMolSupplier(file)
            file = f"{name}.smi"
            with open(file, "w") as f:
                for mol in tqdm(mols, desc="Converting sdf to smi"):
                    sdf_entrys += 1
                    if mol is not None:
                        sdf_valid_entrys += 1
                        f.write(f"{Chem.MolToSmiles(mol)} {mol.GetProp('ID')}\n")
        print(f"Converted {sdf_valid_entrys}/{sdf_entrys} valid sdf entrys to smi")
        process_log.append(["sdf_entrys", sdf_entrys])
        process_log.append(["sdf_valid_entrys", sdf_valid_entrys])
    elif file_type == "smi":
        smiles_entrys = 0
        smiles_valid_entrys = 0
        with open(file, "r") as f:
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
        path = Path(f"{name}_fragments.smi")
        if not path.exists():
            divide_disconnected_smiles(file)
        smiles_entrys = 0
        with open(path, "r") as f:
            smiles_entrys = len(f.readlines())
        print(f"Diveded {smiles_entrys} smiles")
        process_log.append(["divided_smiles", smiles_entrys])

    if args.naomi_check_binary is not None:
        path = Path(f"{name}{'_fragments' if args.divide_disconnected_smiles else ''}.smi")
        naomi_check(path, args.naomi_check_binary)
        smiles_entrys = 0
        with open(path, "r") as f:
            smiles_entrys = len(f.readlines())
        print(f"{smiles_entrys} valid smiles with naomi_check")
        process_log.append(["naomi_valid_smiles", smiles_entrys])


    df_name = f"{name}{'_fragments' if args.divide_disconnected_smiles else ''}.csv"
    if not Path(df_name).exists():
        df = pd.read_csv(f"{name}{'_fragments' if args.divide_disconnected_smiles else ''}.smi", header=None, sep=" ", names=["smiles", "name"])
    else:
        df = pd.read_csv(df_name)

    # load mols
    df["molecule"] = df["smiles"].parallel_apply(lambda x: Chem.MolFromSmiles(x) if x is not None else np.nan)

    if (args.min_heavy_atoms is not None or args.max_heavy_atoms is not None) and "heavy_atoms" not in df.columns:
        df["heavy_atoms"] = df["molecule"].parallel_apply(lambda x: x.GetNumHeavyAtoms() if x is not None else np.nan)

    if (args.min_molecular_weight is not None or args.max_molecular_weight is not None) and "MW" not in df.columns:
        df["MW"] = df["molecule"].parallel_apply(lambda x: ExactMolWt(x) if x is not None else np.nan)
    #save
    df.drop(columns=["molecule"]).to_csv(df_name, index=False)

    # remove molecules with more than max_heavy_atoms or less than min_heavy_atoms
    if args.min_heavy_atoms is not None:
        name = f"{name}_min_heavy_atoms_{args.min_heavy_atoms}"
        df = df[df["MW"] >= args.min_heavy_atoms]
        print(f"{df.shape[0]} molecules with more than {args.min_heavy_atoms} heavy atoms")
        process_log.append(["min_heavy_atoms", df.shape[0]])

    if args.max_heavy_atoms is not None:
        name = f"{name}_max_heavy_atoms_{args.max_heavy_atoms}"
        df = df[df["MW"] <= args.max_heavy_atoms]
        print(f"{df.shape[0]} molecules with less than {args.max_heavy_atoms} heavy atoms")
        process_log.append(["max_heavy_atoms", df.shape[0]])

    if args.max_molecular_weight is not None:
        name = f"{name}_max_molecular_weight_{args.max_molecular_weight}"
        df = df[df["MW"] <= args.max_molecular_weight]
        print(f"{df.shape[0]} molecules with less than {args.max_molecular_weight} molecular weight")
        process_log.append(["max_molecular_weight", df.shape[0]])

    if args.min_molecular_weight is not None:
        name = f"{name}_min_molecular_weight_{args.min_molecular_weight}"
        df = df[df["MW"] >= args.min_molecular_weight]
        print(f"{df.shape[0]} molecules with more than {args.min_molecular_weight} molecular weight")
        process_log.append(["min_molecular_weight", df.shape[0]])

    if args.max_bertz_ct is not None:
        name = f"{name}_max_bertz_ct_{args.max_bertz_ct}"
        df = df[df["molecule"].parallel_apply(lambda x: BertzCT(x) if x is not None else np.nan) <= args.max_bertz_ct]
        print(f"{df.shape[0]} molecules with less than {args.max_bertz_ct} bertz ct")
        process_log.append(["max_bertz_ct", df.shape[0]])

    df["molecule"] = df["molecule"].parallel_apply(standardize)
    df = df.dropna(subset=["molecule"])
    print(f"{df.shape[0]} valid molecules after standardization")
    process_log.append(["standardization", df.shape[0]])
    df["smiles"] = df["molecule"].parallel_apply(lambda x: Chem.MolToSmiles(x) if x is not None else np.nan)

    # combine duplicates by smiles and concatenate names
    size_before = df.shape[0]
    df = df.groupby("smiles").agg({"name": lambda x: " ".join(x)}).reset_index()
    print(f" Removed {size_before - df.shape[0]} duplicates")
    process_log.append(["removed_duplicates", size_before - df.shape[0]])

    print(f"Writing {df.shape[0]} molecules to {name}{'_fragments' if args.divide_disconnected_smiles else ''}_standardized.smi")
    process_log.append(["final_molecules", df.shape[0]])

    df[['smiles', 'name']].to_csv(f"{name}{'_fragments' if args.divide_disconnected_smiles else ''}_standardized.smi", header=None, index=False, sep=" ")

    with open(f"{name}{'_fragments' if args.divide_disconnected_smiles else ''}_standardized_process_log.txt", "w") as f:
        for item in process_log:
            f.write(f"{item[0]}: {item[1]}\n")