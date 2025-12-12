# This material is part of the SAVI-Space created by Malte Korn and Matthias Rarey
# at the Center for Bioinformatics, University of Hamburg, with support from
# Marc Nicklaus (NIH,NCI), Phil Judson, Raphael Klein (BioSolveIT GmbH) and
# Christian Lemmen (BioSolveIT GmbH).
# Philip Judson and Marc Nicklaus provided the LHASA transform rules and assisted with their knowledge about SAVI-2020.
# SAVI-Space and all its components including this file/directory is licensed under CC-BY 4.0.

import json
from pathlib import Path
import re
import logging
import sys

from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

def initialize_logger(log_file, console_level=logging.INFO):
    """
    Initializes and configures the **Root Logger** of the application.

    Args:
        console_level (int): The logging level for the console output 
                             (e.g., logging.INFO, logging.WARNING).
        log_file (str): The name of the file to write logs to.
    """
    # *** WICHTIGE ÄNDERUNG: Aufruf ohne Argumente gibt den Root Logger zurück ***
    logger = logging.getLogger() 
    
    # Der Root Logger wird standardmäßig auf WARNING gesetzt. 
    # Um DEBUG-Meldungen zu verarbeiten, muss das Level hochgesetzt werden.
    logger.setLevel(logging.DEBUG) 

    # Prevent handlers from being added multiple times if called more than once
    if logger.hasHandlers():
        return logger

    formatter = logging.Formatter(
        '%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

    # File Handler
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    # Console Handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(console_level)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)
    
    return logger


# remove tree with option to keep files, when delete_files=False only empty folders are removed
def remove_tree(path, delete_files=True):
    path = Path(path)
    remove_branche = True
    if path.is_dir():
        for child in path.iterdir():
            if child.is_dir():
                remove_branche = remove_tree(child)
            else:
                if delete_files:
                    child.unlink()
                else:
                    remove_branche = False
                    print("File not deleted: ", child)
        if remove_branche:
            path.rmdir()
        else:
            print("Folder not deleted: ", path, " (not empty)")
    else:
        if delete_files:
            path.unlink()
        else:
            print("File not deleted: ", path)
            remove_branche = False

    return remove_branche


# disabel tqdm for submissions to the cluster
class tqdm:
    def __init__(self, iterable, leave=True, desc=None, total=None, file=None, mininterval=0.1, miniters=1, ascii=None, disable=False, unit='it', unit_scale=False, dynamic_ncols=False, smoothing=0.3, bar_format=None, initial=0, position=None, postfix=None, unit_divisor=1000):
        self.iterable = iterable
    def __iter__(self):
        return iter(self.iterable)
    def __len__(self):
        return len(self.iterable)
    def set_description(self, string):
        print(string)
    def write(self, string):
        print(string)


def get_filtered_bbs_path(outpath, reaction_index, reaction_pattern):
    return (outpath / f"processed_bbs/{reaction_index}_{reaction_pattern[reaction_index]['name']}")

def get_pattern_smarts(reactions, reaction_index, unique_matches=None):
    pattern_smarts = dict()
    reactant_patterns = [[r.split(">>")[0].split(".") for r in reaction] for reaction in reactions[reaction_index]]
    for idx1, reactant_pattern in enumerate(reactant_patterns):
        for idx2, r in enumerate(reactant_pattern):
            for reactant_index in range(2):
                pattern_name = f"{reaction_index}_{idx1}_{int_to_char(idx2)}_r{reactant_index+1}"
                if unique_matches:
                    pattern_smarts[pattern_name] = [r[reactant_index], unique_matches[reaction_index][idx1][reactant_index]]
                else:
                    pattern_smarts[pattern_name] = r[reactant_index]
    return pattern_smarts

def load_reactions(filename='smirks.json'):
    script_dir = Path(__file__).resolve().parent
    project_root = script_dir.parent.parent
    filepath_project = project_root / 'data' / filename
    filepath = None
    if filepath_project.is_file():
        filepath = filepath_project
    elif Path(filename).is_file():
        filepath = Path(filename)
    if filepath is None:
        raise FileNotFoundError(f"{filename} could not be found.")
    with open(filepath) as json_file:
        reaction_pattern = json.load(json_file)
    reactions = {}
    unique_matches = {}

    for reaction_index in reaction_pattern:
        if reaction_index.startswith("_"):
                continue
        reactions[reaction_index] = []
        unique_matches[reaction_index] = []
        for idx, patterns in enumerate(reaction_pattern[reaction_index]["smirks"]):
            reactions[reaction_index].append([])
            for pattern in patterns:
                reactions[reaction_index][-1].append(pattern)
            unique_matches[reaction_index].append(reaction_pattern[reaction_index]["unique_matching"][idx])
    return reactions, unique_matches

def get_reaction_name(reaction_idx, filename='smirks.json'):
    with open(f"data/{filename}") as json_file:
        patterns_json = json.load(json_file)
    return patterns_json[str(reaction_idx)]["name"]

def load_pattern_json(filename='smirks.json'):
    script_dir = Path(__file__).resolve().parent
    project_root = script_dir.parent.parent
    filepath_project = project_root / 'data' / filename
    filepath = None
    if filepath_project.is_file():
        filepath = filepath_project
    elif Path(filename).is_file():
        filepath = Path(filename)
    if filepath is None:
        raise FileNotFoundError(f"{filename} could not be found.")
    with open(filepath) as json_file:
        patterns_json = json.load(json_file)
    # delete comment from json
    patterns_json = {key: patterns_json[key] for key in patterns_json if not key.startswith("_")}
    return patterns_json


def load_kill_patterns(filename="KILL.json"):
    script_dir = Path(__file__).resolve().parent
    project_root = script_dir.parent.parent
    filepath_project = project_root / 'data' / filename
    filepath = None
    if filepath_project.is_file():
        filepath = filepath_project
    elif Path(filename).is_file():
        filepath = Path(filename)
    if filepath is None:
        raise FileNotFoundError(f"{filename} could not be found.")
    with open(filepath) as json_file:
        kill_pattern = json.load(json_file)
    # delete comment from json
    kill_pattern = {key: kill_pattern[key] for key in kill_pattern if not key.startswith("_")}
    return kill_pattern

def int_to_char(integer):
    if integer < 0 or integer > 25:
        raise Exception(f"Integer {integer} out of range (0-25)")
    return chr(integer+97)

def indent_multiline_string(string, indent_length=26):
    return string.replace("\n", "\n" + " "*indent_length)


_normalizer_params = rdMolStandardize.CleanupParameters()
mask_exocyclic_doublebonds_pattern = """
//	Name	SMIRKS
Exocyclic doublebonds\t[c:1](=[O,S,N:2])>>[Ge:1](=[*:2])
"""
mask_exocyclic_doublebonds = rdMolStandardize.NormalizerFromData(
    mask_exocyclic_doublebonds_pattern, _normalizer_params
)
unmask_exocyclic_doublebonds_pattern = """
//	Name	SMIRKS
Unmask exocyclic doublebonds\t[Ge:1](=[O,S,N:2])>>[#6:1](=[*:2])
"""
unmask_exocyclic_doublebonds = rdMolStandardize.NormalizerFromData(
    unmask_exocyclic_doublebonds_pattern, _normalizer_params
)

def smarts_exocyclic_doublebonds(smarts):
    return smarts.replace("#6", "#6,Ge").replace("!#6,Ge", "!#6;!Ge").replace("#6,Ge;a", "#6;a")


def re_reaction_replacement(matchobj):
    number = re.findall(r"(?<=:)[0-9]+", matchobj.group(0))[0]
    return f"*;{number}:{number}"


def get_rxn_smarts(reactant_smarts):
    pattern = Chem.MolFromSmarts(reactant_smarts)
    atom_map = [atom.GetAtomMapNum() for atom in pattern.GetAtoms() if atom.GetAtomMapNum() > 0]
    partition_rxn_smarts = reactant_smarts + ">>" + ".".join([f"[*{i}:{i}]" for i in atom_map])

    # remove smarts stuff from pattern
    smiles = Chem.MolToSmiles(pattern)
    smiles = re.sub(r"(?<=\[)[A-Z\*]+[0-9]*:[0-9]+", re_reaction_replacement, smiles)
    smiles = re.sub(r"\[?[A-GI-Z][a-z]?[H]?[0-9]*[a-z]?(?!:)\]?", "[#0]", smiles) # do not match H
    marked_rxn_smarts = reactant_smarts + ">>" + smiles
    return partition_rxn_smarts, marked_rxn_smarts, atom_map


def killing(educts, reactant_idx, sub_reaction_name, kill_name):
    return ~educts[f"{sub_reaction_name}_r{reactant_idx}_{kill_name}"].isna()

def single_matching(educts, reactant_idx, sub_reaction_name, kill_name="no_kill"):
    if kill_name != "no_kill":
        return (educts[f"{sub_reaction_name}_r{reactant_idx}"] == 1) &\
                ((educts[f"{sub_reaction_name}_r{(reactant_idx%2)+1}"] == 0) |\
                ~(educts[f"{sub_reaction_name}_r{(reactant_idx%2)+1}_{kill_name}"].isna()))
    else:
        return (educts[f"{sub_reaction_name}_r{reactant_idx}"] == 1) & (educts[f"{sub_reaction_name}_r{(reactant_idx%2)+1}"] == 0)

def get_valid(enamine_educts, reactant_idx, sub_reaction_name, kill_name="no_kill"):
    single_match = single_matching(enamine_educts, reactant_idx, sub_reaction_name, kill_name)
    if kill_name != "no_kill":
        killed = killing(enamine_educts, reactant_idx, sub_reaction_name, kill_name)
        return(enamine_educts[~killed & single_match])
    return enamine_educts[single_match]