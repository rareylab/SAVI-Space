# This material is part of the SAVI-Space created by Malte Korn and Matthias Rarey
# at the Center for Bioinformatics, University of Hamburg, with support from
# Marc Nicklaus (NIH,NCI), Phil Judson, Raphael Klein (BioSolveIT GmbH) and
# Christian Lemmen (BioSolveIT GmbH).
# Philip Judson and Marc Nicklaus provided the LHASA transform rules and assisted with their knowledge about SAVI-2020.
# SAVI-Space and all its components including this file/directory is licensed under CC-BY 4.0.

import re

from rdkit import Chem
from rdkit.Chem import AllChem
import logging

from utils import mask_exocyclic_doublebonds, smarts_exocyclic_doublebonds, get_rxn_smarts


def kill(smiles, smarts, kill_statements, disallow_exocyclic_doublebonds=False):
    kill_statements = kill_statements["kill"]
    mol = Chem.MolFromSmiles(smiles)
    if disallow_exocyclic_doublebonds:
        mol = mask_exocyclic_doublebonds.normalize(mol)
        smarts = smarts_exocyclic_doublebonds(smarts)
    partition_rxn_smarts, marked_rxn_smarts, atom_map = get_rxn_smarts(smarts)
    # create partitions map based on atom_map
    partition_rxn = AllChem.ReactionFromSmarts(partition_rxn_smarts)
    partitions = {}
    for atom_map_idx in atom_map:
        partitions[str(atom_map_idx)] = []

    try:
        partition_outcome = partition_rxn.RunReactants((mol,))
    except:
        print(f"Error while partitioning {Chem.MolToSmiles(mol)} with {partition_rxn_smarts}")
        logging.error(f"Error while partitioning {Chem.MolToSmiles(mol)} with {partition_rxn_smarts}")
        return {"partition": False}
    if not partition_outcome:
        print(f"Error no partitioning output for {Chem.MolToSmiles(mol)} with {partition_rxn_smarts}")
        logging.error(f"Error no partitioning output for {Chem.MolToSmiles(mol)} with {partition_rxn_smarts}")
        return {"partition": False}
    for part in partition_outcome:
        for atom_map_idx, p in zip(atom_map, part):
            # TODO: workaround because of fragments with aromatic bonds and atoms not in ring
            Chem.SanitizeMol(p,sanitizeOps=Chem.SanitizeFlags.SANITIZE_NONE) # &Chem.SanitizeFlags.SANITIZE_KEKULIZE)
            Chem.GetSSSR(p)
            partitions[str(atom_map_idx)].append(p)

    # mark atoms based on SMARTS atom_map
    marked_rxn = AllChem.ReactionFromSmarts(marked_rxn_smarts)
    try:
        marked_reactant = [m[0] for m in marked_rxn.RunReactants((mol,))]
    except:
        print(f"Error in marking {Chem.MolToSmiles(mol)}")
        logging.error(f"Error in marking {Chem.MolToSmiles(mol)}")
        return {"marking": False}
    for m in marked_reactant:
        Chem.SanitizeMol(m, sanitizeOps=Chem.SanitizeFlags.SANITIZE_NONE) # &Chem.SanitizeFlags.SANITIZE_KEKULIZE)
        Chem.GetSSSR(m)

    result = apply_kill_statements(kill_statements, atom_map, 0, partitions, marked_reactant, disallow_exocyclic_doublebonds=disallow_exocyclic_doublebonds)
    return result


def apply_kill(kill_pattern, atom_map, idx, partitions, marked_reactant, disallow_exocyclic_doublebonds=False):
    success = False
    atom_idx = kill_pattern["atom_idx"]


    if kill_pattern["OFFPATH"] == True:
        if kill_pattern["atom_idx"] == "anywhere":
            atom_idx = atom_map
        for pattern_idx, atom_idx in enumerate(atom_idx):
            if kill_pattern["atom_idx"] == "anywhere":
                pattern_idx = 0
            if not atom_idx in atom_map:
                continue
            for part in partitions[str(atom_idx)]:
                patt_smarts = kill_pattern["smarts"][pattern_idx]
                if disallow_exocyclic_doublebonds:
                    patt_smarts = smarts_exocyclic_doublebonds(patt_smarts)
                patt = Chem.MolFromSmarts(patt_smarts)
                if not patt:
                    raise Exception(f"Invalid SMARTS pattern: {kill_pattern['smarts'][pattern_idx]}; kill_idx: {idx}; pattern_idx: {pattern_idx}")
                if not kill_pattern["atom_idx"] == "anywhere" and not part.HasSubstructMatch(patt):
                    success = False
                    break
                elif part.HasSubstructMatch(patt):
                    success = True
                    if kill_pattern["atom_idx"] == "anywhere":
                        break
            if not kill_pattern["atom_idx"] == "anywhere" and not success:
                break
            elif kill_pattern["atom_idx"] == "anywhere" and success:
                break
    else:
        for reactant in marked_reactant:
            patt_smarts = kill_pattern["smarts"]
            if disallow_exocyclic_doublebonds:
                patt_smarts = smarts_exocyclic_doublebonds(patt_smarts)
            patt = Chem.MolFromSmarts(patt_smarts)
            if not patt:
                raise Exception(f"Invalid SMARTS pattern: {kill_pattern['smarts']}; kill_idx: {idx}")
            try:
                if reactant.HasSubstructMatch(patt):
                    success = True
                    break
            except:
                print(f"Error in matching")
                logging.error(f"Error in matching")
    return success

def apply_kill_statements(kill_statements, atom_map, idx, partitions, marked_reactant, disallow_exocyclic_doublebonds=False):
    if idx >= len(kill_statements):
        return {f"kill_end": False}
    # apply CHMTRN rules on partitions and return True if kill SMARTS pattern matches
    kill_pattern = kill_statements[str(idx)]

    if kill_pattern["atom_idx"] == "anywhere":
        kill_atom_idx = atom_map

    else:
        kill_atom_idx = kill_pattern["atom_idx"]

    if "goto" in kill_pattern:
        if set(kill_atom_idx).issubset(atom_map):
            r = apply_kill(kill_pattern, atom_map, idx, partitions, marked_reactant, disallow_exocyclic_doublebonds=disallow_exocyclic_doublebonds)
            if r == True:
                next_pos = {f"kill_{idx}": None}
                goto_pos = apply_kill_statements(kill_statements, atom_map, kill_pattern["goto"], partitions, marked_reactant, disallow_exocyclic_doublebonds=disallow_exocyclic_doublebonds)
            else:
                next_pos = apply_kill_statements(kill_statements, atom_map, idx+1, partitions, marked_reactant, disallow_exocyclic_doublebonds=disallow_exocyclic_doublebonds)
                goto_pos = apply_kill_statements(kill_statements, atom_map, kill_pattern["goto"], partitions, marked_reactant, disallow_exocyclic_doublebonds=disallow_exocyclic_doublebonds)

        else:
            if set(kill_atom_idx).intersection(atom_map):
                r = apply_kill(kill_pattern, atom_map, idx, partitions, marked_reactant, disallow_exocyclic_doublebonds=disallow_exocyclic_doublebonds)
            else:
                r = None
            next_pos = apply_kill_statements(kill_statements, atom_map, idx+1, partitions, marked_reactant, disallow_exocyclic_doublebonds=disallow_exocyclic_doublebonds)
            goto_pos = apply_kill_statements(kill_statements, atom_map, kill_pattern["goto"], partitions, marked_reactant, disallow_exocyclic_doublebonds=disallow_exocyclic_doublebonds)
        return {f"kill_{idx}":[
                {f"kill_{idx}": r},
                goto_pos,
                next_pos
            ]}


    if not set(kill_atom_idx).issubset(atom_map) and set(kill_atom_idx).intersection(atom_map):
        temp1 = apply_kill(kill_pattern, atom_map, idx, partitions, marked_reactant, disallow_exocyclic_doublebonds=disallow_exocyclic_doublebonds)
        temp2 = apply_kill_statements(kill_statements, atom_map, idx+1, partitions, marked_reactant, disallow_exocyclic_doublebonds=disallow_exocyclic_doublebonds)
        if temp2 == True:
            return {f"kill_{idx}": temp2}
        if temp1 == False and temp2 == False:
            return {f"kill_{idx}": False}
        else:
            return {f"kill_{idx}": [
                {f"kill_{idx}": temp1},
                temp2
            ]}

    success = apply_kill(kill_pattern, atom_map, idx, partitions, marked_reactant, disallow_exocyclic_doublebonds=disallow_exocyclic_doublebonds)
    if success:
        return {f"kill_{idx}": True}
    else:
        return apply_kill_statements(kill_statements, atom_map, idx+1, partitions, marked_reactant, disallow_exocyclic_doublebonds=disallow_exocyclic_doublebonds)


def create_pairs(adata, bdata):
    adata_values = adata[list(adata.keys())[0]]
    bdata_values = bdata[list(bdata.keys())[0]]
    if not list(adata.keys())[0] == list(bdata.keys())[0]:
        # if adata_values == True or bdata_values == True:
        #     return False
        return False
    if adata_values == False or adata_values == None:
        if bdata_values == False or bdata_values == None:
            return list(adata.keys())[0]
        elif isinstance(bdata_values, dict):
            if len(bdata_values) == 2:
                if create_pairs(adata, bdata_values[0]):
                    return list(adata.keys())[0]
                return create_pairs(adata, bdata_values[1])
            elif len(bdata_values) == 3:
                if list(bdata_values[0].values())[0] == True:
                    return create_pairs(adata, bdata_values[1])
                return create_pairs(adata, bdata_values[2])

    elif not isinstance(adata_values, bool):
        if len(adata_values) == 2:
            if bdata_values == False:
                return create_pairs(adata_values[0], bdata)
            elif bdata_values == True:
                return create_pairs(adata_values[1], bdata)
            elif not isinstance(bdata_values, bool):
                if len(bdata_values) == 2:
                    if list(adata_values[0].values())[0] == True and list(bdata_values[0].values())[0] == True:
                        return False
                    return create_pairs(adata_values[1], bdata_values[1])
                if len(bdata_values) == 3:
                    return False
        elif len(adata_values) == 3:
            if bdata_values == False:
                if list(adata_values[0].values())[0] == False:
                    return create_pairs(adata_values[2], bdata)
                return create_pairs(adata_values[1], bdata)
            elif not isinstance(bdata_values, bool):
                if len(bdata_values) == 3:
                    if list(bdata_values[0].values())[0] == True or list(adata_values[0].values())[0] == True:
                        return create_pairs(adata_values[1], bdata_values[1])
                    else:
                        return create_pairs(adata_values[2], bdata_values[2])
                if len(bdata_values) == 2:
                    return False
    return False