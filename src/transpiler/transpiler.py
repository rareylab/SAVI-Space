# This material is part of the SAVI-Space created by Malte Korn and Matthias Rarey
# at the Center for Bioinformatics, University of Hamburg, with support from
# Marc Nicklaus (NIH,NCI), Phil Judson, Raphael Klein (BioSolveIT GmbH) and
# Christian Lemmen (BioSolveIT GmbH).
# Philip Judson and Marc Nicklaus provided the LHASA transform rules and assisted with their knowledge about SAVI-2020.
# SAVI-Space and all its components including this file/directory is licensed under CC-BY-NC 4.0.

import re, json, copy, argparse, itertools

from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import AllChem

# set rdkit log level to error
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

# load keywords global
try:
    with open('src/transpiler/chmtrn_keywords_to_smarts.json') as f:
        keywords = json.load(f)
except:
    try:
        with open('transpiler/chmtrn_keywords_to_smarts.json') as f:
            keywords = json.load(f)
    except:
        with open('chmtrn_keywords_to_smarts.json') as f:
            keywords = json.load(f)


def check_patran_syntax(patran_string):
    def patran_syntax_error(text, patran_string, pos):
        raise SyntaxError(f"{text}\n{patran_string}\n" + "-" * pos + "^")

    def check_for_unknown_value(value, key, patran_string, error_pos):
        if value not in keywords[key]:
            patran_syntax_error(f"Unknown {key}: {value}",
                                            patran_string, error_pos)

    def atom(parser, token): return "ATOM", token.split(","), len(token)
    def atomindex(parser, token): return "INDEX", token[1:], len(token)
    def atomproperty(paser, token):
        output = {}
        for t in token[1:-1].split(";"):
            key, sign, values = re.split(r"(<|>|=|#)",t)
            output[key] = {"sign": sign, "values": values.split(",")}
        return "PROPERTY", output, len(token)
    def fusion(parser, token):
        split = re.split(r"(=|#)", token[1:-1])
        return "FUSION", split[1:3], len(token)
    def end_sidechain(parser, token): return "END_SIDECHAIN", True, len(token)
    def start_sidechain(parser, token): return "START_SIDECHAIN", True, len(token)
    def bond(parser, token): return "BOND", token.split(","), len(token)
    def ring_closure(parser, token): return "RING", int(token[1:]), len(token)

    parser = re.Scanner([
        (r"[A-Za-z,]+", atom),
        (r"\^[0-9]+", atomindex),
        (r"\[(?!FUSION)[A-Za-z0-9\*=,;\<\>]+\]", atomproperty),
        (r"\[FUSION[=#][A-Z]+\]", fusion),
        (r"\(", start_sidechain),
        (r"\)", end_sidechain),
        (r"[-=#%&\+](,[-=#%&\+])*", bond),
        (r"@[0-9]+", ring_closure)
    ])

    syntax_rules = {
        "after":{
            "ATOM": ["INDEX", "PROPERTY", "BOND", "END_SIDECHAIN", "START_SIDECHAIN"],
            "INDEX": ["BOND", "PROPERTY", "START_SIDECHAIN", "END_SIDECHAIN"],
            "PROPERTY": ["BOND", "START_SIDECHAIN", "END_SIDECHAIN"],
            "FUSION": ["ATOM", "RING"],
            "BOND": ["ATOM", "RING", "FUSION"],
            "START_SIDECHAIN": ["BOND"],
            "END_SIDECHAIN": ["START_SIDECHAIN", "BOND"],
            "RING": ["END_SIDECHAIN", "BOND"]
        },
        "before": {
            "ATOM": ["BOND", "FUSION"],
            "INDEX": ["ATOM"],
            "PROPERTY": ["ATOM", "INDEX"],
            "FUSION": ["BOND"],
            "BOND": ["ATOM", "INDEX", "PROPERTY", "START_SIDECHAIN", "END_SIDECHAIN"],
            "START_SIDECHAIN": ["PROPERTY", "ATOM", "END_SIDECHAIN", "INDEX"],
            "END_SIDECHAIN": ["ATOM", "PROPERTY", "RING", "INDEX"],
            "RING": ["BOND","FUSION"]
        },
        "not_end_with": ["BOND", "FUSION"],
        "not_begin_with": ["BOND", "FUSION", "INDEX", "RING",
                        "PROPERTY", "START_SIDECHAIN", "END_SIDECHAIN"]
    }

    tokens, remainder = parser.scan(patran_string.strip())

    # check for unparses
    if remainder:
        expexted_syntax = " or ".join(["<"+ elem + ">"\
                         for elem in syntax_rules["after"][tokens[-1][0]]])
        error_pos = len(patran_string)-len(remainder)
        patran_syntax_error(f"Unknown pattern; expect {expexted_syntax}",
                                         patran_string, error_pos)

    # check start and end of string
    if tokens[0][0] in syntax_rules["not_begin_with"]:
        patran_syntax_error(f"Not start with {tokens[0][0]}", patran_string, 0)
    elif tokens[-1][0] in syntax_rules["not_end_with"]:
        patran_syntax_error(f"Not end with {tokens[-1][0]}",
                        patran_string, len(patran_string)- tokens[-1][2])

    open_sidechains = 0
    error_pos = 0
    for idx, token in enumerate(tokens):
        # check START_SIDECHAIN and END_SIDECHAIN
        if token[0] == "START_SIDECHAIN":
            open_sidechains += 1
        elif token[0] == "END_SIDECHAIN":
            if open_sidechains < 1:
                patran_syntax_error(
                    "there is no sidechain open to close",
                    patran_string, error_pos)
            open_sidechains -= 1

        #check allowed keywords
        elif token[0] == "ATOM":
            for atom in token[1]:
                check_for_unknown_value(atom.upper(), "atom",
                                         patran_string, error_pos)
        elif token[0] == "FUSION":
            check_for_unknown_value(token[1][1], "fusion",
                                        patran_string, error_pos)
            if token[1][0] not in ["=", "#"]:
                patran_syntax_error(f"Unknown sign:\
                                 {token[1][0]}", patran_string, error_pos)
        elif token[0] == "PROPERTY":
            for prop in token[1]:
                check_for_unknown_value(prop, "properties",
                                            patran_string, error_pos)
                #check different props
                if prop == "FGS" or prop == "FGNOT":
                    for groupname in token[1][prop]["values"]:
                        if groupname not in keywords["FG"]\
                            and groupname not in keywords["synonym"]:
                            patran_syntax_error(
                                f"unknown functional group: {groupname}",
                                patran_string, error_pos)
                if prop == "HS" or prop == "HETS":
                    for propvalue in token[1][prop]["values"]:
                        if not propvalue.isdigit() or\
                             int(propvalue) < 0 or int(propvalue) > 4:
                            patran_syntax_error(
                                "value must be integer and between 0 and 4;"\
                                f"current value: {propvalue}",
                                patran_string, error_pos)
                if prop == "ARYL" or prop == "PHENYL" or prop == "RINGS":
                    for typo in token[1][prop]["values"]:
                        if typo not in ["NO", "YES", "EITHER"]:
                            patran_syntax_error(
                                f"unknown value {typo}",
                                patran_string, error_pos)
        elif token[0] == "BOND":
            for bond in token[1]:
                check_for_unknown_value(bond, "bond", patran_string, error_pos)

        #check syntax order
        error_pos += token[2]
        if idx < (len(tokens)-1)\
                and tokens[idx+1][0] not in syntax_rules["after"][token[0]]:
            expexted_syntax = " or ".join(["<"+ elem + ">"\
                             for elem in syntax_rules["after"][token[0]]])
            patran_syntax_error(f"Get <{tokens[idx+1][0]}>,"\
                         f"expect {expexted_syntax}", patran_string, error_pos)

    if open_sidechains > 0:
        raise SyntaxError("There are open sidechains")

    return tokens


def patran_to_graph(parsed_patran):

    def new_node(atom, atom_map, prop, sidechain):
        num_neighbours =  0  # placeholder for number of adj. nodes
        return {"atom": atom,
                "atom_map": atom_map,
                "prop": prop,
                "sidechain": sidechain,
                "num_neighbours": num_neighbours}
    def new_edge(s_node, e_node, bondtype, fusion):
        return{"start": s_node,
               "end": e_node,
               "type": bondtype,
               "fusion": fusion}

    sidechain_indices = [0] # stack of indeces to enumerate the sidechains; 0 -> mainchain
    current_sidechain_index = 0 # number to specify index of the last new sidechain
    node_indices = [0] # stack of node indices to handle start and end node for edges
    current_node_index = 0
    nodes = []
    edges = []
    idx = 0
    while idx < len(parsed_patran):
        # transform node
        if parsed_patran[idx][0] == "ATOM":
            atom = parsed_patran[idx][1]
            idx += 1
            atom_map = None
            if idx < len(parsed_patran) and parsed_patran[idx][0] == "INDEX":
                atom_map = int(parsed_patran[idx][1])
                idx += 1
            prop = []
            if idx < len(parsed_patran) and parsed_patran[idx][0] == "PROPERTY":
                prop = parsed_patran[idx][1]
                idx += 1
            sidechain = sidechain_indices[-1]
            nodes.append(new_node(atom, atom_map, prop, sidechain))
            current_node_index += 1

        # transform edge
        elif parsed_patran[idx][0] == "BOND":
            bondtype = parsed_patran[idx][1]
            s_node = node_indices.pop()
            if bondtype == ["+"]:
                # disconnected bonds are always start at the last node
                s_node = current_node_index-1
            node_indices.append(current_node_index)
            e_node = current_node_index
            idx += 1
            fusion = None
            if idx < len(parsed_patran) and parsed_patran[idx][0] == "FUSION":
                fusion = parsed_patran[idx][1]
                idx += 1
            if idx < len(parsed_patran) and parsed_patran[idx][0] == "RING":
                e_node = parsed_patran[idx][1] - 1
                idx += 1
            edges.append(new_edge(s_node, e_node, bondtype, fusion))

        # handle start of sidechain
        elif parsed_patran[idx][0] == "START_SIDECHAIN":
            #consider whether sidechain only contains ringclosing without additional atom(s)
            # idx+1 is "BOND" and idx+2 or idx+3 might be "RING" if present
            if parsed_patran[idx+2][0] != "RING" and  parsed_patran[idx+3][0] != "RING":
                current_sidechain_index += 1
            sidechain_indices.append(current_sidechain_index)
            node_indices.append(current_node_index-1)
            idx += 1

        elif parsed_patran[idx][0] == "END_SIDECHAIN":
            sidechain_indices.pop()
            node_indices.pop()
            idx += 1

        else:
            raise SyntaxError("Something went wrong")

    # add number of adj. nodes to each node
    for node in list(itertools.chain.from_iterable([[edge["start"],edge["end"]]\
                                                        for edge in edges])):
        nodes[node]["num_neighbours"] += 1

    return {"nodes": nodes, "edges": edges}


def transform_properties_patran_to_smarts(props, element):
    def unknown_value_error(value, keyword):
        raise ValueError(f"'{value}' is not supported "\
                             f"with the keyword '{keyword}'.")

    def unknown_sign_error(sign):
        raise SyntaxError(f"Unknown sign: '{sign}'.")

    smarts_props = []

    for prop in props:

        if prop == "CHARGE":
            for value in props[prop]["values"]:
                if value in list(keywords["charge"]):
                    smarts_props.append(keywords["charge"][value])
                else:
                    unknown_value_error(value, prop)

        elif prop == "HETS":
            values = [int(v) for v in props[prop]["values"]]
            for value in values:
                if value < 0 or value > 4:
                    raise ValueError(f"values({value}) must be between 0 an 4 for {prop}")
            if props[prop]["sign"] == "=":
                temp_props = []
                values.sort()
                for value in values:
                    if value == 0:
                        temp_props.append("!$(*(~[!#6;!#1]))")
                    elif value == values[-1]:
                        temp_props.append("$(*" + "(~[!#6;!#1])" * value +");"\
                                    "!$(*" + "(~[!#6;!#1])" * (value +1) +")")
                    else:
                        temp_props.append("$(*" + "(~[!#6;!#1])" * value +")")
                smarts_props.append(",".join(temp_props))
            elif props[prop]["sign"] == ">":
                if len(values) != 1:
                    raise ValueError(f"Only one value allowed for '{props[prop]['sign']}'")
                smarts_props.append("$(*" + "(~[!#6;!#1])" * (values[0] + 1) + ")")
            elif props[prop]["sign"] == "<":
                if len(values) != 1:
                    raise ValueError(f"Only one value allowed for '{props[prop]['sign']}'")
                smarts_props.append("!$(*" + "(~[!#6;!#1])" * values[0] +")")
            else:
                unknown_sign_error(props[prop]["sign"])

        elif prop == "HS":
            temp_props = []
            for value in props[prop]["values"]:
                value = int(value)
                if value < 0 or value > 4:
                    raise ValueError(f"values({value})"/
                                " must be between 0 an 4 for {prop}")
                if props[prop]["sign"] == "=":
                    temp_props.append(f"H{value}")
                elif props[prop]["sign"] == ">":
                    for i in range(value+1):
                        temp_props.append(f"!H{i}")
                elif props[prop]["sign"] == "<":
                    for i in range(value):
                        temp_props.append(f"H{i}")
                else:
                    unknown_sign_error(props[prop]["sign"])
            smarts_props.append(",".join(temp_props))

        elif prop == "ARYL":
            for value in props[prop]["values"]:
                if value == "NO":
                    smarts_props.append("A")
                elif value == "YES":
                    smarts_props.append("a")
                elif value == "EITHER":
                    pass
                else:
                    unknown_value_error(value, prop)

        elif prop == "EPS":
            if len(element) == 1 and element[0] != "X":
                for value in props[prop]["values"]:
                    if keywords['valence'][element[0]] < (int(value) * 2):
                        raise ValueError(f"Number of electronpairs ({value})"\
                                        "is higher then atom valence"\
                            f"{keywords['valence'][element[0]]}"\
                        f"-> max {int(keywords['valence'][element[0]]/2)} EPS")
                    electron_pairs = keywords['valence'][element[0]] -\
                                                        int(value) * 2
                    smarts_props.append(f"v{electron_pairs}")
            else:
                raise SyntaxError(f"Electronpairs are only valid"\
                                            "for one explicit atom")

        elif prop == "RINGS":
            for value in props[prop]["values"]:
                if value == "NO" or value == "0":
                    smarts_props.append("!R")
                elif value == "YES":
                    smarts_props.append("R")
                elif int(value) > 2:
                    smarts_props.append(f"r{value}")
                elif value == "EITHER":
                    pass
                else:
                    unknown_value_error(value, prop)

        elif prop == "PHENYL":
            for value in props[prop]["values"]:
                if value == "YES":
                    smarts_props.append(f"$({keywords['FG']['PHENYL']})")
                elif value == "NO":
                    smarts_props.append(f"!$({keywords['FG']['PHENYL']})")
                else:
                    unknown_value_error(value, prop)

        elif prop == "T*BUTYL":
            for value in props[prop]["values"]:
                if value == "YES":
                    smarts_props.append(f"$({keywords['FG']['T*BUTYL']})")
                elif value == "NO":
                    smarts_props.append(f"!$({keywords['FG']['T*BUTYL']})")
                else:
                    unknown_value_error(value, prop)

        elif (prop == "FGS" or prop == "FGNOT"):
            temp_props = []
            for value in props[prop]["values"]:
                current_value = [value]
                if value in list(keywords["synonym"]):
                    current_value = keywords["synonym"][value]
                if prop == "FGS":
                    for val in current_value:
                        temp_props.append(f"$({keywords['FG'][val]})")
                else:
                    for val in current_value:
                        temp_props.append(f"!$({keywords['FG'][val]})")
            if prop == "FGS":
                # FGS mustr be connected by 'or' statement
                smarts_props.append(",".join(temp_props))
            else:
                # FGNOT must be connected by 'and' statement
                smarts_props.append(";".join(temp_props))

        else:
            raise ValueError(f"Not supported property key: '{prop}'")

        # remove empty strings from smarts_props
        smarts_props = list(filter(bool, smarts_props))

    return smarts_props


def graph_patran_to_smarts(patran_graph, product=False):
    # create list of a deepcopy of the patran graph to create multiple smarts graphs
    smarts_graph_list = [copy.deepcopy(patran_graph)]

    for idx, node in enumerate(patran_graph["nodes"]):
        temp_atom = []
        if product and len(node["atom"]) > 1:
            smarts_graph_list[0]["nodes"][idx]["atom"] = "*"
        else:
            for atom in node["atom"]:
                temp_atom.append(keywords["atom"][atom.upper()])
            smarts_graph_list[0]["nodes"][idx]["atom"] = temp_atom

        if node["prop"]:
            if product:
                smarts_graph_list[0]["nodes"][idx]["prop"] = []
            else:
                smarts_graph_list[0]["nodes"][idx]["prop"] = \
                    transform_properties_patran_to_smarts(node["prop"], node["atom"])

    for idx, edge in enumerate(patran_graph["edges"]):
        temp_bondtype = []
        if product:
            if len(edge["type"]) > 1:
                smarts_graph_list[0]["edges"][idx]['type'] = ["~"]
                continue
        for bondtype in edge['type']:
            temp_bondtype.append(keywords["bond"][bondtype])

        smarts_graph_list[0]["edges"][idx]['type'] = list(set(temp_bondtype))

        # handle ring fusion
        if edge["fusion"] and not product:
            if edge["fusion"][0] == "=":
                if edge["fusion"][1] != "NO":
                    if "@" not in smarts_graph_list[0]["edges"][idx]['type'] and\
                        not edge["fusion"][1] == "EITHER":
                        smarts_graph_list[0]["edges"][idx]['and_type'] = "@"
                    smarts_graph_list[0]["nodes"][edge["start"]]["prop"]\
                        .append(keywords["fusion"][edge["fusion"][1]])
                    smarts_graph_list[0]["nodes"][edge["end"]]["prop"]\
                        .append(keywords["fusion"][edge["fusion"][1]])
            # negate "NO" is handled as "YES"
            elif edge["fusion"][0] == "#":
                if edge["fusion"][1] == "NO":
                    if "@" not in smarts_graph_list[0]["edges"][idx]['type']:
                        smarts_graph_list[0]["edges"][idx]['and_type'] = "@"
                    smarts_graph_list[0]["nodes"][edge["start"]]["prop"]\
                        .append(keywords["fusion"]["YES"])
                    smarts_graph_list[0]["nodes"][edge["end"]]["prop"]\
                        .append(keywords["fusion"]["YES"])
                    #other negate fusions are handled at the end
            else:
                raise SyntaxError("Something went wrong...")

    if not product:
        # !!!WARNING: handling negate fusions have to happend at the and
        # because of creating copys of the smartsgraph until now. (except remove duplicates)
        for edge_idx, edge in enumerate(patran_graph["edges"]):
            if edge["fusion"]:
                if edge["fusion"][0] == "#" and edge["fusion"][1] != "NO":
                    # create a deep copy of the current smarts_graph_list
                    # because there are two versions possible when handling negate fusion
                    smarts_graph_list += copy.deepcopy(smarts_graph_list)
                    for idx in range(len(smarts_graph_list)):
                        if idx < len(smarts_graph_list) / 2:
                            smarts_graph_list[idx]["nodes"][edge["start"]]["prop"]\
                                .append(keywords["no_fusion"][edge["fusion"][1]])
                        else:
                            smarts_graph_list[idx]["nodes"][edge["end"]]["prop"]\
                                .append(keywords["no_fusion"][edge["fusion"][1]])
                # "NO" is handles as negate "YES"
                if edge["fusion"][0] == "=" and edge["fusion"][1] == "NO":
                    # create a deep copy of the current smarts_graph_list
                    # because there are two versions possible when handling negate fusion
                    smarts_graph_list += copy.deepcopy(smarts_graph_list)
                    for idx in range(len(smarts_graph_list)):
                        if idx < len(smarts_graph_list) / 2:
                            smarts_graph_list[idx]["nodes"][edge["start"]]["prop"]\
                                .append(keywords["no_fusion"]["YES"])
                        else:
                            smarts_graph_list[idx]["nodes"][edge["end"]]["prop"]\
                                .append(keywords["no_fusion"]["YES"])

    # remove duplicates and empty strings in props
    for idx, smarts_graph in enumerate(smarts_graph_list):
        for node_idx, node in enumerate(smarts_graph["nodes"]):
            if node["prop"]:
                smarts_graph["nodes"][node_idx]["prop"] = list(filter(bool, list(set(node["prop"]))))

    return smarts_graph_list


def smarts_graph_to_string(smarts_graphs_list, product=False, add_atomindices=False, used_atom_indices=None):
    smarts_strings = []
    atomindices = []
    for pattern_idx, smarts_graph in enumerate(smarts_graphs_list):
        atoms = [None for _ in smarts_graph["nodes"]]
        for idx, node in enumerate(smarts_graph["nodes"]):
            # join atom types with 'or' (',') and properties with 'and' (';')
            atom_string = ",".join(node["atom"])
            if node["prop"]:
                atom_string += ";" + ";".join(node["prop"])
            # add atom map index
            if node["atom_map"] and not product:
                if add_atomindices and pattern_idx == 0:
                    atomindices.append(node["atom_map"])
                atom_string += ":" + str(node["atom_map"])
            elif product:
                # in CHMTRN/PATRAN product atom map indices are not given explicitly, but by order
                if used_atom_indices and idx+1 not in used_atom_indices:
                    pass
                else:
                    if add_atomindices and pattern_idx == 0:
                        atomindices.append(idx+1)
                    atom_string += ":" + str(idx+1)

            atoms[idx] = atom_string

        # finish atoms with square brackets '[...]'
        for idx, atom in enumerate(atoms):
            atoms[idx] = "[" + atom + "]"

        # add bonds
        ring_closure = 1 # counter for closing rings
        for bond in smarts_graph["edges"]:
            bond_num = len(bond["type"])
            if "and_type" in bond:
                bond_num += len(bond["and_type"])
            if product and bond_num > 1:
                bond_string = "~"
            else:
                bond_string = ",".join(bond["type"])
                if "and_type" in bond:
                    bond_string += ";" + ";".join(bond["and_type"])
            if bond["end"] > bond["start"]:
                atoms[bond["end"]] = bond_string + atoms[bond["end"]]
            else:
                atoms[bond["start"]] += bond_string + str(ring_closure)
                atoms[bond["end"]] += bond_string + str(ring_closure)
                ring_closure += 1

        # handle sidechains
        sidechain = [node["sidechain"] for node in smarts_graph["nodes"]]
        start_end = [[None, None] for _ in range(max(sidechain))]
        for idx, elem in enumerate(sidechain):
            if elem != 0:
                # handle new start of sidechain
                if start_end[elem-1][0] == None:
                    start_end[elem-1][0] = idx
                # handle end of sidechain
                start_end[elem-1][1] = idx
        start = [s[0] for s in start_end]
        end = [e[1] for e in start_end]
        for i in range(len(atoms)):
            if i in start:
                atoms[i] = "(" + atoms[i]
            if i in end:
                atoms[i] += ")"

        smarts_strings.append("".join(atoms))

    return smarts_strings, atomindices


def patran_to_smarts(patran, product=False, add_atomindices=False, used_atom_indices=None):
    parsed_patran = check_patran_syntax(patran)
    patran_graph = patran_to_graph(parsed_patran)
    smarts_graphs_list =  graph_patran_to_smarts(patran_graph, product=product)
    return smarts_graph_to_string(smarts_graphs_list, product=product,
                                    add_atomindices=add_atomindices, used_atom_indices=used_atom_indices)


def patran_to_reactionsmarts(patran_string, validate_atomindices=False):
    patran_string = patran_string.replace(" ", "")
    split = patran_string.split("=>")
    if len(split) != 2:
        raise SyntaxError("invalid patran reaction pattern; "\
            "there must be a product followed by a retro-reaction arrow '=>' "\
            "and at least one reactant. Multiple reactants are separated by plus-signs '+'")
    products, reactants = split
    reactants_smarts, reactant_atomindices = patran_to_smarts(reactants,
                    product=False,
                    add_atomindices=True)

    products_smarts, products_atomindices = patran_to_smarts(products,
                    product=False,
                    add_atomindices=True,
                    used_atom_indices=set(reactant_atomindices))

    # validation of atomindices in reactants and product
    if validate_atomindices:
        reactant_atomindices_set = set(reactant_atomindices)
        if len(reactant_atomindices_set) != len(reactant_atomindices):
            raise ValueError("There are duplicates of atom indices in the reactants")
        product_atomindices_set = set(products_atomindices)
        diff = reactant_atomindices_set.symmetric_difference(product_atomindices_set)
        if diff:
            raise SyntaxError(f"There are different atomindces in reactants: "\
            f"({', '.join([str(x) for x in reactant_atomindices_set-product_atomindices_set])}) "\
                f"and product: ({', '.join([str(x) for x in product_atomindices_set-reactant_atomindices_set])})" \
                    f"pattern: {', '.join([str(x) for x in diff])}")

    return [smarts + ">>" + products_smarts[0] for smarts in reactants_smarts]


def translate(patran_string, validate_atomindices=False, pbar=None):
    '''Transpile PATRAN/CHMTRN reaction pattern to SMIRKS and PATRAN/CHMTRN pattern to SMARTS'''
    try:
        if patran_string.count("=>") == 1:
            result = patran_to_reactionsmarts(patran_string, validate_atomindices=validate_atomindices)

        elif patran_string.count("=>") == 0:
            result = patran_to_smarts(patran_string)

        else:
            raise SyntaxError("invalid patran reaction pattern; "\
                "there must be a product followed by a retro-reaction arrow '=>' "\
                "and at least one reactant. Multiple reactants are separated by plus-signs '+'")

        return result

    except Exception as e:
        if pbar:
            pbar.write(f"Error in pattern {patran_string}: {e}")
        else:
            print(f"Error in pattern {patran_string}: {e}")
        return []



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Transpile PATRAN/CHMTRN reaction pattern to SMIRKS')
    parser.add_argument('-f', '--file', type=str, help='PATRAN/CHMTRN reaction pattern json file to translate')
    parser.add_argument('-o', '--output', type=str, help='output JSON file to save the translated SMIRKS')
    parser.add_argument('reaction_pattern', type=str, nargs='?', default=None,
                            help='PATRAN/CHMTRN reaction pattern string to translate')

    args = parser.parse_args()

    if args.file and args.output:
        with open(args.file) as f:
            data = json.load(f)
        result = {"_comment": ["This material is part of the SAVI-Space created by Malte Korn and Matthias Rarey",
                "at the Center for Bioinformatics, University of Hamburg, with support from",
                "Marc Nicklaus (NIH,NCI), Phil Judson, Raphael Klein (BioSolveIT GmbH) and",
                "Christian Lemmen (BioSolveIT GmbH).",
                "Philip Judson and Marc Nicklaus provided the LHASA transform rules and assisted with their knowledge about SAVI-2020.",
                "SAVI-Space and all its components including this file/directory is licensed under CC-BY-NC 4.0."]}
        pbar = tqdm(data, desc="Transpiling patran to smarts")
        for key in pbar:
            if key.startswith("_"):
                continue
            number, name = key.split("_", 1)
            result[number] = {"name": name, "unique_matching": [], "ring": None, "smirks": []}
            for reaction in tqdm(data[key], leave=False, desc=f"Transpiling {key}"):

                result[number]["smirks"].append(translate(reaction, pbar=pbar))
                result[number]["unique_matching"].append([[] for _ in range(len(result[number]["smirks"][-1]))])

            # validate smirks with rdkit
            for smirks in result[number]["smirks"]:
                for smarts in smirks:
                    if not Chem.MolFromSmarts(smarts) and not AllChem.ReactionFromSmarts(smarts):
                        print(Chem.MolFromSmarts(smarts), AllChem.ReactionFromSmarts(smarts))
                        pbar.write(f"Invalid smirks: {smarts}\n{key}")

        with open(args.output, 'w') as f:
            json.dump(result, f, indent=4)

    elif args.reaction_pattern and not args.file and not args.output:
        result = translate(args.reaction_pattern)
        for smarts in result:
            for s in smarts:
                print(s)

    else:
        raise ValueError("Either provide a patran string:\n" +
                         "\t<patran_string>\n" +\
                         "or a file and an output file:\n" +\
                            "\t-f <patran_file> -o <output_file>\n" +\
                            "-h for help")
