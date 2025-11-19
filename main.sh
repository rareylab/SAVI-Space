
# This script is used to create the SAVI-Space for the two datasets SAVI-Space-2020 and SAVI-Space-2024.
# The script is divided into two parts: one for the creation of the space without applying the kill and one for the creation of the space with the kill.

export SCRIPT_GENERATE_SPACELIGHT_DB_FILE=./generate_spacelight_db_file.py
export SCRIPT_GENERATE_SINGLE_FRAGMENT_SPACES=./generate_single_fragment_spaces.py
export SCRIPT_MERGE_SINGLE_FRAGMENT_SPACES=./merge_single_fragment_spaces.py

export TOOL_SPACELIGHT_DB_CREATOR=./spacelight_db_creator
export TOOL_EXECUTABLE_REACTION_SYNTHESIZER=./reaction_synthesizer
export TOOL_FRAGSPACE_MERGER=./fragspace_merger

export FUNCTIONAL_GROUP_LABEL=./FunctionalGroupLabel.txt
export PROTECTING_GROUPS=./ProtectingGroups.txt

export SPACELIGHT=./spacelight

############################################################################################################################################################################

# without kill

#################################
echo "SAVI-Space-2024 no_kill" ##
#################################

echo "reactants validation"
###########################
python src/space_builder/reactants_validation.py -n SAVI-Space-2024-no_kill  --reaction_path smirks_Colibri.json  -b /work/korn/SAVISpace/data/building_blocks/Enamine_Building_Blocks_Stock_288748cmpd_20240717_max_molecular_weight_700.0_fragments_standardized.smi --unique_matches --cluster  --overwrite

echo "space creation"
#####################

echo "NAOMI"
python src/space_builder/naomi_space_creation.py -n SAVI-Space-2024-no_kill --no_kills --n_cores 36 --reaction_path smirks_Naomi.json --deprotection_path /work/korn/SAVISpace/data/space_creation_data/protectionpattern.smarts --overwrite
bash SAVI-Space-2024-no_kill/create_naomi_space.sh

echo "Colibri"
python src/space_builder/Colibri_space_creation.py -n SAVI-Space-2024-no_kill --no_kills --n_cores 36 --reaction_path smirks_Colibri.json --deprotection_path /work/korn/SAVISpace/data/space_creation_data/protectionpattern.smarts --overwrite
bash SAVI-Space-2024-no_kill/create_Colibri_space.sh

#################################
echo "SAVI-Space-2020 no_kill" ##
#################################

echo "reactants validation"
###########################
python src/space_builder/reactants_validation.py -n SAVI-Space-2020-no_kill --reaction_path smirks_Colibri.json -b /work/korn/SAVISpace/data/building_blocks/Dec2019_instock_BBs_155k-sdf_max_molecular_weight_700.0_fragments_standardized.smi --unique_matches --cluster --overwrite

echo "space creation"
#####################

echo "NAOMI"
python src/space_builder/naomi_space_creation.py -n SAVI-Space-2020-no_kill --no_kills --n_cores 36 --reaction_path smirks_Naomi.json --deprotection_path /work/korn/SAVISpace/data/space_creation_data/protectionpattern.smarts --overwrite
bash SAVI-Space-2020-no_kill/create_naomi_space.sh

echo "Colibri"
python src/space_builder/Colibri_space_creation.py -n SAVI-Space-2020-no_kill --no_kills --n_cores 36 --reaction_path smirks_Colibri.json --deprotection_path /work/korn/SAVISpace/data/space_creation_data/protectionpattern.smarts --overwrite
bash SAVI-Space-2020-no_kill/create_Colibri_space.sh


############################################################################################################################################################################

# with kill

#########################
echo "SAVI-Space-2020" ##
#########################

echo "reactants validation"
###########################
python src/space_builder/reactants_validation.py -n SAVI-Space-2020 --reaction_path smirks_Colibri.json -b /work/korn/SAVISpace/data/building_blocks/Dec2019_instock_BBs_155k-sdf_max_molecular_weight_700.0_fragments_standardized.smi --apply_kill --unique_matches --cluster --overwrite

echo "space creation"
#####################

echo "NAOMI"
python src/space_builder/naomi_space_creation.py -n SAVI-Space-2020 --n_cores 36 --reaction_path smirks_Naomi.json --deprotection_path /work/korn/SAVISpace/data/space_creation_data/protectionpattern.smarts --overwrite
bash SAVI-Space-2020/create_naomi_space.sh

echo "Colibri"
python src/space_builder/Colibri_space_creation.py -n SAVI-Space-2020 --n_cores 36 --reaction_path smirks_Colibri.json --deprotection_path /work/korn/SAVISpace/data/space_creation_data/protectionpattern.smarts --overwrite
bash SAVI-Space-2020/create_Colibri_space.sh

##################################
echo "SAVI-Space-2020-Librules" ##
##################################

echo "reactants validation"
###########################
python src/space_builder/reactants_validation.py -n SAVI-Space-2020-Librules --b /work/korn/SAVISpace/data/building_blocks/Dec2019_instock_BBs_155k-sdf_max_molecular_weight_700.0_fragments_standardized.smi --apply_kill --uniquify  --disallow_exocyclic_doublebonds --overwrite

echo "space creation"
#####################

echo "NAOMI"
python src/space_builder/naomi_space_creation.py -n SAVI-Space-2020-Librules --n_cores 36 --deprotection_path /work/korn/SAVISpace/data/space_creation_data/protectionpattern.smarts -c --overwrite
bash SAVI-Space-2020-Librules/create_naomi_space.sh

#########################
echo "SAVI-Space-2024" ##
#########################

echo "reactants validation"
###########################
python src/space_builder/reactants_validation.py -n SAVI-Space-2024 --reaction_path smirks_Colibri.json  -b /work/korn/SAVISpace/data/building_blocks/Enamine_Building_Blocks_Stock_288748cmpd_20240717_max_molecular_weight_700.0_fragments_standardized.smi --apply_kill --unique_matches --cluster  --overwrite

echo "space creation"
#####################

echo "NAOMI"
python src/space_builder/naomi_space_creation.py -n SAVI-Space-2024 --n_cores 36 --reaction_path smirks_Naomi.json --deprotection_path /work/korn/SAVISpace/data/space_creation_data/protectionpattern.smarts --overwrite
bash SAVI-Space-2024/create_naomi_space.sh

echo "Colibri"
python src/space_builder/Colibri_space_creation.py -n SAVI-Space-2024 --n_cores 36 --reaction_path smirks_Colibri.json --deprotection_path /work/korn/SAVISpace/data/space_creation_data/protectionpattern.smarts --overwrite
bash SAVI-Space-2024/create_Colibri_space.sh
