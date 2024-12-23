echo "SAVI-Space-2024"
echo "reactants validation"
# python src/space_builder/reactants_validation.py -n SAVI-Space-2024 --reaction_path smirks_Naomi.json -b /work/korn/SAVISpace/data/building_blocks/Enamine_Building_Blocks_Stock_288748cmpd_20240717_max_molecular_weight_700.0_fragments_standardized.smi --apply_kill --n_cores 36 --unique_matches --cluster  --overwrite
echo "space creation"
python src/space_builder/naomi_space_creation.py -n SAVI-Space-2024 --n_cores 36 --reaction_path smirks_Naomi.json --deprotection_path /work/korn/SAVISpace/data/space_creation_data/protectionpattern.smarts --overwrite
python src/space_builder/Colibri_space_creation.py -n SAVI-Space-2024 --n_cores 36 --reaction_path smirks_Colibri.json --deprotection_path /work/korn/SAVISpace/data/space_creation_data/protectionpattern.smarts --overwrite

echo "SAVI-Space-2020"
echo "reactants validation"
# python src/space_builder/reactants_validation.py -n SAVI-Space-2020 --reaction_path smirks_Naomi.json -b /work/korn/SAVISpace/data/building_blocks/Dec2019_instock_BBs_155k-sdf_max_molecular_weight_700.0_fragments_standardized.smi --apply_kill --n_cores 36 --unique_matches --cluster --overwrite
echo "space creation"
python src/space_builder/naomi_space_creation.py -n SAVI-Space-2020 --n_cores 36 --reaction_path smirks_Naomi.json --deprotection_path /work/korn/SAVISpace/data/space_creation_data/protectionpattern.smarts --overwrite
python src/space_builder/Colibri_space_creation.py -n SAVI-Space-2020 --n_cores 36 --reaction_path smirks_Colibri.json --deprotection_path /work/korn/SAVISpace/data/space_creation_data/protectionpattern.smarts --overwrite

echo "SAVI-Space-2020-Librules"
echo "reactants validation"
# python src/space_builder/reactants_validation.py -n SAVI-Space-2020-Librules -b /work/korn/SAVISpace/data/building_blocks/Dec2019_instock_BBs_155k-sdf_max_molecular_weight_700.0_fragments_standardized.smi --apply_kill --n_cores 36 --uniquify  --disallow_exocyclic_doublebonds --overwrite
echo "space creation"
python src/space_builder/naomi_space_creation.py -n SAVI-Space-2020-Librules --n_cores 36 --deprotection_path /work/korn/SAVISpace/data/space_creation_data/protectionpattern.smarts -c --overwrite
