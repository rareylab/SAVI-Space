#!/bin/bash

echo 'Creating Space...'

# check if the variables already set in the environment
if [ -z "$SCRIPT_GENERATE_SPACELIGHT_DB_FILE" ]; then
echo "Setting up the environment variables..."
export SCRIPT_GENERATE_SPACELIGHT_DB_FILE=/work/korn/SAVISpace/bin/BSIT/script/generate_spacelight_db_file.py
export SCRIPT_GENERATE_SINGLE_FRAGMENT_SPACES=/work/korn/SAVISpace/bin/BSIT/script/generate_single_fragment_spaces.py
export SCRIPT_MERGE_SINGLE_FRAGMENT_SPACES=/work/korn/SAVISpace/bin/BSIT/script/merge_single_fragment_spaces.py

export TOOL_SPACELIGHT_DB_CREATOR=/work/korn/SAVISpace/bin/BSIT/spacelight_db_creator
export TOOL_EXECUTABLE_REACTION_SYNTHESIZER=/work/korn/SAVISpace/bin/BSIT/reaction_synthesizer
export TOOL_FRAGSPACE_MERGER=/work/korn/SAVISpace/bin/BSIT/fragspace_merger

export FUNCTIONAL_GROUP_LABEL=/work/korn/SAVISpace/bin/BSIT/FunctionalGroupLabel.txt
export PROTECTING_GROUPS=/work/korn/SAVISpace/bin/BSIT/ProtectingGroups.txt
fi
SCRIPT_GENERATE_SPACELIGHT_DB_FILE=/work/korn/SAVISpace/bin/BSIT/script/generate_spacelight_db_file.py
SCRIPT_GENERATE_SINGLE_FRAGMENT_SPACES=/work/korn/SAVISpace/bin/BSIT/script/generate_single_fragment_spaces.py
SCRIPT_MERGE_SINGLE_FRAGMENT_SPACES=/work/korn/SAVISpace/bin/BSIT/script/merge_single_fragment_spaces.py

TOOL_SPACELIGHT_DB_CREATOR=/work/korn/SAVISpace/bin/BSIT/spacelight_db_creator
TOOL_EXECUTABLE_REACTION_SYNTHESIZER=/work/korn/SAVISpace/bin/BSIT/reaction_synthesizer
TOOL_FRAGSPACE_MERGER=/work/korn/SAVISpace/bin/BSIT/fragspace_merger

FUNCTIONAL_GROUP_LABEL=/work/korn/SAVISpace/bin/BSIT/FunctionalGroupLabel.txt
PROTECTING_GROUPS=/work/korn/SAVISpace/bin/BSIT/ProtectingGroups.txt

# python $SCRIPT_GENERATE_SPACELIGHT_DB_FILE -t $TOOL_SPACELIGHT_DB_CREATOR -i /work/korn/SAVI-Space/SAVI-Space-2020/space_tool_input/Colibri/SAVI-Space-2020.list -o /work/korn/SAVI-Space/SAVI-Space-2020/spaces/Colibri -f SAVI-Space-2020 --tool_executable_reaction_synthesizer $TOOL_EXECUTABLE_REACTION_SYNTHESIZER -s $FUNCTIONAL_GROUP_LABEL -p $PROTECTING_GROUPS

# python $SCRIPT_GENERATE_SINGLE_FRAGMENT_SPACES -t $TOOL_EXECUTABLE_REACTION_SYNTHESIZER -i /work/korn/SAVI-Space/SAVI-Space-2020/space_tool_input/Colibri/SAVI-Space-2020.list -o /work/korn/SAVI-Space/SAVI-Space-2020/spaces/Colibri/singleFragSpaces -s $FUNCTIONAL_GROUP_LABEL -p $PROTECTING_GROUPS

python $SCRIPT_MERGE_SINGLE_FRAGMENT_SPACES -t $TOOL_FRAGSPACE_MERGER -i /work/korn/SAVI-Space/SAVI-Space-2020/space_tool_input/Colibri/SAVI-Space-2020.list -d /work/korn/SAVI-Space/SAVI-Space-2020/spaces/Colibri/singleFragSpaces -o /work/korn/SAVI-Space/SAVI-Space-2020/spaces/Colibri/mergedFragSpace -f Merged_SAVI-Space-2020 -a /work/korn/SAVI-Space/SAVI-Space-2020/spaces/Colibri/SAVI-Space-2020.tfsdb