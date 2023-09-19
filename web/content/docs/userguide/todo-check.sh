#!/bin/bash

# This script greps all occurances of "*TODO:" in md files in subfolders (recursively)

# Can be used to update list of missing content in the userguide

# Run it from the folder where it is stored with sh todo-check.sh

echo "## TODOs in userguide/basics\n" > todo.md

grep -R -B 2 -n 'TODO:' basics/ >> todo.md

#echo -e "\n" >> todo.md

echo "## TODOs in userguide/blocks\n" >> todo.md

grep -R -B 2 -n 'TODO:' blocks/ >> todo.md

#echo -e "\n" >> todo.md

echo "## TODOs in userguide/features\n" >> todo.md

grep -R -B 2 -n 'TODO:' features/ >> todo.md

#echo -e "\n" >> todo.md

echo "## TODOs in userguide/troubleshooting\n" >> todo.md

grep -R -B 2 -n 'TODO:' troubleshooting/ >> todo.md

#echo -e "\n" >> todo.md
