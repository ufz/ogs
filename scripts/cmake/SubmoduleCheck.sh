# Checks for unitialized git submodules
# Parameter 1: The path to the submodule
# Returns 0 if the module is already initialized
# Returns 1 if there is no such submodule
# Returns 2 if the checked out submodule is out of date

# Check for modified
result=$(git submodule status $1 | grep '^+' | wc -l);
if [[ result -eq 0 ]]; then
	result=$(git submodule status $1 | grep '^-' | wc -l);
	exit $result;
else
	exit 2;
fi
