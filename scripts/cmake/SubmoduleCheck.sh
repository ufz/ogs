# Checks for unitialized git submodules
# Parameter 1: The path to the submodule
# Returns 0 if there is no such submodule
# Returns 1 if the module is already initialized
result=$(git submodule status $1 | grep '^-' | wc -l); exit $result