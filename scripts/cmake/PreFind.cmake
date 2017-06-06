### Find Git
# Check for cmder git installed via chocolatey
find_program(GIT_EXECUTABLE
  NAMES git
  PATHS C:/tools/cmder/vendor/git-for-windows
  PATH_SUFFIXES cmd bin
  DOC "Git command line client"
)

find_package(Git REQUIRED)
string(REPLACE "mingw64/" "" GIT_EXECUTABLE ${GIT_EXECUTABLE}) # Windows git submodule fix
set(GIT_TOOL_PATH ${GIT_EXECUTABLE} CACHE FILEPATH "The git command line interface" FORCE)
### End Find Git
