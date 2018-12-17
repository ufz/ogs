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

execute_process(COMMAND ${GIT_EXECUTABLE} status
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    RESULT_VARIABLE IS_GIT_REPO
    OUTPUT_QUIET)
if(IS_GIT_REPO GREATER 0)
    message(FATAL_ERROR "Source code at ${PROJECT_SOURCE_DIR} is not a git "
        "repository! Please use git to obtain the source code! See "
        "https://www.opengeosys.org/docs/devguide/getting-started/get-the-source-code/")
endif()
