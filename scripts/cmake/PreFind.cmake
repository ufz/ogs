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
### End Find Git

execute_process(COMMAND ${GIT_EXECUTABLE} status
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    RESULT_VARIABLE IS_GIT_REPO
    OUTPUT_QUIET)
if(IS_GIT_REPO GREATER 0)
    set(IS_GIT_REPO FALSE CACHE INTERNAL "")
    if(DEFINED OGS_VERSION)
        message(WARNING "Using user-provided OGS_VERSION: ${OGS_VERSION}!")
        message(WARNING "Submodule setup is skipped!")
    else()
        message(FATAL_ERROR "No git repository found at ${PROJECT_SOURCE_DIR}! "
            "Please use git to obtain the source code! See "
            "https://www.opengeosys.org/docs/devguide/getting-started/get-the-source-code/"
            " OR manually set the OGS_VERSION variable.")
    endif()
else()
set(IS_GIT_REPO TRUE CACHE INTERNAL "")
endif()
