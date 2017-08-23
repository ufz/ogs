# Setting (local) git include path to .gitconfig
execute_process(
    COMMAND ${GIT_TOOL_PATH} config --local include.path ../.gitconfig
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)