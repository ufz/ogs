# Prepend directory to each expected file
set(EXPECTED_FILES_ABS)
foreach(file IN LISTS EXPECTED_FILES)
    list(APPEND EXPECTED_FILES_ABS ${DIR_TO_CHECK}/${file})
endforeach()

# Collect all files starting with the prefix
file(GLOB FILES_CREATED "${DIR_TO_CHECK}/${FILE_PREFIX}*")

# Sort both lists
list(SORT EXPECTED_FILES_ABS)
list(SORT FILES_CREATED)

# Check that only expected files are created.
if(NOT FILES_CREATED STREQUAL EXPECTED_FILES_ABS)
    message(WARNING "Expected only the following files with prefix ${FILE_PREFIX} to be created in the ${DIR_TO_CHECK} directory:")
    message(WARNING "${EXPECTED_FILES_ABS}")
    message(WARNING "but other files were created:")
    message(WARNING "${FILES_CREATED}")
    message(FATAL_ERROR "Unexpected files created.")
endif()
