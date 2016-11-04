# Check requirements / supported configurations
if(MSVC AND NOT HAVE_64_BIT AND NOT OGS_32_BIT)
    message(FATAL_ERROR "Building OGS on Windows with 32-bit is not supported! \
Either use the correct generator, e.g. 'Visual Studio 14 2015 Win64' or define \
'-DOGS_32_BIT=ON' if you know what you are doing.")
endif()

# Set build directories
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

# Enable Visual Studio project folder grouping
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

include(GetGitRevisionDescription)
GET_GIT_HEAD_REVISION(GIT_REFSPEC GIT_SHA1)
string(SUBSTRING ${GIT_SHA1} 0 8 GIT_SHA1_SHORT)

if(IS_SUBPROJECT)
    set(OGS_VERSION x.x.x)
else()
    GIT_GET_TAG(GIT_DESCRIBE)
    if(GIT_DESCRIBE)
        string(REGEX MATCH ^[0-9|\\.]+ GIT_TAG ${GIT_DESCRIBE})
        set(OGS_VERSION ${GIT_DESCRIBE})

        if(GIT_DESCRIBE MATCHES ".*-.*-.*")
            # Commit is not a tag
            string(REGEX MATCH "-([0-9]+)-" GIT_COMMITS_AFTER_TAG ${GIT_DESCRIBE})
        else()
            set(OGS_VERSION ${GIT_TAG})
        endif()
        message(STATUS "OGS version: ${OGS_VERSION}")
    else()
        message(WARNING "Git repository contains no tags! Please run: git fetch --tags")
    endif()
endif()
