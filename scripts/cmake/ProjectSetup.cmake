# Check requirements / supported configurations
if(MSVC AND NOT HAVE_64_BIT AND NOT OGS_32_BIT)
    message(FATAL_ERROR "Building OGS on Windows with 32-bit is not supported! \
Either use the correct generator, e.g. 'Visual Studio 14 2015 Win64' or define \
'-DOGS_32_BIT=ON' if you know what you are doing.")
endif()

# Set build directories
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
if(USE_CONAN AND MSVC)
    foreach(OUTPUTCONFIG ${CMAKE_CONFIGURATION_TYPES})
        string(TOUPPER ${OUTPUTCONFIG} OUTPUTCONFIG)
        set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
        set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
        set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
    endforeach(OUTPUTCONFIG CMAKE_CONFIGURATION_TYPES)
endif()

set(Data_SOURCE_DIR ${PROJECT_SOURCE_DIR}/Tests/lfs-data CACHE INTERNAL "")
set(Data_BINARY_DIR ${PROJECT_BINARY_DIR}/Tests/lfs-data CACHE INTERNAL "")

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
