option(OGS_DOWNLOAD_CPM_CACHE "Download cpm cache into build directory" OFF)
if(NOT OGS_DOWNLOAD_CPM_CACHE)
    return()
endif()

file(READ ${PROJECT_SOURCE_DIR}/web/data/versions.json jsonFileString)
string(JSON _cpm_package_file_id GET "${jsonFileString}" cpm package_file_id)
string(JSON _cpm_package_file_hash GET "${jsonFileString}" cpm
       package_file_sha256
)
set(_cpm_cache_dir ${PROJECT_BINARY_DIR}/.cpm-${_cpm_package_file_id})
if(NOT EXISTS ${_cpm_cache_dir})
    message(STATUS "Downloading cpm cache with file id ${_cpm_package_file_id}")
    file(
        DOWNLOAD
        https://gitlab.opengeosys.org/ogs/ogs/-/package_files/${_cpm_package_file_id}/download
        cpm.tar.gz
        SHOW_PROGRESS
        EXPECTED_HASH SHA256=${_cpm_package_file_hash}
    )
    file(ARCHIVE_EXTRACT INPUT cpm.tar.gz)
    file(RENAME ${PROJECT_BINARY_DIR}/cpm ${_cpm_cache_dir})
endif()
set(ENV{CPM_SOURCE_CACHE} ${_cpm_cache_dir})
