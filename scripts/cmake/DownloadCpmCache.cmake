option(OGS_DOWNLOAD_CPM_CACHE "Download cpm cache into build directory" OFF)
if(NOT OGS_DOWNLOAD_CPM_CACHE)
    return()
endif()

file(READ ${PROJECT_SOURCE_DIR}/web/data/versions.json jsonFileString)
string(JSON _cpm_package_id GET "${jsonFileString}" cpm package_id)
string(JSON _cpm_package_hash GET "${jsonFileString}" cpm sha256)
if(NOT EXISTS ${PROJECT_BINARY_DIR}/cpm)
    message(STATUS "Downloading cpm cache with file id ${_cpm_package_id}")
    file(
        DOWNLOAD
        https://gitlab.opengeosys.org/ogs/ogs/-/package_files/${_cpm_package_id}/download
        cpm.tar.gz
        SHOW_PROGRESS
        EXPECTED_HASH SHA256=${_cpm_package_hash}
    )
    file(ARCHIVE_EXTRACT INPUT cpm.tar.gz)
endif()
set(ENV{CPM_SOURCE_CACHE} ${PROJECT_BINARY_DIR}/cpm)
