set(BASE_URL "http://ogsstorage.blob.core.windows.net/jenkins/ogs6-lib-sources")

set(OGS_BOOST_VERSION 1.56.0)
string(REPLACE "." "_" OGS_BOOST_VERSION_UNDERLINE ${OGS_BOOST_VERSION})
set(OGS_BOOST_URL "${BASE_URL}/boost_${OGS_BOOST_VERSION_UNDERLINE}.tar.bz2")
set(OGS_BOOST_MD5 "a744cf167b05d72335f27c88115f211d")

set(OGS_EIGEN_URL "${BASE_URL}/eigen-3.2.9.tar.gz")
set(OGS_EIGEN_MD5 "6a578dba42d1c578d531ab5b6fa3f741")

set(OGS_VTK_VERSION 6.3.0)
set(OGS_VTK_URL "${BASE_URL}/vtk-6.3.0.tar.gz")
set(OGS_VTK_MD5 "0231ca4840408e9dd60af48b314c5b6d")

set(OGS_TIFF_URL "${BASE_URL}/tiff-4.0.1.zip")
set(OGS_TIFF_MD5 "d4dc85296f790c159fc8fab443dd697c")
set(OGS_GEOTIFF_URL "${BASE_URL}/libgeotiff-1.4.0.zip")
set(OGS_GEOTIFF_MD5 "6f2583a0ec88eeb589f1cd326b2ed49c")
