set(OGS_BOOST_VERSION 1.56.0)
string(REPLACE "." "_" OGS_BOOST_VERSION_UNDERLINE ${OGS_BOOST_VERSION})
set(OGS_BOOST_URL "http://opengeosys.s3.amazonaws.com/ogs6-lib-sources/boost_${OGS_BOOST_VERSION_UNDERLINE}.tar.bz2")
set(OGS_BOOST_MD5 "a744cf167b05d72335f27c88115f211d")

set(OGS_EIGEN_URL "http://opengeosys.s3.amazonaws.com/ogs6-lib-sources/eigen-3.2.8.tar.gz")
set(OGS_EIGEN_MD5 "135d8d43aaee5fb54cf5f3e981b1a6db")

set(OGS_VTK_VERSION 7.0.0)
set(OGS_VTK_URL "http://www.vtk.org/files/release/7.0/VTK-7.0.0.tar.gz")
set(OGS_VTK_MD5 "5fe35312db5fb2341139b8e4955c367d")

set(OGS_TIFF_URL "https://github.com/ufz/tiff/archive/4.0.1.zip")
set(OGS_TIFF_MD5 "8d5c18654bda9c731d8a6e2dc958751b")
set(OGS_GEOTIFF_URL "https://github.com/ufz/geotiff/archive/1.4.0.zip")
set(OGS_GEOTIFF_MD5 "5c86853b96f0edca7b9d9e894a8dd93b")
