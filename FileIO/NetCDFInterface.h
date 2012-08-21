/**
 * \file NetCDFInterface.h
 * 29/07/2010 YW Initial implementation
 */

#ifndef NetCDFInterface_H
#define NetCDFInterface_H

#include "GEOObjects.h"
#include <string>
#include <vector>

namespace GEOLIB
{
class GEOObjects;
}

namespace MeshLib
{
class Mesh;
}

namespace FileIO
{
class NetCDFInterface
{
public:
	/// Import climate data from a NetCDF file.
	static void readNetCDFData(std::string &fname,
	                           std::vector<GeoLib::Point*>* points_vec,
	                           GEOLIB::GEOObjects* geo_obj,
	                           size_t &NRLAT,
	                           size_t &NRLON);
	/// Convert imported point data to an CFEMesh.
	static MeshLib::Mesh* createMeshFromPoints(std::vector<GeoLib::Point*>* points_vec,
	                                              size_t &NRLAT,
	                                              size_t &NRLON);

private:
};
} // end namespace

#endif /* NetCDFInterface_H */
