/*
 * Gmsh2GeoIO.h
 *
 *  Created on: Aug 18, 2011
 *      Author: TF
 */

#ifndef GMSH2GEOIO_H_
#define GMSH2GEOIO_H_

#include <string>

namespace GeoLib
{
class GEOObjects;
}

namespace FileIO
{
class Gmsh2GeoIO
{
public:
	/**
	 * load a surface mesh (gmsh format) as a geometric surface
	 * @param fname file name of the surface mesh
	 * @param geo the object that manages all geometries,
	 * new surface will be put into this container
	 */
	static void loadMeshAsGeometry (std::string & fname, GeoLib::GEOObjects* geo);
};
} // end namespace FileIO

#endif /* GMSH2GEOIO_H_ */
