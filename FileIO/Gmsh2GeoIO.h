/**
 * \file
 * \author Thomas Fischer
 * \date   2011-08-18
 * \brief  Definition of the Gmsh2GeoIO class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
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
