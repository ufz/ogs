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
 * \file Gmsh2GeoIO.h
 *
 *  Created on 2011-08-18 by Thomas Fischer
 */

#ifndef GMSH2GEOIO_H_
#define GMSH2GEOIO_H_

#include <string>

namespace GeoLib
{
class GEOObjects;
}

namespace MeshLib
{
	class Mesh;
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

	/**
	 * Converts a 2D mesh into a geometry.
	 * A new geometry with the name of the mesh will be inserted into geo_objects, consisting
	 * of points identical with mesh nodes and one surface representing the mesh. Triangles are
	 * converted to geometric triangles, quads are split into two triangles, all other elements
	 * are ignored.
	 */
	static bool convertMeshToGeo(const MeshLib::Mesh &mesh, GeoLib::GEOObjects* geo_objects);

};
} // end namespace FileIO

#endif /* GMSH2GEOIO_H_ */
