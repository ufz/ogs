/*
 * GMSHInterface.h
 *
 *  Created on: Apr 29, 2010
 *      Author: TF
 */

#ifndef GMSHINTERFACE_H_
#define GMSHINTERFACE_H_

#include <string>
#include <list>

// FileIO
#include "Writer.h"
#include "GMSHPoint.h"
#include "GMSHPolygonTree.h"
#include "GMSHMeshDensityStrategy.h"

// GEOLIB
#include "GEOObjects.h"
#include "Polygon.h"

namespace MeshLib
{
class CFEMesh;
}

namespace FileIO
{

namespace GMSH {

enum MeshDensityAlgorithm {
	NoMeshDensity = 0, //!< do not set the parameter
	FixedMeshDensity, //!< set the parameter with a fixed value
	AdaptiveMeshDensity //!< computing the mesh density employing a QuadTree
};

}

/**
 * \brief Reads and writes GMSH-files to and from OGS data structures.
 */
class GMSHInterface : public Writer
{
public:

	/**
	 *
	 * @param geo_objs reference tp instance of class GEOObject that maintains the geometries.
	 * 	The instance is used for preparation geometries for writing them to the gmsh file format.
	 * @param include_stations_as_constraints switch to enable writing stations as constraints
	 * @param mesh_density_algorithm one of the mesh density algorithms (\@see enum MeshDensityAlgorithm)
	 * @param param1 parameter that can be used for the mesh density algorithm
	 * @param param2 parameter that can be used for the mesh density algorithm
	 * @param param3 parameter that can be used for the mesh density algorithm
	 * @param selected_geometries vector of names of geometries, that should be employed for mesh generation.
	 * @return
	 */
	GMSHInterface (GEOLIB::GEOObjects & geo_objs,
					bool include_stations_as_constraints,
					GMSH::MeshDensityAlgorithm mesh_density_algorithm,
					double param1, double param2, size_t param3,
					std::vector<std::string> & selected_geometries);

	/**
	 * checks if there is a GMSH mesh file header
	 * @param fname the file name of the mesh (including the path)
	 * @return true, if the file seems to be a valid GMSH file, else false
	 */
	static bool isGMSHMeshFile (const std::string& fname);
	/**
	 * reads a mesh created by GMSH - this implementation is based on the former function GMSH2MSH
	 * @param fname the file name of the mesh (including the path)
	 * @param mesh the new mesh
	 * @return
	 */
	static void readGMSHMesh (std::string const& fname, MeshLib::CFEMesh* mesh);

protected:
	int write(std::ostream& stream);

private:
	/**
	 * 1. get and merge data from _geo_objs
	 * 2. compute topological hierarchy
	 * @param out
	 */
	void writeGMSHInputFile(std::ostream & out);

	void writePoints(std::ostream& out) const;

	size_t _n_lines;
	size_t _n_plane_sfc;

	GEOLIB::GEOObjects & _geo_objs;
	std::vector<std::string>& _selected_geometries;
	std::string _gmsh_geo_name;
	std::list<GMSHPolygonTree*> _polygon_tree_list;

	bool _include_stations_as_constraints;

	std::vector<FileIO::GMSHPoint*> _gmsh_pnts;

	GMSHMeshDensityStrategy *_mesh_density_strategy;
};
}

#endif /* GMSHINTERFACE_H_ */
