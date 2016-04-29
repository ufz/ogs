/**
 * \file
 * \author Thomas Fischer
 * \date   2010-04-29
 * \brief  Definition of the GMSHInterface class.
 *
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef GMSHINTERFACE_H_
#define GMSHINTERFACE_H_

#include <list>
#include <map>
#include <string>

#include "BaseLib/IO/Writer.h"
#include "FileIO/GmshIO/GMSHPoint.h"
#include "FileIO/GmshIO/GMSHPolygonTree.h"
#include "FileIO/GmshIO/GMSHMeshDensityStrategy.h"

#include "MathLib/LinAlg/Dense/DenseMatrix.h"

namespace GeoLib
{
	class GEOObjects;
	class Polygon;
}

namespace MeshLib
{
	class Mesh;
	class Element;
	class Node;
}

namespace FileIO
{

namespace GMSH {

enum class MeshDensityAlgorithm {
	FixedMeshDensity, //!< set the parameter with a fixed value
	AdaptiveMeshDensity //!< computing the mesh density employing a QuadTree
};

}

/**
 * \brief Reads and writes GMSH-files to and from OGS data structures.
 */
class GMSHInterface : public BaseLib::IO::Writer
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
	GMSHInterface (GeoLib::GEOObjects & geo_objs,
					bool include_stations_as_constraints,
					GMSH::MeshDensityAlgorithm mesh_density_algorithm,
					double param1, double param2, std::size_t param3,
					std::vector<std::string> & selected_geometries);

	~GMSHInterface();

	/**
	 * checks if there is a GMSH mesh file header
	 * @param fname the file name of the mesh (including the path)
	 * @return true, if the file seems to be a valid GMSH file, else false
	 */
	static bool isGMSHMeshFile (const std::string& fname);
	/**
	 * reads a mesh created by GMSH - this implementation is based on the former function GMSH2MSH
	 * @param fname the file name of the mesh (including the path)
	 * @return
	 */
	static MeshLib::Mesh* readGMSHMesh (std::string const& fname);

	/**
	 * Export script for writing geo files.
	 * To do this, all geometries currently loaded are merged, the merged result is written to a
	 * file and then the merged geometry is removed again.
	 * @return error code, i.e. 0 = okay, 1 = geo_objects is empty, 2 = error while merging, 3 = error writing file
	 */
	static int writeGeoFile(GeoLib::GEOObjects &geo_objects, std::string const& file_name);

protected:
	bool write();

private:
	/// Reads a mesh element from the input stream
	static std::pair<MeshLib::Element*, int> readElement(std::ifstream &in,
		std::vector<MeshLib::Node*> const& nodes,
		std::map<unsigned, unsigned> const& id_map);

	/**
	 * 1. get and merge data from _geo_objs
	 * 2. compute topological hierarchy
	 * @param out
	 */
	void writeGMSHInputFile(std::ostream & out);

	static void readNodeIDs(std::ifstream &in, unsigned n_nodes, std::vector<unsigned> &node_ids, std::map<unsigned, unsigned> const& id_map);

	void writePoints(std::ostream& out) const;

	std::size_t _n_lines;
	std::size_t _n_plane_sfc;

	GeoLib::GEOObjects & _geo_objs;
	std::vector<std::string>& _selected_geometries;
	std::string _gmsh_geo_name;
	std::list<GMSH::GMSHPolygonTree*> _polygon_tree_list;

	std::vector<FileIO::GMSH::GMSHPoint*> _gmsh_pnts;

	GMSH::GMSHMeshDensityStrategy *_mesh_density_strategy;
	/// Holds the inverse rotation matrix. The matrix is used in writePoints() to
	/// revert the rotation done in writeGMSHInputFile().
	MathLib::DenseMatrix<double> _inverse_rot_mat = MathLib::DenseMatrix<double>(3,3);
};
}

#endif /* GMSHINTERFACE_H_ */
