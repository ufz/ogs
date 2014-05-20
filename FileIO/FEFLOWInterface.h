/**
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#ifndef FEFLOWINTERFACE_H_
#define FEFLOWINTERFACE_H_

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "GEOObjects.h"
#include "MeshEnums.h"

class QDomElement;
class QString;

namespace MeshLib
{
class Mesh;
class Element;
class Node;
}

namespace FileIO
{

/**
 * Read FEFLOW files into OGS data structure
 */
class FEFLOWInterface
{
public:
	/// Constructor
	explicit FEFLOWInterface(GeoLib::GEOObjects* obj)
			: _geoObjects(obj)
	{
	}

	/**
	 * read a FEFLOW Model file (*.fem) in ASCII format (Version 5.4)
	 *
	 * This function reads mesh data in addition to geometry data given in Supermesh.
	 *
	 * @param filename  FEFLOW file name
	 * @return a pointer to a created OGS mesh
	 */
	MeshLib::Mesh* readFEFLOWFile(const std::string &filename);

private:
	// CLASS
	struct FEM_CLASS
	{
		unsigned problem_class;
		unsigned time_mode;
		unsigned orientation;
		unsigned dimension;
		unsigned n_layers3d;
		unsigned saturation_flag;
		unsigned save_fsize_rreal;
		unsigned save_fsize_creal;

		FEM_CLASS()
		{
			problem_class = 0;
			time_mode = 0;
			orientation = 0;
			dimension = 0;
			n_layers3d = 0;
			saturation_flag = 0;
			save_fsize_rreal = 0;
			save_fsize_creal = 0;
		}
	};

	// DIMENSION
	struct FEM_DIM
	{
		size_t n_nodes;
		size_t n_elements;
		unsigned n_nodes_of_element;
		unsigned n_steps;
		unsigned icrank;
		unsigned upwind;
		size_t obs;
		unsigned optim;
		unsigned aquifer_type;
		unsigned nwca;
		size_t np_cor;
		unsigned adaptive_mesh;
		unsigned sp_fem_pcs_id;
		unsigned sorption_type;
		unsigned reaction_type;
		unsigned dispersion_type;

		FEM_DIM()
		{
			n_nodes = 0;
			n_elements = 0;
			n_nodes_of_element = 0;
			n_steps = 0;
			icrank = 0;
			upwind = 0;
			obs = 0;
			optim = 0;
			aquifer_type = 0;
			nwca = 0;
			np_cor = 0;
			adaptive_mesh = 0;
			sp_fem_pcs_id = 0;
			sorption_type = 0;
			reaction_type = 0;
			dispersion_type = 0;
		}
	};

	/// read node index and create a mesh element
	MeshLib::Element* readElement(const FEM_DIM &fem_dim, const MeshElemType elem_type, const std::string& line, const std::vector<MeshLib::Node*> &nodes);

	/// read node coordinates
	void readNodeCoordinates(std::ifstream &in, const FEM_CLASS &fem_class, const FEM_DIM &fem_dim, std::vector<MeshLib::Node*> &nodes);

	/// read elevation data
	void readELEV(std::ifstream &in, const FEM_CLASS &fem_class, const FEM_DIM &fem_dim, std::vector<MeshLib::Node*> &vec_nodes);

	/// read Supermesh data
	void readSuperMesh(std::ifstream &feflow_file, const FEM_CLASS &fem_class, std::vector<GeoLib::Point*>** points, std::vector<GeoLib::Polyline*>** lines);

	//// read point data in Supermesh
	void readPoints(QDomElement &nodesEle, const std::string &tag, int dim, std::vector<GeoLib::Point*> &points);

	/// set element material IDs
	void setMaterialID(std::vector<MeshLib::Element*> &vec_elements, std::vector<GeoLib::Polyline*>* lines);

	//// Geometric objects
	GeoLib::GEOObjects* _geoObjects;
};
} // end namespace FileIO

#endif /* FEFLOWINTERFACE_H_ */
