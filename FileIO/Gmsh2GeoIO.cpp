/**
 * \file
 * \author Thomas Fischer
 * \date   2011-08-18
 * \brief  Implementation of the Gmsh2GeoIO class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file Gmsh2GeoIO.cpp
 *
 *  Created on 2011-08-18 by Thomas Fischer
 */

#include <fstream>
#include <vector>

// ThirdParty/logog
#include "logog/include/logog.hpp"

// BaseLib
#include "StringTools.h"

#include "GEOObjects.h"
#include "Gmsh2GeoIO.h"

#include "Mesh.h"
#include "Elements/Tri.h"
#include "Elements/Quad.h"
#include "Node.h"


namespace FileIO
{
void Gmsh2GeoIO::loadMeshAsGeometry (std::string & fname, GeoLib::GEOObjects* geo)
{
	// open file
	std::ifstream ins (fname.c_str());
	if (!ins)
	{
		WARN("Gmsh2GeoIO::loadMeshAsGeometry(): could not open file %s", fname.c_str());
		return;
	}

	std::string line;
	// read gmsh header
	getline (ins, line); // $MeshFormat
	getline (ins, line);
	getline (ins, line); // $EndMeshFormat

	// read nodes tag
	getline (ins, line);
	// read number of nodes
	getline (ins, line);
	const size_t n_pnts (BaseLib::str2number<size_t>(line));
	std::vector<GeoLib::Point*>* pnts (new std::vector<GeoLib::Point*>);
	for (size_t k(0); k < n_pnts; k++)
	{
		getline (ins, line);
		// parse id
		size_t pos_beg(0);
		size_t pos_end (line.find(" "));
		// the sub string line.substr(pos_beg, pos_end-pos_beg) represents the id
		// parse x coordinate
		pos_beg = pos_end + 1;
		pos_end = line.find(" ", pos_beg);
		double x (BaseLib::str2number<double>(line.substr(pos_beg, pos_end - pos_beg)));
		// parse y coordinate
		pos_beg = pos_end + 1;
		pos_end = line.find(" ", pos_beg);
		double y (BaseLib::str2number<double>(line.substr(pos_beg, pos_end - pos_beg)));
		// parse z coordinate
		pos_beg = pos_end + 1;
		pos_end = line.find("\n", pos_beg);
		double z (BaseLib::str2number<double>(line.substr(pos_beg, pos_end - pos_beg)));

		pnts->push_back (new GeoLib::Point (x,y,z));
	}
	// read end nodes tag
	getline (ins, line);

	geo->addPointVec (pnts, fname);

	std::vector<size_t> const& pnt_id_map (geo->getPointVecObj(fname)->getIDMap());
	// read element tag
	getline (ins, line);
	// read number of elements
	getline (ins, line);
	const size_t n_elements (BaseLib::str2number<size_t>(line));
	GeoLib::Surface* sfc (new GeoLib::Surface (*pnts));
	for (size_t k(0); k < n_elements; k++)
	{
		getline (ins, line);
		// parse id
		size_t pos_beg(0);
		size_t pos_end (line.find(" "));
		// the sub string line.substr(pos_beg, pos_end-pos_beg) represents the id
		// parse element type
		pos_beg = pos_end + 1;
		pos_end = line.find(" ", pos_beg);
		size_t ele_type (BaseLib::str2number<size_t>(line.substr(pos_beg, pos_end - pos_beg)));
		if (ele_type == 2) // read 3 node triangle
		{ // parse number of tags
			pos_beg = pos_end + 1;
			pos_end = line.find(" ", pos_beg);
			const size_t n_tags (BaseLib::str2number<size_t>(line.substr(pos_beg,
			                                                    pos_end - pos_beg)));
			// (over) read tags
			for (size_t j(0); j < n_tags; j++)
			{
				pos_beg = pos_end + 1;
				pos_end = line.find(" ", pos_beg);
			}
			// parse first id of triangle
			pos_beg = pos_end + 1;
			pos_end = line.find(" ", pos_beg);
			const size_t id0 (BaseLib::str2number<size_t>(line.substr(pos_beg,
			                                                 pos_end - pos_beg)) - 1); // shift -1!
			// parse second id of triangle
			pos_beg = pos_end + 1;
			pos_end = line.find(" ", pos_beg);
			const size_t id1 (BaseLib::str2number<size_t>(line.substr(pos_beg,
			                                                 pos_end - pos_beg)) - 1); // shift -1!
			// parse third id of triangle
			pos_beg = pos_end + 1;
			pos_end = line.find(" ", pos_beg);
			const size_t id2 (BaseLib::str2number<size_t>(line.substr(pos_beg,
			                                                 pos_end - pos_beg)) - 1); // shift -1!
			sfc->addTriangle (pnt_id_map[id0], pnt_id_map[id1], pnt_id_map[id2]);
		}
	}
	// read end element tag
	getline (ins, line);

	std::vector<GeoLib::Surface*>* sfcs (new std::vector<GeoLib::Surface*>);
	sfcs->push_back(sfc);
	geo->addSurfaceVec (sfcs, fname);
}


bool Gmsh2GeoIO::convertMeshToGeo(const MeshLib::Mesh &mesh, GeoLib::GEOObjects* geo_objects)
{
	if (mesh.getDimension() != 2)
	{
		ERR ("Mesh to geometry conversion is only working for 2D meshes.");
		return false;
	}
	
	// nodes to points conversion
	const unsigned nNodes (mesh.getNNodes());
	std::vector<GeoLib::Point*> *points = new std::vector<GeoLib::Point*>(nNodes);
	const std::vector<MeshLib::Node*> nodes = mesh.getNodes();

	for (unsigned i=0; i<nNodes; ++i)
		(*points)[i] = new GeoLib::Point(static_cast<GeoLib::Point>(*nodes[i]));


	// elements to surface triangles conversion
	const std::vector<MeshLib::Element*> elements = mesh.getElements();
	GeoLib::Surface* sfc = new GeoLib::Surface(*points);
	const unsigned nElems (mesh.getNElements());

	for (unsigned i=0; i<nElems; ++i)
	{
		MeshLib::Element* e (elements[i]);
		if (e->getGeomType() == MshElemType::TRIANGLE)
			sfc->addTriangle(e->getNodeIndex(0), e->getNodeIndex(1), e->getNodeIndex(2));
		if (e->getGeomType() == MshElemType::QUAD)
		{
			sfc->addTriangle(e->getNodeIndex(0), e->getNodeIndex(1), e->getNodeIndex(2));
			sfc->addTriangle(e->getNodeIndex(0), e->getNodeIndex(2), e->getNodeIndex(3));
		}
		// all other element types are ignored (i.e. lines)
	}
	
	std::vector<GeoLib::Surface*> *sfcs = new std::vector<GeoLib::Surface*>(1);
	(*sfcs)[0] = sfc;

	std::string mesh_name (mesh.getName());
	geo_objects->addPointVec(points, mesh_name);
	geo_objects->addSurfaceVec(sfcs, mesh_name);
	return true;
}


} // end namespace FileIO
