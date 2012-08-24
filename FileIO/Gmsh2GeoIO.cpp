/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 * \file Gmsh2GeoIO.cpp
 *
 *  Created on 2011-08-18 by Thomas Fischer
 */

#include <fstream>
#include <vector>

// Base
#include "StringTools.h"

#include "GEOObjects.h"
#include "Gmsh2GeoIO.h"

namespace FileIO
{
void Gmsh2GeoIO::loadMeshAsGeometry (std::string & fname, GeoLib::GEOObjects* geo)
{
	// open file
	std::ifstream ins (fname.c_str());
	if (!ins)
	{
		std::cout << "could not open file " << fname << std::endl;
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
	const size_t n_pnts (str2number<size_t>(line));
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
		double x (str2number<double>(line.substr(pos_beg, pos_end - pos_beg)));
		// parse y coordinate
		pos_beg = pos_end + 1;
		pos_end = line.find(" ", pos_beg);
		double y (str2number<double>(line.substr(pos_beg, pos_end - pos_beg)));
		// parse z coordinate
		pos_beg = pos_end + 1;
		pos_end = line.find("\n", pos_beg);
		double z (str2number<double>(line.substr(pos_beg, pos_end - pos_beg)));

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
	const size_t n_elements (str2number<size_t>(line));
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
		size_t ele_type (str2number<size_t>(line.substr(pos_beg, pos_end - pos_beg)));
		if (ele_type == 2) // read 3 node triangle
		{ // parse number of tags
			pos_beg = pos_end + 1;
			pos_end = line.find(" ", pos_beg);
			const size_t n_tags (str2number<size_t>(line.substr(pos_beg,
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
			const size_t id0 (str2number<size_t>(line.substr(pos_beg,
			                                                 pos_end - pos_beg)) - 1); // shift -1!
			// parse second id of triangle
			pos_beg = pos_end + 1;
			pos_end = line.find(" ", pos_beg);
			const size_t id1 (str2number<size_t>(line.substr(pos_beg,
			                                                 pos_end - pos_beg)) - 1); // shift -1!
			// parse third id of triangle
			pos_beg = pos_end + 1;
			pos_end = line.find(" ", pos_beg);
			const size_t id2 (str2number<size_t>(line.substr(pos_beg,
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
} // end namespace FileIO
