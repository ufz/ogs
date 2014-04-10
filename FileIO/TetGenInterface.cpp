/**
 * \file
 * \author Thomas Fischer
 * \date   2011-09-12
 * \brief  Implementation of the TetGenInterface class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <cstddef>
#include <string>
#include <fstream>

// BaseLib
#include "FileTools.h"
#include "StringTools.h"

// ThirdParty/logog
#include "logog/include/logog.hpp"

// FileIO
#include "TetGenInterface.h"

// MeshLib
#include "Mesh.h"
#include "Node.h"
#include "Elements/Element.h"
#include "Elements/Tet.h"

namespace FileIO
{
TetGenInterface::TetGenInterface() :
	_zero_based_idx (false)
{
}

TetGenInterface::~TetGenInterface()
{
}

bool TetGenInterface::readTetGenPoly (std::string const& poly_fname,
	                                  GeoLib::GEOObjects &geo_objects)
{
	std::ifstream poly_stream (poly_fname.c_str());

	if (!poly_stream)
	{
		ERR ("TetGenInterface::readTetGenPoly failed to open %s", poly_fname.c_str());
		return false;
	}

	std::vector<MeshLib::Node*> nodes;
	if (!readNodesFromStream (poly_stream, nodes)) 
	{
		// remove nodes read until now
		for (std::size_t k(0); k<nodes.size(); ++k)
			delete nodes[k];
		return false;
	}
	const std::size_t nNodes (nodes.size());
	std::vector<GeoLib::Point*> *points = new std::vector<GeoLib::Point*>;
	points->reserve(nNodes);
	for (std::size_t k(0); k<nNodes; ++k) 
	{
		points->push_back(new GeoLib::Point(nodes[k]->getCoords()));
		delete nodes[k];
	}
	std::vector<GeoLib::Surface*> *surfaces = new std::vector<GeoLib::Surface*>;
	if (!parseFacets(poly_stream, *surfaces, *points))
	{
		// remove surfaces read until now but keep the points
		for (std::size_t k=0; k<surfaces->size(); k++)
			delete (*surfaces)[k];
		delete surfaces;
		surfaces = nullptr;
	}

	std::string geo_name (BaseLib::extractBaseNameWithoutExtension(poly_fname));
	geo_objects.addPointVec(points, geo_name);
	if (surfaces)
		geo_objects.addSurfaceVec(surfaces, geo_name);
	return true;
}

std::size_t TetGenInterface::getNFacets(std::ifstream &input) const
{
	std::string line;
	while (!input.fail())
	{
		getline (input, line);
		if (input.fail())
		{
			ERR("TetGenInterface::getNFacets(): Error reading number of facets.");
			return false;
		}
		
		BaseLib::simplify(line);
		if (line.empty() || line.compare(0,1,"#") == 0)
			continue;

		const std::list<std::string> fields = BaseLib::splitString(line, ' ');
		return BaseLib::str2number<size_t> (*fields.begin());
		// here this line also includes a flag for boundary markers which we ignore for now
	}
	return false;
}

bool TetGenInterface::parseFacets(std::ifstream &input, 
                                  std::vector<GeoLib::Surface*> &surfaces,
                                  std::vector<GeoLib::Point*> &points)
{
	const std::size_t nFacets (this->getNFacets(input));
	std::size_t nMultPolys (0);
	std::string line;
	surfaces.reserve(nFacets);
	std::list<std::string>::const_iterator it;

	const unsigned offset = (_zero_based_idx) ? 0 : 1;
	for (std::size_t k(0); k<nFacets && !input.fail(); k++)
	{
		getline (input, line);
		if (input.fail())
		{
			ERR("TetGenInterface::parseFacets(): Error reading facet %d.", k);
			return false;
		}

		BaseLib::simplify(line);
		if (line.empty() || line.compare(0,1,"#") == 0)
		{
			k--;
			continue;
		}
		
		// read facets
		const std::list<std::string> poly_def_fields = BaseLib::splitString(line, ' ');
		it = poly_def_fields.begin();
		const std::size_t nPolys     = BaseLib::str2number<std::size_t>(*it);
		const std::size_t nPolyHoles = (poly_def_fields.size()>1) ? BaseLib::str2number<std::size_t>(*(++it)) : 0;
		// here this line also potentially includes a boundary marker which we ignore for now
		nMultPolys += (nPolys-1);

		// read polys
		for (std::size_t i(0); i<nPolys && !input.fail(); ++i)
		{
			getline (input, line);
			BaseLib::simplify(line);
			if (line.empty() || line.compare(0,1,"#") == 0)
			{
				i--;
				continue;
			}

			const std::list<std::string> point_fields = BaseLib::splitString(line, ' ');
			it = point_fields.begin();
			const std::size_t nPoints = BaseLib::str2number<std::size_t>(*it);
			if (point_fields.size() > nPoints)
			{
				GeoLib::Polyline polyline(points);
				for (std::size_t j(0); j<nPoints; ++j)
					polyline.addPoint(BaseLib::str2number<std::size_t>(*(++it))-offset);
				
				polyline.closePolyline();
				surfaces.push_back(GeoLib::Surface::createSurface(polyline));
			}
			else
			{
				ERR("TetGenInterface::parseFacets(): Error reading points for polygon %d of facet %d.", i, k);
				return false;
			}
		}
		for (std::size_t j(0); j<nPolyHoles && !input.fail(); ++j)
			getline(input, line);
			// Here are points defined which are located in holes within the surface. We ignore these as they are not part of the actual geometry.
	}
	// here the poly-file potentially defines a number of points to mark holes within the volumes defined by the facets, these are ignored for now
	// here the poly-file potentially defines a number of region attributes, these are ignored for now

	if (surfaces.size() == nFacets+nMultPolys)
		return true;

	ERR ("TetGenInterface::parseFacets(): Number of expected surfaces (%d) does not match number of found surfaces (%d).", nFacets+nMultPolys, surfaces.size());
	return false;
}

MeshLib::Mesh* TetGenInterface::readTetGenMesh (std::string const& nodes_fname,
                                                std::string const& ele_fname)
{
	std::ifstream ins_nodes (nodes_fname.c_str());
	std::ifstream ins_ele (ele_fname.c_str());

	if (!ins_nodes || !ins_ele)
	{
		if (!ins_nodes)
			ERR ("TetGenInterface::readTetGenMesh failed to open %s", nodes_fname.c_str());
		if (!ins_ele)
			ERR ("TetGenInterface::readTetGenMesh failed to open %s", ele_fname.c_str());
		return nullptr;
	}

	std::vector<MeshLib::Node*> nodes;
	if (!readNodesFromStream (ins_nodes, nodes)) {
		// remove nodes read until now
		for (std::size_t k(0); k<nodes.size(); k++) {
			delete nodes[k];
		}
		return nullptr;
	}

	std::vector<MeshLib::Element*> elements;
	if (!readElementsFromStream (ins_ele, elements, nodes)) {
		// remove elements read until now
		for (std::size_t k(0); k<elements.size(); k++) {
			delete elements[k];
		}
		// remove nodes
		for (std::size_t k(0); k<nodes.size(); k++) {
			delete nodes[k];
		}
		return nullptr;
	}

	const std::string mesh_name (BaseLib::extractBaseNameWithoutExtension(nodes_fname));
	return new MeshLib::Mesh(mesh_name, nodes, elements);
}

bool TetGenInterface::readNodesFromStream (std::ifstream &ins, 
                                           std::vector<MeshLib::Node*> &nodes)
{
	std::string line;
	getline (ins, line);
	size_t pos_beg (line.find_first_not_of(" "));
	size_t n_nodes, dim, n_attributes;
	bool boundary_markers;

	while (!ins.fail())
	{
		line = line.substr(pos_beg);
		if (line.compare(0,1,"#") == 0)
		{
			// this line is a comment - skip
			getline (ins, line);
			pos_beg = line.find_first_not_of(" ");
			continue;
		}
		// read header line
		bool header_okay = parseNodesFileHeader(line, n_nodes, dim, n_attributes, boundary_markers);
		if (!header_okay)
			return false;
		if (!parseNodes(ins, nodes, n_nodes, dim))
			return false;
		return true;
	}
	return false;	
}

bool TetGenInterface::parseNodesFileHeader(std::string &line, 
                                           std::size_t &n_nodes, 
                                           std::size_t &dim,
                                           std::size_t &n_attributes, 
                                           bool &boundary_markers) const
{
	std::size_t pos_beg, pos_end;

	// number of nodes
	pos_beg = line.find_first_not_of (" ");
	pos_end = line.find_first_of(" ", pos_beg);
	if (pos_beg != std::string::npos && pos_end != std::string::npos)
		n_nodes = BaseLib::str2number<size_t> (line.substr(pos_beg, pos_end - pos_beg));
	else
	{
		ERR("TetGenInterface::parseNodesFileHeader(): could not read number of nodes specified in header.");
		return false;
	}
	// dimension
	pos_beg = line.find_first_not_of (" ", pos_end);
	pos_end = line.find_first_of(" ", pos_beg);
	dim = BaseLib::str2number<size_t> (line.substr(pos_beg, pos_end - pos_beg));
	// number of attributes
	pos_beg = line.find_first_not_of (" ", pos_end);
	pos_end = line.find_first_of(" ", pos_beg);
	n_attributes = BaseLib::str2number<size_t> (line.substr(pos_beg, pos_end - pos_beg));
	// boundary marker at nodes?
	pos_beg = line.find_first_not_of (" ", pos_end);
	pos_end = line.find_first_of(" ", pos_beg);
	if (pos_end == std::string::npos)
		pos_end = line.size();
	if ((line.substr(pos_beg, pos_end - pos_beg)).compare("1") == 0)
		boundary_markers = true;
	else
		boundary_markers = false;

	return true;
}

bool TetGenInterface::parseNodes(std::ifstream &ins, 
                                 std::vector<MeshLib::Node*> &nodes, 
                                 std::size_t n_nodes, 
                                 std::size_t dim)
{
	std::size_t pos_beg, pos_end, id;
	std::string line;
	double* coordinates (static_cast<double*> (alloca (sizeof(double) * dim)));
	nodes.reserve(n_nodes);

	for (std::size_t k(0); k < n_nodes && !ins.fail(); k++) 
	{
		getline(ins, line);
		if (ins.fail()) 
		{
			ERR("TetGenInterface::parseNodes(): Error reading node %d.", k);
			return false;
		}

		pos_end = 0;
		pos_beg = line.find_first_not_of(" ", pos_end);
		pos_end = line.find_first_of(" \n", pos_beg);

		if (line.empty() || pos_beg==pos_end || line.compare(pos_beg,1,"#") == 0)
		{
			k--;
			continue;
		}

		if (pos_beg != std::string::npos && pos_end != std::string::npos) {
			id = BaseLib::str2number<size_t> (line.substr(pos_beg, pos_end - pos_beg));
			if (k == 0 && id == 0)
				_zero_based_idx = true;
		} else {
			ERR("TetGenInterface::parseNodes(): Error reading ID of node %d.", k);
			return false;
		}
		// read coordinates
		const unsigned offset = (_zero_based_idx) ? 0 : 1;
		for (std::size_t i(0); i < dim; i++) {
			pos_beg = line.find_first_not_of(" ", pos_end);
			pos_end = line.find_first_of(" \n", pos_beg);
			if (pos_end == std::string::npos) pos_end = line.size();
			if (pos_beg != std::string::npos)
				coordinates[i] = BaseLib::str2number<double>(line.substr(pos_beg, pos_end-pos_beg));
			else {
				ERR("TetGenInterface::parseNodes(): error reading coordinate %d of node %d.", i, k);
				return false;
			}
		}

		nodes.push_back(new MeshLib::Node(coordinates, id-offset));
		// read attributes and boundary markers ... - at the moment we do not use this information
	}

	return true;
}

bool TetGenInterface::readElementsFromStream(std::ifstream &ins, 
                                             std::vector<MeshLib::Element*> &elements, 
                                             const std::vector<MeshLib::Node*> &nodes)
{
	std::string line;
	getline (ins, line);
	std::size_t pos_beg (line.find_first_not_of(" "));
	std::size_t n_tets, n_nodes_per_tet;
	bool region_attributes;

	while (!ins.fail())
	{
		line = line.substr(pos_beg);
		if (line.compare(0,1,"#") == 0)
		{
			// this line is a comment - skip
			getline (ins, line);
			pos_beg = line.find_first_not_of(" ");
			continue;
		}
		
		// read header line
		bool header_okay = parseElementsFileHeader(line, n_tets, n_nodes_per_tet, region_attributes);
		if (!header_okay)
			return false;
		if (!parseElements(ins, elements, nodes, n_tets, n_nodes_per_tet, region_attributes))
			return false;
		return true;
	}
	return false;
}

bool TetGenInterface::parseElementsFileHeader(std::string &line,
                                              std::size_t& n_tets,
                                              std::size_t& n_nodes_per_tet,
                                              bool& region_attribute) const
{
	std::size_t pos_beg, pos_end;

	// number of tetrahedras
	pos_beg = line.find_first_not_of (" ");
	pos_end = line.find_first_of(" ", pos_beg);
	if (pos_beg != std::string::npos && pos_end != std::string::npos)
		n_tets = BaseLib::str2number<size_t> (line.substr(pos_beg, pos_end - pos_beg));
	else {
		ERR("TetGenInterface::parseElementsFileHeader(): Could not read number of tetrahedra specified in header.");
		return false;
	}
	// nodes per tet - either 4 or 10
	pos_beg = line.find_first_not_of (" \t", pos_end);
	pos_end = line.find_first_of(" \t", pos_beg);
	n_nodes_per_tet = BaseLib::str2number<size_t> (line.substr(pos_beg, pos_end - pos_beg));
	// region attribute at tetrahedra?
	pos_beg = line.find_first_not_of (" \t", pos_end);
	pos_end = line.find_first_of(" \t\n", pos_beg);
	if (pos_end == std::string::npos)
		pos_end = line.size();
	if ((line.substr(pos_beg, pos_end - pos_beg)).compare("1") == 0)
		region_attribute = true;
	else
		region_attribute = false;

	return true;
}

bool TetGenInterface::parseElements(std::ifstream& ins, 
                                    std::vector<MeshLib::Element*> &elements, 
                                    const std::vector<MeshLib::Node*> &nodes, 
                                    std::size_t n_tets, 
                                    std::size_t n_nodes_per_tet,
                                    bool region_attribute)
{
	std::size_t pos_beg, pos_end, id;
	std::string line;
	std::size_t* ids (static_cast<size_t*>(alloca (sizeof (size_t) * n_nodes_per_tet)));
	elements.reserve(n_tets);

	const unsigned offset = (_zero_based_idx) ? 0 : 1;
	for (std::size_t k(0); k < n_tets && !ins.fail(); k++)
	{
		getline (ins, line);
		if (ins.fail())
		{
			ERR("TetGenInterface::parseElements(): Error reading tetrahedron %d.", k);
			return false;
		}

		pos_end = 0;
		pos_beg = line.find_first_not_of(" ", pos_end);
		pos_end = line.find_first_of(" \n", pos_beg);

		if (line.empty() || pos_beg==pos_end || line.compare(pos_beg,1,"#") == 0)
		{
			k--;
			continue;
		}
				
		if (pos_beg != std::string::npos && pos_end != std::string::npos)
			id = BaseLib::str2number<size_t>(line.substr(pos_beg, pos_end - pos_beg));
		else {
			ERR("TetGenInterface::parseElements(): Error reading id of tetrahedron %d.", k);
			return false;
		}
		// read node ids
		for (std::size_t i(0); i < n_nodes_per_tet; i++)
		{
			pos_beg = line.find_first_not_of(" ", pos_end);
			pos_end = line.find_first_of(" ", pos_beg);
			if (pos_end == std::string::npos)
				pos_end = line.size();
			if (pos_beg != std::string::npos && pos_end != std::string::npos)
				ids[i] = BaseLib::str2number<std::size_t>(line.substr(pos_beg, pos_end - pos_beg)) - offset;
			else
			{
				ERR("TetGenInterface::parseElements(): Error reading node %d of tetrahedron %d.", i, k);
				return false;
			}
		}

		// read region attribute - this is something like material group
		unsigned region (0);
		if (region_attribute) {
			pos_beg = line.find_first_not_of(" ", pos_end);
			pos_end = line.find_first_of(" ", pos_beg);
			if (pos_end == std::string::npos) pos_end = line.size();
			if (pos_beg != std::string::npos && pos_end != std::string::npos)
				region = BaseLib::str2number<unsigned> (line.substr(pos_beg, pos_end - pos_beg));
			else {
				ERR("TetGenInterface::parseElements(): Error reading region attribute of tetrahedron %d.", k);
				return false;
			}
		}
		// insert new element into vector
		MeshLib::Node** tet_nodes = new MeshLib::Node*[4];
		for (unsigned k(0); k<4; k++) {
			tet_nodes[k] = nodes[ids[k]];
		}
		elements.push_back (new MeshLib::Tet(tet_nodes, region));
	}
	return true;
}

bool TetGenInterface::writeTetGenPoly(const std::string &file_name, 
                                      const GeoLib::GEOObjects &geo_objects, 
                                      const std::string &geo_name) const
{
	std::vector<GeoLib::Point*> const*const points = geo_objects.getPointVec(geo_name);
	std::vector<GeoLib::Surface*> const*const surfaces = geo_objects.getSurfaceVec(geo_name);

	if (points==nullptr)
	{
		ERR ("Geometry %s not found.", geo_name.c_str());
		return false;
	}
	if (surfaces==nullptr)
		WARN ("No surfaces found for geometry %s. Writing points only.", geo_name.c_str());

	std::ofstream out( file_name.c_str(), std::ios::out );
	// the points header
	const std::size_t nPoints (points->size());
	out << nPoints << " 3 0 0\n";
	// the point list
	for (std::size_t i=0; i<nPoints; ++i)
		out << i << "  " << (*(*points)[i])[0] << " " << (*(*points)[i])[1] << " " << (*(*points)[i])[2] << "\n";
	// the surfaces header
	const std::size_t nSurfaces = (surfaces) ? surfaces->size() : 0;
	out << nSurfaces << " 0\n";
	// the facets list
	for (std::size_t i=0; i<nSurfaces; ++i)
	{
		// the number of polys per facet
		const std::size_t nTriangles ((*surfaces)[i]->getNTriangles());
		out << nTriangles << "\n";
		// the poly list
		for (std::size_t j=0; j<nTriangles; ++j)
		{
			const GeoLib::Triangle &tri = *(*(*surfaces)[i])[j];
			out << "3  " << tri[0] << " " << tri[1] << " " << tri[2] << "\n";			
		}
	}
	out << "0\n"; // the polygon holes list
	out << "0\n"; // the region attribues list
	INFO ("TetGenInterface::writeTetGenPoly() - %d points and %d surfaces successfully written.", nPoints, nSurfaces);
	out.close();
	return true;
}

} // end namespace FileIO
