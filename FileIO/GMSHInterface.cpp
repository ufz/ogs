/**
 * \file
 * \author Thomas Fischer
 * \date   2010-04-29
 * \brief  Implementation of the GMSHInterface class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * @file GMSHInterface.cpp
 * @date 2010-04-29
 * @author Thomas Fischer
 */

#include <fstream>
#include <vector>

// ThirdParty/logog
#include "logog/include/logog.hpp"

// BaseLib
#include "BaseLib/BuildInfo.h"
#include "FileTools.h"
#include "StringTools.h"

// FileIO
#include "GMSHInterface.h"
#include "GmshIO/GMSHAdaptiveMeshDensity.h"
#include "GmshIO/GMSHFixedMeshDensity.h"
#include "GmshIO/GMSHNoMeshDensity.h"

// GeoLib
#include "Point.h"
#include "Polygon.h"
#include "Polyline.h"
#include "PolylineWithSegmentMarker.h"
#include "QuadTree.h"

// MSH
#include "Elements/Line.h"
#include "Elements/Hex.h"
#include "Elements/Prism.h"
#include "Elements/Pyramid.h"
#include "Elements/Quad.h"
#include "Elements/Tet.h"
#include "Elements/Tri.h"
#include "Mesh.h"
#include "Node.h"

namespace FileIO
{
GMSHInterface::GMSHInterface(GeoLib::GEOObjects & geo_objs,
                             bool /*include_stations_as_constraints*/,
                             GMSH::MeshDensityAlgorithm mesh_density_algorithm,
                             double param1,
                             double param2,
                             size_t param3,
                             std::vector<std::string>& selected_geometries) :
	_n_lines(0), _n_plane_sfc(0), _geo_objs(geo_objs), _selected_geometries(selected_geometries)
{
	switch (mesh_density_algorithm) {
	case GMSH::MeshDensityAlgorithm::NoMeshDensity:
		_mesh_density_strategy = new GMSH::GMSHNoMeshDensity;
		break;
	case GMSH::MeshDensityAlgorithm::FixedMeshDensity:
		_mesh_density_strategy = new GMSH::GMSHFixedMeshDensity(param1);
		break;
	case GMSH::MeshDensityAlgorithm::AdaptiveMeshDensity:
		_mesh_density_strategy = new GMSH::GMSHAdaptiveMeshDensity(param1, param2, param3);
		break;
	}
}

bool GMSHInterface::isGMSHMeshFile(const std::string& fname)
{
	std::ifstream input(fname.c_str());

	if (!input) {
		ERR("GMSHInterface::isGMSHMeshFile(): Could not open file %s.", fname.c_str());
		return false;
	}

	std::string header_first_line;
	input >> header_first_line;
	if (header_first_line.find("$MeshFormat") != std::string::npos) {
		// read version
		std::string version;
		getline(input, version);
		getline(input, version);
		INFO("GMSHInterface::isGMSHMeshFile(): Found GMSH mesh file version: %s.",
		     version.c_str());
		input.close();
		return true;
	}

	return false;
}

MeshLib::Mesh* GMSHInterface::readGMSHMesh(std::string const& fname)
{
	std::string line;
	std::ifstream in(fname.c_str(), std::ios::in);
	getline(in, line); // Node keyword
	std::vector<MeshLib::Node*> nodes;
	std::vector<MeshLib::Element*> elements;

	if (line.find("$MeshFormat") != std::string::npos)
	{
		getline(in, line); // version-number file-type data-size
		if (line.substr(0,3).compare("2.2") != 0) {
			WARN("Wrong gmsh file format version.");
			return nullptr;
		}
		if (line.substr(4,1).compare("0") != 0) {
			WARN("Currently reading gmsh binary file type is not supported.");
			return nullptr;
		}
		getline(in, line); //$EndMeshFormat
		getline(in, line); //$Nodes Keywords

		size_t n_nodes(0);
		size_t n_elements(0);
		while (line.find("$EndElements") == std::string::npos)
		{
			// Node data
			long id;
			double x, y, z;
			in >> n_nodes >> std::ws;
			nodes.resize(n_nodes);
			std::map<unsigned, unsigned> id_map;
			for (size_t i = 0; i < n_nodes; i++) {
				in >> id >> x >> y >> z >> std::ws;
				id_map.insert(std::map<unsigned, unsigned>::value_type(id, i));
				nodes[i] = new MeshLib::Node(x,y,z,id);
			}
			getline(in, line); // End Node keyword $EndNodes

			// Element data
			getline(in, line); // Element keyword $Elements
			in >> n_elements >> std::ws; // number-of-elements
			elements.reserve(n_elements);
			unsigned idx, type, n_tags, dummy, mat_id;
			for (size_t i = 0; i < n_elements; i++)
			{
				MeshLib::Element* elem (NULL);
				std::vector<unsigned> node_ids;
				std::vector<MeshLib::Node*> elem_nodes;
				in >> idx >> type >> n_tags >> dummy >> mat_id;

				// skip tags
				for (size_t j = 2; j < n_tags; j++)
					in >> dummy;

				switch (type)
				{
				case 1: {
					readNodeIDs(in, 2, node_ids, id_map);
					// edge_nodes array will be deleted from Line object
					MeshLib::Node** edge_nodes = new MeshLib::Node*[2];
					edge_nodes[0] = nodes[node_ids[0]];
					edge_nodes[1] = nodes[node_ids[1]];
					elem = new MeshLib::Line(edge_nodes, 0);
					break;
				}
				case 2: {
					readNodeIDs(in, 3, node_ids, id_map);
					MeshLib::Node** tri_nodes = new MeshLib::Node*[3];
					tri_nodes[0] = nodes[node_ids[2]];
					tri_nodes[1] = nodes[node_ids[1]];
					tri_nodes[2] = nodes[node_ids[0]];
					elem = new MeshLib::Tri(tri_nodes, mat_id);
					break;
				}
				case 3: {
					readNodeIDs(in, 4, node_ids, id_map);
					MeshLib::Node** quad_nodes = new MeshLib::Node*[4];
					for (unsigned k(0); k < 4; k++)
						quad_nodes[k] = nodes[node_ids[k]];
					elem = new MeshLib::Quad(quad_nodes, mat_id);
					break;
				}
				case 4: {
					readNodeIDs(in, 4, node_ids, id_map);
					MeshLib::Node** tet_nodes = new MeshLib::Node*[5];
					for (unsigned k(0); k < 4; k++)
						tet_nodes[k] = nodes[node_ids[k]];
					elem = new MeshLib::Tet(tet_nodes, mat_id);
					break;
				}
				case 5: {
					readNodeIDs(in, 8, node_ids, id_map);
					MeshLib::Node** hex_nodes = new MeshLib::Node*[8];
					for (unsigned k(0); k < 8; k++)
						hex_nodes[k] = nodes[node_ids[k]];
					elem = new MeshLib::Hex(hex_nodes, mat_id);
					break;
				}
				case 6: {
					readNodeIDs(in, 6, node_ids, id_map);
					MeshLib::Node** prism_nodes = new MeshLib::Node*[6];
					for (unsigned k(0); k < 6; k++)
						prism_nodes[k] = nodes[node_ids[k]];
					elem = new MeshLib::Prism(prism_nodes, mat_id);
					break;
				}
				case 7: {
					readNodeIDs(in, 5, node_ids, id_map);
					MeshLib::Node** pyramid_nodes = new MeshLib::Node*[5];
					for (unsigned k(0); k < 5; k++)
						pyramid_nodes[k] = nodes[node_ids[k]];
					elem = new MeshLib::Pyramid(pyramid_nodes, mat_id);
					break;
				}
				case 15:
					in >> dummy; // skip rest of line
					continue;
					break;
				default:
						WARN("GMSHInterface::readGMSHMesh(): Unknown element type %d.", type);
				}
				in >> std::ws;

				if (type > 0 && type < 8)
					elements.push_back(elem);
			}

			getline(in, line); // END keyword
		}
	}
	in.close();
	if (elements.empty()) {
		for (auto it(nodes.begin()); it != nodes.end(); ++it) {
			delete *it;
		}
		return nullptr;
	}
	return new MeshLib::Mesh(BaseLib::extractBaseNameWithoutExtension(fname), nodes, elements);
}

void GMSHInterface::readNodeIDs(std::ifstream &in,
                                unsigned n_nodes,
                                std::vector<unsigned> &node_ids,
                                std::map<unsigned, unsigned> &id_map)
{
	unsigned idx;
	for (unsigned i = 0; i < n_nodes; i++)
	{
		in >> idx;
		node_ids.push_back(id_map[idx]);
	}
}

bool GMSHInterface::write()
{
	_out << "// GMSH input file created by OpenGeoSys " << BaseLib::BuildInfo::ogs_version_and_persons;
#ifdef BUILD_TIMESTAMP
	_out << " built on " << BaseLib::BuildInfo::build_timestamp;
#endif
	_out << "\n\n";

	writeGMSHInputFile(_out);
	return true;
}

void GMSHInterface::writeGMSHInputFile(std::ostream& out)
{
	DBUG("GMSHInterface::writeGMSHInputFile(): get data from GEOObjects.");

	if (_selected_geometries.empty())
		return;

	bool remove_geometry(false);
	// *** get and merge data from _geo_objs
	if (_selected_geometries.size() > 1) {
		_gmsh_geo_name = "GMSHGeometry";
		remove_geometry = true;
		_geo_objs.mergeGeometries(_selected_geometries, _gmsh_geo_name);
	} else {
		_gmsh_geo_name = _selected_geometries[0];
		remove_geometry = false;
	}
	std::vector<GeoLib::Point*> * merged_pnts(const_cast<std::vector<GeoLib::Point*> *>(_geo_objs.getPointVec(_gmsh_geo_name)));
	if (! merged_pnts) {
		ERR("GMSHInterface::writeGMSHInputFile(): Did not found any points.");
		return;
	} else {
		const size_t n_pnts(merged_pnts->size());
		for (size_t k(0); k<n_pnts; k++) {
			(*((*merged_pnts)[k]))[2] = 0.0;
		}
	}
	std::vector<GeoLib::Polyline*> const* merged_plys(_geo_objs.getPolylineVec(_gmsh_geo_name));
	DBUG("GMSHInterface::writeGMSHInputFile(): \t ok.");

	// *** compute topological hierarchy of polygons
	if (merged_plys) {
		for (std::vector<GeoLib::Polyline*>::const_iterator it(merged_plys->begin());
			it!=merged_plys->end(); ++it) {
			if ((*it)->isClosed()) {
				_polygon_tree_list.push_back(new GMSH::GMSHPolygonTree(new GeoLib::Polygon(*(*it), true), NULL, _geo_objs, _gmsh_geo_name, _mesh_density_strategy));
			}
		}
		DBUG("GMSHInterface::writeGMSHInputFile(): Compute topological hierarchy - detected %d polygons.", _polygon_tree_list.size());
		GeoLib::createPolygonTrees<GMSH::GMSHPolygonTree>(_polygon_tree_list);
		DBUG("GMSHInterface::writeGMSHInputFile(): Compute topological hierarchy - calculated %d polygon trees.", _polygon_tree_list.size());
	} else {
		return;
	}

	// *** insert stations and polylines (except polygons) in the appropriate object of
	//     class GMSHPolygonTree
	// *** insert stations
	const size_t n_geo_names(_selected_geometries.size());
	for (size_t j(0); j < n_geo_names; j++) {
		const std::vector<GeoLib::Point*>* stations (_geo_objs.getStationVec(_selected_geometries[j]));
		if (stations) {
			const size_t n_stations(stations->size());
			for (size_t k(0); k < n_stations; k++) {
				bool found(false);
				for (std::list<GMSH::GMSHPolygonTree*>::iterator it(_polygon_tree_list.begin());
					it != _polygon_tree_list.end() && !found; ++it) {
					if ((*it)->insertStation((*stations)[k])) {
						found = true;
					}
				}
			}
		}
	}
	// *** insert polylines
	const size_t n_plys(merged_plys->size());
	for (size_t k(0); k<n_plys; k++) {
		if (! (*merged_plys)[k]->isClosed()) {
			for (std::list<GMSH::GMSHPolygonTree*>::iterator it(_polygon_tree_list.begin());
				it != _polygon_tree_list.end(); ++it) {
				(*it)->insertPolyline(new GeoLib::PolylineWithSegmentMarker(*(*merged_plys)[k]));
			}
		}
	}

	// *** init mesh density strategies
	for (std::list<GMSH::GMSHPolygonTree*>::iterator it(_polygon_tree_list.begin());
		it != _polygon_tree_list.end(); ++it) {
		(*it)->initMeshDensityStrategy();
	}

	// *** create GMSH data structures
	const size_t n_merged_pnts(merged_pnts->size());
	_gmsh_pnts.resize(n_merged_pnts);
	for (size_t k(0); k<n_merged_pnts; k++) {
		_gmsh_pnts[k] = NULL;
	}
	for (std::list<GMSH::GMSHPolygonTree*>::iterator it(_polygon_tree_list.begin());
		it != _polygon_tree_list.end(); ++it) {
		(*it)->createGMSHPoints(_gmsh_pnts);
	}

	// *** finally write data :-)
	writePoints(out);
	size_t pnt_id_offset(_gmsh_pnts.size());
	for (std::list<GMSH::GMSHPolygonTree*>::iterator it(_polygon_tree_list.begin());
		it != _polygon_tree_list.end(); ++it) {
		(*it)->writeLineLoop(_n_lines, _n_plane_sfc, out);
		(*it)->writeSubPolygonsAsLineConstraints(_n_lines, _n_plane_sfc-1, out);
		(*it)->writeLineConstraints(_n_lines, _n_plane_sfc-1, out);
		(*it)->writeStations(pnt_id_offset, _n_plane_sfc-1, out);
		(*it)->writeAdditionalPointData(pnt_id_offset, _n_plane_sfc-1, out);
	}

	if (remove_geometry) {
		_geo_objs.removeSurfaceVec(_gmsh_geo_name);
		_geo_objs.removePolylineVec(_gmsh_geo_name);
		_geo_objs.removePointVec(_gmsh_geo_name);
	}
}

void GMSHInterface::writePoints(std::ostream& out) const
{
	const size_t n_gmsh_pnts(_gmsh_pnts.size());
	for (size_t k(0); k<n_gmsh_pnts; k++) {
		if (_gmsh_pnts[k]) {
			out << *(_gmsh_pnts[k]) << "\n";
		}
	}
}

} // end namespace FileIO
