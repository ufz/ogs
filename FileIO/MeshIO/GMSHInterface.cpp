/*
 * GMSHInterface.cpp
 *
 *  Created on: Apr 29, 2010
 *      Author: TF
 */

#include <fstream>
#include <vector>

// Base
#include "swap.h"
#include "Configure.h"
#include "BuildInfo.h"

// FileIO
#include "GMSHInterface.h"
#include "GMSHNoMeshDensity.h"
#include "GMSHFixedMeshDensity.h"
#include "GMSHAdaptiveMeshDensity.h"

// GEOLIB
#include "Point.h"
#include "Polygon.h"
#include "Polyline.h"
#include "QuadTree.h"
#include "PolylineWithSegmentMarker.h"

// MSH
#include "msh_elem.h"
#include "msh_mesh.h"

namespace FileIO {
GMSHInterface::GMSHInterface(GEOLIB::GEOObjects & geo_objs, bool include_stations_as_constraints,
				GMSH::MeshDensityAlgorithm mesh_density_algorithm, double param1, double param2,
				size_t param3, std::vector<std::string>& selected_geometries) :
	_n_lines(0), _n_plane_sfc(0), _geo_objs(geo_objs), _selected_geometries(selected_geometries),
	_include_stations_as_constraints(include_stations_as_constraints)
{
	switch (mesh_density_algorithm) {
	case GMSH::NoMeshDensity:
		_mesh_density_strategy = new GMSHNoMeshDensity;
		break;
	case GMSH::FixedMeshDensity:
		_mesh_density_strategy = new GMSHFixedMeshDensity(param1);
		break;
	case GMSH::AdaptiveMeshDensity:
		_mesh_density_strategy = new GMSHAdaptiveMeshDensity(param1, param2, param3);
		break;
	}
}

bool GMSHInterface::isGMSHMeshFile(const std::string& fname)
{
	std::ifstream input(fname.c_str());

	if (!input) {
		std::cerr << "GMSHInterface::isGMSHMeshFile could not open file " << fname << std::endl;
		return false;
	}

	std::string header_first_line;
	input >> header_first_line;
	if (header_first_line.find("$MeshFormat") != std::string::npos) {
		// read version
		std::string version;
		getline(input, version);
		getline(input, version);
		std::cerr << "found GMSH mesh file version: " << version << std::endl;
		input.close();
		return true;
	}

	return false;
}

void GMSHInterface::readGMSHMesh(std::string const& fname, MeshLib::CFEMesh* mesh)
{
	std::string line;
	std::ifstream in(fname.c_str(), std::ios::in);
	getline(in, line); // Node keyword

	if (line.find("$MeshFormat") != std::string::npos) {
		getline(in, line); // version-number file-type data-size
		getline(in, line); //$EndMeshFormat
		getline(in, line); //$Nodes Keywords

		size_t n_nodes(0);
		size_t n_elements(0);
		while (line.find("$EndElements") == std::string::npos) {
			// Node data
			long id;
			double x, y, z;
			in >> n_nodes >> std::ws;
			for (size_t i = 0; i < n_nodes; i++) {
				in >> id >> x >> y >> z >> std::ws;
				mesh->nod_vector.push_back(new MeshLib::CNode(id, x, y, z));
			}
			getline(in, line); // End Node keyword $EndNodes

			// Element data
			getline(in, line); // Element keyword $Elements
			in >> n_elements >> std::ws; // number-of-elements
			for (size_t i = 0; i < n_elements; i++) {
				MeshLib::CElem* elem(new MeshLib::CElem(i));
				elem->Read(in, 7);
				if (elem->GetElementType() != MshElemType::INVALID) mesh->ele_vector.push_back(elem);
			}
			getline(in, line); // END keyword

			// correct indices TF
			const size_t n_elements(mesh->ele_vector.size());
			for (size_t k(0); k < n_elements; k++)
				mesh->ele_vector[k]->SetIndex(k);

			// ordering nodes and closing gaps TK
			std::vector<size_t> gmsh_id;
			size_t counter(0);
			for (size_t i = 0; i < mesh->nod_vector.size(); i++) {
				const size_t diff = mesh->nod_vector[i]->GetIndex() - counter;
				if (diff == 0) {
					gmsh_id.push_back(i);
					counter++;
				} else {
					for (size_t j = 0; j < diff; j++) {
						gmsh_id.push_back(i);
						counter++;
					}
					i--;
				}
			}

			for (size_t i = 0; i < mesh->ele_vector.size(); i++)
				for (long j = 0; j < mesh->ele_vector[i]->GetVertexNumber(); j++)
					mesh->ele_vector[i]->getNodeIndices()[j]
									= gmsh_id[mesh->ele_vector[i]->GetNodeIndex(j) + 1];

			for (size_t i = 0; i < mesh->nod_vector.size(); i++)
				mesh->nod_vector[i]->SetIndex(i);
			// END OF: ordering nodes and closing gaps TK
		} /*End while*/
	}
	in.close();
}

int GMSHInterface::write(std::ostream& out)
{
	out << "// GMSH input file created by OpenGeoSys " << OGS_VERSION << " built on ";
#ifdef BUILD_TIMESTAMP
	out << BUILD_TIMESTAMP;
#endif
	out << std::endl << std::endl;

	writeGMSHInputFile(out);
	return 1;
}

void GMSHInterface::writeGMSHInputFile(std::ostream& out)
{
#ifndef NDEBUG
	std::cerr << "[GMSHInterface::writeGMSHInputFile] get data from GEOObjects ... " << std::flush;
#endif
	// *** get and merge data from _geo_objs
	_gmsh_geo_name = "GMSHGeometry";
	_geo_objs.mergeGeometries(_selected_geometries, _gmsh_geo_name);
	std::vector<GEOLIB::Point*> * merged_pnts(const_cast<std::vector<GEOLIB::Point*> *>(_geo_objs.getPointVec(_gmsh_geo_name)));
	if (! merged_pnts) {
		std::cerr << "[GMSHInterface::writeGMSHInputFile] did not found any points" << std::endl;
		return;
	} else {
		const size_t n_pnts(merged_pnts->size());
		for (size_t k(0); k<n_pnts; k++) {
			(*((*merged_pnts)[k]))[2] = 0.0;
		}
	}
	std::vector<GEOLIB::Polyline*> const* merged_plys(_geo_objs.getPolylineVec(_gmsh_geo_name));
#ifndef NDEBUG
	std::cerr << "ok" << std::endl;
#endif

	// *** compute topological hierarchy of polygons
	if (merged_plys) {
		for (std::vector<GEOLIB::Polyline*>::const_iterator it(merged_plys->begin());
			it!=merged_plys->end(); it++) {
			if ((*it)->isClosed()) {
				_polygon_tree_list.push_back(new GMSHPolygonTree(new GEOLIB::Polygon(*(*it), true), NULL, _geo_objs, _gmsh_geo_name, _mesh_density_strategy));
			}
		}
		std::cout << "[GMSHInterface::writeGMSHInputFile] compute topological hierarchy - detected "
						<< _polygon_tree_list.size() << " polygons" << std::endl;
		GEOLIB::createPolygonTrees<FileIO::GMSHPolygonTree>(_polygon_tree_list);
		std::cout << "[GMSHInterface::writeGMSHInputFile] compute topological hierarchy - calculated "
								<< _polygon_tree_list.size() << " polygon trees" << std::endl;
	} else {
		return;
	}

	// *** insert stations and polylines (except polygons) in the appropriate object of
	//     class GMSHPolygonTree
	// *** insert stations
	const size_t n_geo_names(_selected_geometries.size());
	for (size_t j(0); j < n_geo_names; j++) {
		const std::vector<GEOLIB::Point*>* stations (_geo_objs.getStationVec(_selected_geometries[j]));
		if (stations) {
			const size_t n_stations(stations->size());
			for (size_t k(0); k < n_stations; k++) {
				bool found(false);
				for (std::list<GMSHPolygonTree*>::iterator it(_polygon_tree_list.begin());
					it != _polygon_tree_list.end() && !found; it++) {
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
			for (std::list<GMSHPolygonTree*>::iterator it(_polygon_tree_list.begin());
				it != _polygon_tree_list.end(); it++) {
				(*it)->insertPolyline(new GEOLIB::PolylineWithSegmentMarker(*(*merged_plys)[k]));
			}
		}
	}

	// *** init mesh density strategies
	for (std::list<GMSHPolygonTree*>::iterator it(_polygon_tree_list.begin());
		it != _polygon_tree_list.end(); it++) {
		(*it)->initMeshDensityStrategy();
	}

	// *** create GMSH data structures
	const size_t n_merged_pnts(merged_pnts->size());
	_gmsh_pnts.resize(n_merged_pnts);
	for (size_t k(0); k<n_merged_pnts; k++) {
		_gmsh_pnts[k] = NULL;
	}
	for (std::list<GMSHPolygonTree*>::iterator it(_polygon_tree_list.begin());
		it != _polygon_tree_list.end(); it++) {
		(*it)->createGMSHPoints(_gmsh_pnts);
	}

	// *** finally write data :-)
	writePoints(out);
	size_t pnt_id_offset(_gmsh_pnts.size());
	for (std::list<GMSHPolygonTree*>::iterator it(_polygon_tree_list.begin());
		it != _polygon_tree_list.end(); it++) {
		(*it)->writeLineLoop(_n_lines, _n_plane_sfc, out);
		(*it)->writeSubPolygonsAsLineConstraints(_n_lines, _n_plane_sfc-1, out);
		(*it)->writeLineConstraints(_n_lines, _n_plane_sfc-1, out);
		(*it)->writeStations(pnt_id_offset, _n_plane_sfc-1, out);
		(*it)->writeAdditionalPointData(pnt_id_offset, _n_plane_sfc-1, out);
	}

	_geo_objs.removeSurfaceVec(_gmsh_geo_name);
	_geo_objs.removePolylineVec(_gmsh_geo_name);
	_geo_objs.removePointVec(_gmsh_geo_name);
}

void GMSHInterface::writePoints(std::ostream& out) const
{
	const size_t n_gmsh_pnts(_gmsh_pnts.size());
	for (size_t k(0); k<n_gmsh_pnts; k++) {
		if (_gmsh_pnts[k]) {
			out << *(_gmsh_pnts[k]) << std::endl;
		}
	}
}

} // end namespace FileIO
