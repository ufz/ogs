/**
 * @file ComputeSurfaceNodeIDsInPolygonalRegion.cpp
 * @brief Computes mesh node ids of mesh nodes within a polygonal region, that resides on the surface.
 *
 * @copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include <algorithm>
#include <memory>
#include <string>
#include <vector>

#include "tclap/CmdLine.h"
#include "logog/include/logog.hpp"

#include "Applications/ApplicationsLib/LogogSetup.h"

#include "BaseLib/LogogSimpleFormatter.h"
#include "BaseLib/StringTools.h"
#include "BaseLib/FileTools.h"

#include "FileIO/readMeshFromFile.h"
#include "FileIO/readGeometryFromFile.h"

#include "GeoLib/GEOObjects.h"
#include "GeoLib/Polygon.h"

#include "MathLib/Vector3.h"

#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshLib/MeshSurfaceExtraction.h"

void writeToFile(std::string const& id_area_fname, std::string const& csv_fname,
	std::vector<std::pair<std::size_t, double>> const& ids_and_areas,
	std::vector<MeshLib::Node*> const& mesh_nodes)
{
	std::ofstream ids_and_area_out(id_area_fname);
	if (!ids_and_area_out) {
		ERR("Unable to open the file \"%s\" - aborting.", id_area_fname.c_str());
		std::abort();
	}
	std::ofstream csv_out(csv_fname);
	if (!csv_out) {
		ERR("Unable to open the file \"%s\" - aborting.", csv_fname.c_str());
		std::abort();
	}

	ids_and_area_out << std::setprecision(20);
	csv_out << std::setprecision(20);

	ids_and_area_out << ids_and_areas[0].first << " " << ids_and_areas[0].second;
	csv_out << "ID x y z area node_id\n"; // CSV header
	csv_out << 0 << " " << *mesh_nodes[ids_and_areas[0].first]
			<< ids_and_areas[0].second << " " << ids_and_areas[0].first;
	for (std::size_t k(1); k<ids_and_areas.size(); k++) {
		ids_and_area_out << "\n"
			<< ids_and_areas[k].first << " " << ids_and_areas[k].second;
		csv_out << "\n" << k << " " << *mesh_nodes[ids_and_areas[k].first]
			<< ids_and_areas[k].second << " " << ids_and_areas[k].first;
	}
	ids_and_area_out << "\n";
	csv_out << "\n";
}

int main (int argc, char* argv[])
{
	ApplicationsLib::LogogSetup logog_setup;

	TCLAP::CmdLine cmd("Computes ids of mesh nodes that are in polygonal "
		"regions and resides on the top surface. The polygonal regions have to "
		"be given in a gml- or gli-file. The found mesh nodes and the associated"
		" area are written as txt and csv data.", ' ', "0.1");
	TCLAP::ValueArg<std::string> mesh_in("m", "mesh-input-file",
		"the name of the file containing the input mesh", true,
		"", "file name of input mesh");
	cmd.add(mesh_in);
	TCLAP::ValueArg<std::string> geo_in("g", "geo-file",
		"the name of the gml file containing the polygons", true,
		"", "file name of input geometry");
	cmd.add(geo_in);

	cmd.parse(argc, argv);

	std::unique_ptr<MeshLib::Mesh const> mesh(FileIO::readMeshFromFile(mesh_in.getValue()));
	INFO("Mesh read: %ul nodes, %ul elements.", mesh->getNNodes(), mesh->getNElements());

	GeoLib::GEOObjects geo_objs;
	FileIO::readGeometryFromFile(geo_in.getValue(), geo_objs);
	std::vector<std::string> geo_names;
	geo_objs.getGeometryNames(geo_names);
	INFO("Geometry \"%s\" read: %ul points, %ul polylines.",
		geo_names[0].c_str(),
		geo_objs.getPointVec(geo_names[0])->size(),
		geo_objs.getPolylineVec(geo_names[0])->size());

	MathLib::Vector3 const dir(0.0, 0.0, -1.0);
	double angle(90);

	auto computeElementTopSurfaceAreas = [](MeshLib::Mesh const& mesh,
		MathLib::Vector3 const& d, double angle)
	{
		std::unique_ptr<MeshLib::Mesh> surface_mesh(
			MeshLib::MeshSurfaceExtraction::getMeshSurface(mesh, d, angle));
		return MeshLib::MeshSurfaceExtraction::getSurfaceAreaForNodes(
			*surface_mesh.get());
	};

	std::vector<double> areas(computeElementTopSurfaceAreas(*mesh, dir, angle));
	std::vector<GeoLib::Point*> all_sfc_pnts(
		MeshLib::MeshSurfaceExtraction::getSurfaceNodes(*mesh, dir, angle)
	);

	std::for_each(all_sfc_pnts.begin(), all_sfc_pnts.end(), [](GeoLib::Point* p) { (*p)[2] = 0.0; });

	std::vector<MeshLib::Node*> const& mesh_nodes(mesh->getNodes());
	GeoLib::PolylineVec const* ply_vec(
		geo_objs.getPolylineVecObj(geo_names[0])
	);
	std::vector<GeoLib::Polyline*> const& plys(*(ply_vec->getVector()));

	for (std::size_t j(0); j<plys.size(); j++) {
		if (! plys[j]->isClosed()) {
			continue;
		}
		std::string polygon_name;
		ply_vec->getNameOfElement(plys[j], polygon_name);
		if (polygon_name.empty())
			polygon_name = "Polygon-" + std::to_string(j);
		// create Polygon from Polyline
		GeoLib::Polygon const& polygon(*(plys[j]));
		// ids of mesh nodes on surface that are within the given polygon
		std::vector<std::pair<std::size_t, double>> ids_and_areas;
		for (std::size_t k(0); k<all_sfc_pnts.size(); k++) {
			GeoLib::Point const& pnt(*(all_sfc_pnts[k]));
			if (polygon.isPntInPolygon(pnt)) {
				ids_and_areas.push_back(std::make_pair(pnt.getID(), areas[k]));
			}
		}
		if (ids_and_areas.empty()) {
			ERR("Polygonal part of surface \"%s\" doesn't contains nodes. No "
				"output will be generated.", polygon_name.c_str());
			continue;
		}

		std::string const out_path(BaseLib::extractPath(geo_in.getValue()));
		std::string id_and_area_fname(out_path + polygon_name);
		std::string csv_fname(out_path + polygon_name);
		id_and_area_fname += std::to_string(j) + ".txt";
		csv_fname += std::to_string(j) + ".csv";
		INFO("Polygonal part of surface \"%s\" contains %ul nodes. Writting to"
			" files \"%s\" and \"%s\".",
			polygon_name.c_str(),
			ids_and_areas.size(),
			id_and_area_fname.c_str(),
			csv_fname.c_str()
		);
		writeToFile(id_and_area_fname, csv_fname, ids_and_areas, mesh_nodes);
	}

	std::for_each(all_sfc_pnts.begin(), all_sfc_pnts.end(), std::default_delete<GeoLib::Point>());

	return 0;
}
