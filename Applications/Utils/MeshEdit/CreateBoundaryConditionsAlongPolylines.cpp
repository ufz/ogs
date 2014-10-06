/*
 * \date 2014-09-30
 * \brief Create BoundaryConditions from a polylines.
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <map>
#include <string>
#include <vector>
#include <fstream>

// TCLAP
#include "tclap/CmdLine.h"

// ThirdParty/logog
#include "logog/include/logog.hpp"

// BaseLib
#include "BaseLib/LogogSimpleFormatter.h"

// FileIO
#include "FileIO/readMeshFromFile.h"
#include "FileIO/XmlIO/Boost/BoostXmlGmlInterface.h"
#include "FileIO/XmlIO/Boost/BoostVtuInterface.h"

// GeoLib
#include "GeoLib/GEOObjects.h"

// MeshLib
#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"

// MeshGeoToolsLib
#include "MeshGeoToolsLib/MeshNodeSearcher.h"

void convertMeshNodesToGeometry(std::vector<MeshLib::Node*> const& nodes,
	std::vector<std::size_t> const& node_ids,
	std::string & geo_name,
	GeoLib::GEOObjects & geometry_sets)
{
	// copy data
	std::vector<GeoLib::Point*> * pnts(new std::vector<GeoLib::Point*>);
	std::map<std::string, std::size_t>* pnt_names(
		new std::map<std::string, std::size_t>);
	std::size_t cnt(0);
	for (std::size_t id: node_ids) {
		pnts->push_back(new GeoLib::Point(nodes[id]->getCoords()));
		pnt_names->insert(std::pair<std::string, std::size_t>(
			geo_name+"-PNT-"+std::to_string(cnt), cnt));
		cnt++;
	}

	// create data structures for geometry
	geometry_sets.addPointVec(pnts, geo_name, pnt_names);
}

void writeBCsAndGML(GeoLib::GEOObjects & geometry_sets,
	std::string & geo_name)
{
	INFO("write points to \"%s.gml\".", geo_name.c_str());
	FileIO::BoostXmlGmlInterface xml_io(geometry_sets);
	xml_io.setNameForExport(geo_name);
	xml_io.writeToFile(geo_name+".gml");

	GeoLib::PointVec const* pnt_vec_objs(geometry_sets.getPointVecObj(geo_name));
	std::vector<GeoLib::Point*> const& pnts(*(pnt_vec_objs->getVector()));
	std::string fname("UnstrutCatchment.gli");
	std::ofstream out (fname.c_str());
	out << "#POINTS\n";
	out.precision(20);
	std::ofstream bc_out ("UnstrutCatchment.bc");
	for (std::size_t k(0); k<pnts.size(); k++) {
		out << k << " " << *(pnts[k]);
		std::string pnt_name;
		if (pnt_vec_objs->getNameOfElementByID(k, pnt_name)) {
			out << "$NAME " << pnt_name;
			bc_out << "#BOUNDARY_CONDITION\n";
			bc_out << "  $PCS_TYPE\n";
			bc_out << "    LIQUID_FLOW\n";
			bc_out << "  $PRIMARY_VARIABLE\n";
			bc_out << "    PRESSURE1\n";
			bc_out << "  $GEO_TYPE\n";
			bc_out << "    POINT " << pnt_name << "\n";
			bc_out << "  $DIS_TYPE\n";
			bc_out << "    CONSTANT 0.0\n";
		}
		out << "\n";
	}
	out << "#STOP\n";
	bc_out << "#STOP\n";
	bc_out.close();
	out.close();
}

int main (int argc, char* argv[])
{
	LOGOG_INITIALIZE();
	logog::Cout* logog_cout (new logog::Cout);
	BaseLib::LogogSimpleFormatter *custom_format (new BaseLib::LogogSimpleFormatter);
	logog_cout->SetFormatter(*custom_format);

	TCLAP::CmdLine cmd(
		"Creates boundary conditions for mesh nodes along polylines.",
		' ',
		"0.1");
	TCLAP::ValueArg<std::string> mesh_arg("m", "mesh-file",
		"the name of the file containing the mesh", true,
		"", "file name");
	cmd.add(mesh_arg);
	TCLAP::ValueArg<std::string> geometry_fname("g", "geometry",
		"the name of the file containing the input geometry", true,
		"", "file name");
	cmd.add(geometry_fname);
	TCLAP::ValueArg<bool> vis_arg("v", "visualize",
		"write gml file with found mesh nodes for visualization", false, 0, "bool");
	cmd.add(vis_arg);
	TCLAP::ValueArg<double> search_length_arg("s", "search-length",
		"The size of the search length. The default value is "
		"std::numeric_limits<double>::epsilon()", false,
		std::numeric_limits<double>::epsilon(), "floating point number");
	cmd.add(search_length_arg);
	cmd.parse(argc, argv);

	// *** read mesh
	INFO("Reading mesh \"%s\" ... ", mesh_arg.getValue().c_str());
	MeshLib::Mesh * mesh(FileIO::readMeshFromFile(mesh_arg.getValue()));
	INFO("done.");

	// *** read geometry
	GeoLib::GEOObjects geometries;
	{
		FileIO::BoostXmlGmlInterface xml_io(geometries);
		if (xml_io.readFile(geometry_fname.getValue())) {
			INFO("Read geometry from file \"%s\".",
				geometry_fname.getValue().c_str());
		} else {
			ERR("Problems to read geometry from file \"%s\".",
				geometry_fname.getValue().c_str());
			delete mesh;
			return -1;
		}
	}

	std::string geo_name;
	{
		std::vector<std::string> geo_names;
		geometries.getGeometryNames(geo_names);
		geo_name = geo_names[0];
	}

	// *** check if the data is usable
	// *** get vector of polylines
	std::vector<GeoLib::Polyline*> const* plys(geometries.getPolylineVec(geo_name));
	if (!plys) {
		ERR("Could not get vector of polylines out of geometry \"%s\".",
			geo_name.c_str());
		delete mesh;
		return -1;
	}

	MeshGeoToolsLib::SearchLength search_length_strategy;
	if (search_length_arg.isSet()) {
		search_length_strategy =
			MeshGeoToolsLib::SearchLength(search_length_arg.getValue());
	}

	GeoLib::GEOObjects geometry_sets;
	MeshGeoToolsLib::MeshNodeSearcher mesh_searcher(*mesh,
		search_length_strategy);
	for(std::size_t k(0); k<plys->size(); k++) {
		std::vector<std::size_t> ids
			(mesh_searcher.getMeshNodeIDsAlongPolyline(*((*plys)[k])));
		if (ids.empty())
			continue;
		std::string geo_name("Polyline-"+std::to_string(k));
		convertMeshNodesToGeometry(mesh->getNodes(), ids, geo_name,
			geometry_sets);
	}

	// merge all together
	std::vector<std::string> geo_names;
	std::string merge_name("TM");
	geometry_sets.getGeometryNames(geo_names);
	geometry_sets.mergeGeometries(geo_names, merge_name);

	// write the BCs and the merged geometry set to file
	writeBCsAndGML(geometry_sets, merge_name);

	return 0;
}
