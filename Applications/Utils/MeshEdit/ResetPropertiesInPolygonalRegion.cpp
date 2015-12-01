/*
 * \brief Reset material properties in meshes in a polygonal region.
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <algorithm>
#include <cstdlib>
#include <vector>

// TCLAP
#include "tclap/CmdLine.h"

// ThirdParty/logog
#include "logog/include/logog.hpp"

#include "Applications/ApplicationsLib/LogogSetup.h"

// FileIO
#include "FileIO/readMeshFromFile.h"
#include "FileIO/writeMeshToFile.h"
#include "FileIO/readGeometryFromFile.h"

// GeoLib
#include "GeoLib/GEOObjects.h"
#include "GeoLib/Polygon.h"
#include "GeoLib/AnalyticalGeometry.h"

// MathLib
#include "MathLib/Vector3.h"
#include "MathLib/LinAlg/Dense/DenseMatrix.h"

// MeshLib
#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshLib/Elements/Element.h"

std::vector<bool> markNodesOutSideOfPolygon(
	std::vector<MeshLib::Node*> const& nodes,
	GeoLib::Polygon const& polygon)
{
	// *** rotate polygon to xy_plane
	MathLib::Vector3 normal;
	GeoLib::Polygon rot_polygon(GeoLib::rotatePolygonToXY(polygon, normal));

	// *** rotate mesh nodes to xy-plane
	// 1 copy all mesh nodes to GeoLib::Points
	std::vector<GeoLib::Point*> rotated_nodes;
	for (std::size_t j(0); j < nodes.size(); j++)
		rotated_nodes.push_back(new GeoLib::Point(*nodes[j], nodes[j]->getID()));
	// 2 rotate the Points
	MathLib::DenseMatrix<double> rot_mat(3,3);
	GeoLib::computeRotationMatrixToXY(normal, rot_mat);
	GeoLib::rotatePoints(rot_mat, rotated_nodes);
	// 3 set z coord to zero
	std::for_each(rotated_nodes.begin(), rotated_nodes.end(),
		[] (GeoLib::Point* p) { (*p)[2] = 0.0; }
	);

	// *** mark rotated nodes
	std::vector<bool> outside(rotated_nodes.size(), true);
	for (std::size_t k(0); k<rotated_nodes.size(); k++) {
		if (rot_polygon.isPntInPolygon(*(rotated_nodes[k]))) {
			outside[k] = false;
		}
	}

	for (std::size_t j(0); j < rotated_nodes.size(); j++)
		delete rotated_nodes[j];

	std::vector<GeoLib::Point*> & rot_polygon_pnts(
		const_cast<std::vector<GeoLib::Point*> &>(
			rot_polygon.getPointsVec()
		)
	);
	for (std::size_t k(0); k < rot_polygon_pnts.size(); k++)
		delete rot_polygon_pnts[k];

	return outside;
}

template <typename PT>
void resetMeshElementProperty(MeshLib::Mesh &mesh, GeoLib::Polygon const& polygon,
	std::string const& property_name, PT new_property_value)
{
	boost::optional<MeshLib::PropertyVector<PT> &> opt_pv(
		mesh.getProperties().getPropertyVector<PT>(property_name)
	);
	if (!opt_pv) {
		WARN("Did not find a PropertyVector with name \"%s\".",
			property_name.c_str());
		return;
	}
	MeshLib::PropertyVector<PT> & pv(opt_pv.get());
	if (pv.getMeshItemType() != MeshLib::MeshItemType::Cell) {
		ERR("Values of the PropertyVector are not assigned to cells.");
		return;
	}

	std::vector<bool> outside(markNodesOutSideOfPolygon(mesh.getNodes(),
		polygon));

	for(std::size_t j(0); j<mesh.getElements().size(); ++j) {
		bool elem_out(true);
		MeshLib::Element const*const elem(mesh.getElements()[j]);
		for(std::size_t k(0); k<elem->getNNodes() && elem_out; ++k) {
			if (! outside[elem->getNode(k)->getID()]) {
				elem_out = false;
			}
		}
		if (elem_out) {
			pv[j] = new_property_value;
		}
	}
}

int main (int argc, char* argv[])
{
	ApplicationsLib::LogogSetup logog_setup;

	TCLAP::CmdLine cmd("Sets the property id of an mesh element to a given new "
		"id iff at least one node of the element is within a polygonal region "
		"that is given by a polygon read from a gml file.", ' ', "0.1");
	TCLAP::ValueArg<std::string> mesh_out("o", "mesh-output-file",
		"the name of the file the mesh will be written to", true,
		"", "file name");
	cmd.add(mesh_out);
	TCLAP::ValueArg<std::string> polygon_name_arg("p", "polygon-name",
		"name of polygon in the geometry", true, "", "string");
	cmd.add(polygon_name_arg);
	TCLAP::ValueArg<std::string> geometry_fname("g", "geometry",
		"the name of the file containing the input geometry", true,
		"", "file name");
	cmd.add(geometry_fname);
	TCLAP::ValueArg<char> new_property_arg("v", "new-property-value",
		"new property value", false, 'A', "number");
	cmd.add(new_property_arg);
	TCLAP::ValueArg<std::string> property_name_arg("n", "property-name",
		"name of property in the mesh", false, "MaterialIDs", "string");
	cmd.add(property_name_arg);
	TCLAP::ValueArg<std::string> mesh_in("i", "mesh-input-file",
		"the name of the file containing the input mesh", true,
		"", "file name");
	cmd.add(mesh_in);
	cmd.parse(argc, argv);

	// *** read geometry
	GeoLib::GEOObjects geometries;
	FileIO::readGeometryFromFile(geometry_fname.getValue(), geometries);

	std::string geo_name;
	{
		std::vector<std::string> geo_names;
		geometries.getGeometryNames(geo_names);
		geo_name = geo_names[0];
	}

	// *** check if the data is usable
	// *** get vector of polylines
	GeoLib::PolylineVec const* plys(geometries.getPolylineVecObj(geo_name));
	if (!plys) {
		ERR("Could not get vector of polylines out of geometry \"%s\".",
			geo_name.c_str());
		return EXIT_FAILURE;
	}

	// *** get polygon
	GeoLib::Polyline const* ply(
		plys->getElementByName(polygon_name_arg.getValue())
	);
	if (! ply) {
		ERR("Polyline \"%s\" not found.", polygon_name_arg.getValue().c_str());
		return EXIT_FAILURE;
	}

	// *** check if the polyline is closed (i.e. is a polygon)
	bool closed (ply->isClosed());
	if (!closed)
	{
		ERR("Polyline \"%s\" is not closed, i.e. does not describe a\
			region.", polygon_name_arg.getValue().c_str());
		return EXIT_FAILURE;
	}

	GeoLib::Polygon polygon(*(ply));

	// *** read mesh
	MeshLib::Mesh * mesh(FileIO::readMeshFromFile(mesh_in.getValue()));
	std::vector<std::string> property_names(
		mesh->getProperties().getPropertyVectorNames());
	INFO("Mesh contains %d property vectors:", property_names.size());
	for (auto name : property_names) {
		INFO("- %s", name.c_str());
	}
	std::string const& property_name(property_name_arg.getValue());

	{
		char new_property_val(static_cast<char>(new_property_arg.getValue()));

		// check if PropertyVector exists
		boost::optional<MeshLib::PropertyVector<char> &> opt_pv(
			mesh->getProperties().getPropertyVector<char>(property_name)
		);
		if (!opt_pv) {
			opt_pv = mesh->getProperties().createNewPropertyVector<char>(
				property_name, MeshLib::MeshItemType::Cell, 1);
			opt_pv.get().resize(mesh->getElements().size());
			INFO("Created PropertyVector with name \"%s\".", property_name.c_str());
		}
		resetMeshElementProperty(*mesh, polygon, property_name, new_property_val);
	}

	{
		unsigned new_property_val(new_property_arg.getValue());
		resetMeshElementProperty(*mesh, polygon, property_name, new_property_val);
	}

	{
		bool new_property_val(static_cast<bool>(new_property_arg.getValue()));
		resetMeshElementProperty(*mesh, polygon, property_name, new_property_val);
	}

	FileIO::writeMeshToFile(*mesh, mesh_out.getValue());

	return EXIT_SUCCESS;
}
