/*
 * \date 2014-09-25
 * \brief Reset material properties in meshes in a polygonal region.
 */

#include <algorithm>
#include <cstdlib>
#include <vector>

// TCLAP
#include "tclap/CmdLine.h"

// ThirdParty/logog
#include "logog/include/logog.hpp"

// BaseLib
#include "BaseLib/LogogSimpleFormatter.h"

// FileIO
#include "FileIO/readMeshFromFile.h"
#include "FileIO/writeMeshToFile.h"
#include "FileIO/XmlIO/Boost/BoostXmlGmlInterface.h"

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
#include "MeshLib/MeshEditing/ElementValueModification.h"

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
		rotated_nodes.push_back(new GeoLib::Point(nodes[j]->getCoords()));
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

void resetProperty(MeshLib::Mesh &mesh, GeoLib::Polygon const& polygon,
	std::size_t new_property)
{
	std::vector<bool> outside(markNodesOutSideOfPolygon(mesh.getNodes(),
		polygon));

	for(MeshLib::Element * elem : mesh.getElements()) {
		bool elem_out(true);
		for(std::size_t k(0); k<elem->getNNodes() && elem_out; ++k) {
			if (! outside[elem->getNode(k)->getID()]) {
				elem_out = false;
			}
		}
		if (elem_out) {
			elem->setValue(new_property);
		}
	}
}

int main (int argc, char* argv[])
{
	LOGOG_INITIALIZE();
	logog::Cout* logog_cout (new logog::Cout);
	BaseLib::LogogSimpleFormatter *custom_format (new BaseLib::LogogSimpleFormatter);
	logog_cout->SetFormatter(*custom_format);

	TCLAP::CmdLine cmd("Sets the property id of an mesh element to a given new "
		"id iff at least one node of the element is within a polygonal region "
		"that is given by a polygon read from a gml file.", ' ', "0.1");
	TCLAP::ValueArg<std::string> mesh_in("i", "mesh-input-file",
		"the name of the file containing the input mesh", true,
		"", "file name");
	cmd.add(mesh_in);
	TCLAP::ValueArg<std::string> mesh_out("o", "mesh-output-file",
		"the name of the file the mesh will be written to", true,
		"", "file name");
	cmd.add(mesh_out);
	TCLAP::ValueArg<std::string> geometry_fname("g", "geometry",
		"the name of the file containing the input geometry", true,
		"", "file name");
	cmd.add(geometry_fname);
	TCLAP::ValueArg<std::string> polygon_name_arg("p", "polygon-name",
		"name of polygon in the geometry", true, "", "string");
	cmd.add(polygon_name_arg);
	TCLAP::ValueArg<unsigned> new_material_id_arg("n", "new-material-id",
		"new material id", false, 0, "number");
	cmd.add(new_material_id_arg);
	cmd.parse(argc, argv);

	// *** read mesh
	MeshLib::Mesh * mesh(FileIO::readMeshFromFile(mesh_in.getValue()));

	// *** read geometry
	GeoLib::GEOObjects geometries;
	{
		FileIO::BoostXmlGmlInterface xml_io(geometries);
		if (xml_io.readFile(geometry_fname.getValue())) {
			INFO("Read geometry from file \"%s\".",
				geometry_fname.getValue().c_str());
		} else {
			delete mesh;
			return EXIT_FAILURE;
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
	GeoLib::PolylineVec const* plys(geometries.getPolylineVecObj(geo_name));
	if (!plys) {
		ERR("Could not get vector of polylines out of geometry \"%s\".",
			geo_name.c_str());
		delete mesh;
		return EXIT_FAILURE;
	}

	// *** get polygon
	GeoLib::Polyline const* ply(
		plys->getElementByName(polygon_name_arg.getValue())
	);
	if (! ply) {
		ERR("Polyline \"%s\" not found.", polygon_name_arg.getValue().c_str());
		delete mesh;
		return EXIT_FAILURE;
	}

	// *** check if the polyline is closed (i.e. is a polygon)
	bool closed (ply->isClosed());
	if (!closed)
	{
		ERR("Polyline \"%s\" is not closed, i.e. does not describe a\
			region.", polygon_name_arg.getValue().c_str());
		delete mesh;
		return EXIT_FAILURE;
	}

	std::size_t new_property(new_material_id_arg.getValue());

	GeoLib::Polygon polygon(*(ply));

	resetProperty(*mesh, polygon, new_property);

	FileIO::writeMeshToFile(*mesh, mesh_out.getValue());

	return EXIT_SUCCESS;
}
