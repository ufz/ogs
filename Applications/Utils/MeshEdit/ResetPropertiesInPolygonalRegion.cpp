/*
 * \date 2014-09-25
 * \brief Reset material properties in meshes in a polygonal region.
 */

// TCLAP
#include "tclap/CmdLine.h"

// ThirdParty/logog
#include "logog/include/logog.hpp"

// BaseLib
#include "BaseLib/LogogSimpleFormatter.h"

// FileIO
#include "FileIO/readMeshFromFile.h"
#include "FileIO/XmlIO/Boost/BoostXmlGmlInterface.h"
#include "FileIO/VtkIO/VtuInterface.h"

// GeoLib
#include "GeoLib/GEOObjects.h"
#include "GeoLib/Polygon.h"
#include "GeoLib/AnalyticalGeometry.h"

// MathLib
#include "MathLib/Vector3.h"

// MeshLib
#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/MeshEditing/ElementValueModification.h"

GeoLib::Polygon rotatePolygonToXY(GeoLib::Polygon const& polygon_in,
	MathLib::Vector3 & plane_normal)
{
	// 1 copy all points
	std::vector<GeoLib::Point*> *polygon_pnts(new std::vector<GeoLib::Point*>);
	for (std::size_t k(0); k < polygon_in.getNumberOfPoints(); k++)
		polygon_pnts->push_back (new GeoLib::Point (*(polygon_in.getPoint(k))));

	// 2 rotate points
	double d_polygon (0.0);
	GeoLib::getNewellPlane (*polygon_pnts, plane_normal, d_polygon);
	MathLib::DenseMatrix<double> rot_mat(3,3);
	GeoLib::computeRotationMatrixToXY(plane_normal, rot_mat);
	GeoLib::rotatePoints(rot_mat, *polygon_pnts);

	// 3 set z coord to zero
	std::for_each(polygon_pnts->begin(), polygon_pnts->end(),
		[] (GeoLib::Point* p) { (*p)[2] = 0.0; }
	);

	// 4 create new polygon
	GeoLib::Polyline rot_polyline(*polygon_pnts);
	for (std::size_t k(0); k < polygon_in.getNumberOfPoints(); k++)
		rot_polyline.addPoint (k);
	rot_polyline.addPoint (0);
	return GeoLib::Polygon(rot_polyline);
}

std::vector<bool> markNodesOutSideOfPolygon(
	std::vector<MeshLib::Node*> const& nodes,
	GeoLib::Polygon const& polygon)
{
	// *** rotate polygon to xy_plane
	MathLib::Vector3 normal;
	GeoLib::Polygon rot_polygon(rotatePolygonToXY(polygon, normal));

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

	TCLAP::CmdLine cmd("desc", ' ', "0.1");
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
	TCLAP::ValueArg<unsigned> polygon_id_arg("p", "polygon-id",
		"id of polygon in the geometry", true, 0, "number");
	cmd.add(polygon_id_arg);
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
			ERR("Problems to read geometry from file \"%s\".",
				geometry_fname.getValue().c_str());
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
	std::vector<GeoLib::Polyline*> const* plys(geometries.getPolylineVec(geo_name));
	if (!plys) {
		ERR("Could not get vector of polylines out of geometry \"%s\".",
			geo_name.c_str());
		delete mesh;
		return EXIT_FAILURE;
	}

	// *** get polygon id
	std::size_t polygon_id(polygon_id_arg.getValue());
	if (plys->size() <= polygon_id) {
		ERR("Polyline for id %d not found.", polygon_id);
		delete mesh;
		return EXIT_FAILURE;
	}

	// *** check if the polyline is closed (i.e. is a polygon)
	bool closed ((*plys)[polygon_id]->isClosed());
	if (!closed)
	{
		ERR("Polyline with id %d is not closed, i.e. does not describe a\
			region.", polygon_id);
		delete mesh;
		return EXIT_FAILURE;
	}

	std::size_t new_property(new_material_id_arg.getValue());

	GeoLib::Polygon polygon(*((*plys)[polygon_id]));

	resetProperty(*mesh, polygon, new_property);

	FileIO::VtuInterface mesh_io(mesh);
	mesh_io.writeToFile (mesh_out.getValue());

	return EXIT_SUCCESS;
}
