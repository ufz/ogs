/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <algorithm>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

// ThirdParty
#include <tclap/CmdLine.h>

#include <QCoreApplication>

#include "Applications/FileIO/Gmsh/GMSHInterface.h"
#include "Applications/FileIO/Gmsh/GmshReader.h"
#include "BaseLib/FileTools.h"
#include "BaseLib/IO/readStringListFromFile.h"
#include "GeoLib/AABB.h"
#include "GeoLib/AnalyticalGeometry.h"
#include "GeoLib/GEOObjects.h"
#include "GeoLib/IO/XmlIO/Qt/XmlGmlInterface.h"
#include "GeoLib/Point.h"
#include "GeoLib/Polygon.h"
#include "GeoLib/Polyline.h"
#include "InfoLib/GitInfo.h"
#include "MathLib/MathTools.h"
#include "MathLib/Point3d.h"
#include "MeshGeoToolsLib/GeoMapper.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/IO/VtkIO/VtuInterface.h"
#include "MeshLib/IO/readMeshFromFile.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshEditing/MeshRevision.h"
#include "MeshLib/MeshEditing/RemoveMeshComponents.h"
#include "MeshLib/MeshSearch/ElementSearch.h"
#include "MeshLib/MeshSurfaceExtraction.h"
#include "MeshLib/Node.h"

/// creates a vector of sampling points based on the specified resolution
std::unique_ptr<std::vector<GeoLib::Point*>> createPoints(
    MathLib::Point3d const& start,
    MathLib::Point3d const& end,
    std::size_t const n_intervals)
{
    auto points = std::make_unique<std::vector<GeoLib::Point*>>();
    points->push_back(new GeoLib::Point(start, 0));
    double const length = std::sqrt(MathLib::sqrDist(start, end));
    double const interval = length / n_intervals;

    GeoLib::Point const vec((end[0] - start[0]), (end[1] - start[1]),
                            (end[2] - start[2]));
    for (std::size_t i = 1; i < n_intervals; ++i)
    {
        points->push_back(new GeoLib::Point(
            start[0] + ((i * interval) / length * vec[0]),
            start[1] + ((i * interval) / length * vec[1]), 0, i));
    }
    points->push_back(new GeoLib::Point(end, n_intervals));
    return points;
}

/// creates a polyline to be mapped on a mesh layer
GeoLib::Polyline* createPolyline(std::vector<GeoLib::Point*> const& points)
{
    GeoLib::Polyline* line = new GeoLib::Polyline(points);
    std::size_t const length = points.size();
    for (std::size_t i = 0; i < length; ++i)
    {
        line->addPoint(i);
    }
    return line;
}

/// creates a mapped line for each of the mesh layers
std::vector<std::string> createGeometries(
    GeoLib::GEOObjects& geo,
    std::vector<std::string> const& layer_names,
    MathLib::Point3d const& pnt_start,
    MathLib::Point3d const& pnt_end,
    double const resolution)
{
    std::vector<std::string> geo_name_list;
    std::size_t const n_layers = layer_names.size();
    for (std::size_t i = 0; i < n_layers; ++i)
    {
        std::unique_ptr<MeshLib::Mesh> const layer(
            MeshLib::IO::readMeshFromFile(layer_names[i]));
        if (layer == nullptr)
        {
            ERR("Could not read file {:s}. Skipping layer...", layer_names[i]);
            continue;
        }
        if (layer->getDimension() != 2)
        {
            ERR("Layer {:d} is not a 2D mesh. Skipping layer...", i);
            continue;
        }

        std::string geo_name(std::to_string(i));
        std::unique_ptr<std::vector<GeoLib::Point*>> points(
            createPoints(pnt_start, pnt_end, resolution));
        geo.addPointVec(std::move(points), geo_name);

        auto lines = std::make_unique<std::vector<GeoLib::Polyline*>>();
        lines->push_back(createPolyline(*geo.getPointVec(geo_name)));
        geo.addPolylineVec(std::move(lines), geo_name);

        MeshGeoToolsLib::GeoMapper mapper(geo, geo_name);
        mapper.mapOnMesh(layer.get());
        geo_name_list.push_back(geo_name);
    }
    return geo_name_list;
}

/// Merges all layer geometries into one. Each layer is specified by one
/// polygon. This will ensure that GMSH deals with the subsequent assignment of
/// material groups automatically.
void mergeGeometries(GeoLib::GEOObjects& geo,
                     std::vector<std::string> const& geo_names,
                     std::string& merged_geo_name)
{
    auto points = std::make_unique<std::vector<GeoLib::Point*>>();
    auto lines = std::make_unique<std::vector<GeoLib::Polyline*>>();

    auto layer_pnts = *geo.getPointVec(geo_names[0]);
    std::size_t const pnts_per_line = layer_pnts.size();
    std::size_t const n_layers = geo_names.size();
    std::vector<std::size_t> last_line_idx(pnts_per_line, 0);

    for (std::size_t i = 0; i < pnts_per_line; ++i)
    {
        std::size_t const idx = pnts_per_line - i - 1;
        points->push_back(new GeoLib::Point(*layer_pnts[i], idx));
        last_line_idx[i] = idx;
    }
    for (std::size_t j = 1; j < n_layers; ++j)
    {
        GeoLib::Polyline* line = new GeoLib::Polyline(*points);
        for (std::size_t i = 0; i < pnts_per_line; ++i)
        {
            line->addPoint(last_line_idx[i]);
        }
        layer_pnts = *geo.getPointVec(geo_names[j]);
        for (std::size_t i = 0; i < pnts_per_line; ++i)
        {
            // check if for current point the lower layer boundary is actually
            // located below the upper boundary
            std::size_t idx = last_line_idx[pnts_per_line - i - 1];
            if ((*(*points)[idx])[2] > (*layer_pnts[i])[2])
            {
                idx = points->size();
                points->push_back(new GeoLib::Point(*layer_pnts[i], idx));
                last_line_idx[pnts_per_line - i - 1] = idx;
            }
            line->addPoint(idx);
        }
        // close polygon
        line->addPoint(line->getPointID(0));
        lines->push_back(line);
    }

    geo.addPointVec(std::move(points), merged_geo_name);
    geo.addPolylineVec(std::move(lines), merged_geo_name);
}

/// rotates the merged geometry into the XY-plane
std::pair<Eigen::Matrix3d, double> rotateGeometryToXY(
    std::vector<GeoLib::Point*>& points)
{
    // compute the plane normal
    auto const [plane_normal, d] =
        GeoLib::getNewellPlane(points.begin(), points.end());
    // rotate points into x-y-plane
    Eigen::Matrix3d const rotation_matrix =
        GeoLib::computeRotationMatrixToXY(plane_normal);
    GeoLib::rotatePoints(rotation_matrix, points.begin(), points.end());

    GeoLib::AABB aabb(points.begin(), points.end());
    double const z_shift =
        (aabb.getMinPoint()[2] + aabb.getMaxPoint()[2]) / 2.0;
    std::for_each(points.begin(), points.end(),
                  [z_shift](GeoLib::Point* p) { (*p)[2] -= z_shift; });
    return {rotation_matrix, z_shift};
}

/// This encapsulates a workaround:
/// For unknown reasons, the GML->GEO converter will not work correctly when
/// inputting the merged geometry directly. However, if the geometry is saved to
/// a file and immedeately loaded again, everything works fine.
void consolidateGeometry(GeoLib::GEOObjects& geo,
                         std::string const& output_name,
                         std::string& merged_geo_name,
                         bool const keep_gml_file)
{
    std::string const filename(output_name + ".gml");
    GeoLib::IO::XmlGmlInterface xml(geo);
    xml.export_name = merged_geo_name;
    BaseLib::IO::writeStringToFile(xml.writeToString(), filename);

    geo.removePolylineVec(merged_geo_name);
    geo.removePointVec(merged_geo_name);

    xml.readFile(filename);

    if (!keep_gml_file)
    {
        BaseLib::removeFile(filename);
        BaseLib::removeFile(filename + ".md5");
    }
}

/// converts geometry into GMSH format and creates mesh
MeshLib::Mesh* generateMesh(GeoLib::GEOObjects& geo,
                            std::string const& geo_name,
                            std::string const& output_name, double res)
{
    std::string const gmsh_geo_name(output_name + ".geo");
    std::vector<std::string> gmsh_geo;
    gmsh_geo.push_back(geo_name);
    FileIO::GMSH::GMSHInterface gmsh_io(
        geo, true, FileIO::GMSH::MeshDensityAlgorithm::FixedMeshDensity, res, 0,
        0, gmsh_geo, false, false);
    gmsh_io.writePhysicalGroups(true);
    if (!BaseLib::IO::writeStringToFile(gmsh_io.writeToString(), gmsh_geo_name))
    {
        ERR("Writing gmsh geo file '{:s}' failed.", gmsh_geo_name);
    }

    std::string const gmsh_mesh_name = output_name + ".msh";
    std::string gmsh_command = "gmsh -2 -algo meshadapt " + gmsh_geo_name;
    gmsh_command += " -o " + gmsh_mesh_name + " -format msh22";
    int const return_value = std::system(gmsh_command.c_str());
    if (return_value != 0)
    {
        ERR("Execution of gmsh command returned non-zero status, %d",
            return_value);
    }
    return FileIO::GMSH::readGMSHMesh(gmsh_mesh_name);
}

/// inverse rotation of the mesh, back into original position
void rotateMesh(MeshLib::Mesh& mesh, Eigen::Matrix3d const& rot_mat,
                double const z_shift)
{
    std::vector<MeshLib::Node*> const& nodes = mesh.getNodes();
    std::for_each(nodes.begin(), nodes.end(),
                  [z_shift](MeshLib::Node* n) { (*n)[2] += z_shift; });
    GeoLib::rotatePoints(rot_mat.transpose(), nodes.begin(), nodes.end());
}

/// removes line elements from mesh such that only triangles remain
MeshLib::Mesh* removeLineElements(MeshLib::Mesh const& mesh)
{
    std::vector<std::size_t> line_idx;
    std::vector<MeshLib::Element*> const& elems = mesh.getElements();
    std::for_each(elems.begin(), elems.end(), [&](auto e) {
        if (e->getGeomType() == MeshLib::MeshElemType::LINE)
        {
            line_idx.push_back(e->getID());
        }
    });
    if (line_idx.size() == mesh.getNumberOfElements())
    {
        return nullptr;
    }
    return MeshLib::removeElements(mesh, line_idx, "mesh");
}

void writeBoundary(MeshLib::Mesh const& mesh,
                   std::vector<std::size_t> const& idx_array,
                   std::string const& file_name)
{
    std::unique_ptr<MeshLib::Mesh> const boundary(
        MeshLib::removeElements(mesh, idx_array, "mesh"));
    if (boundary == nullptr)
    {
        ERR("Error extracting boundary '{:s}'", file_name);
        return;
    }
    MeshLib::IO::VtuInterface vtu(boundary.get());
    vtu.writeToFile(file_name + ".vtu");
}

void extractBoundaries(MeshLib::Mesh const& mesh,
                       std::string const& output_name,
                       MathLib::Point3d const& pnt_start)
{
    double const eps = mesh.getMinEdgeLength() / 100.0;
    std::unique_ptr<MeshLib::Mesh> boundary_mesh(
        MeshLib::BoundaryExtraction::getBoundaryElementsAsMesh(
            mesh, "bulk_node_ids", "bulk_element_ids", "bulk_face_ids"));

    auto const& elems = boundary_mesh->getElements();
    std::vector<std::size_t> left_bound_idx, right_bound_idx, top_bound_idx,
        bottom_bound_idx;
    Eigen::Vector2d const anchor(pnt_start[0], pnt_start[1]);
    for (auto e : elems)
    {
        Eigen::Vector2d const n1((*e->getNode(0))[0], (*e->getNode(0))[1]);
        Eigen::Vector2d const n2((*e->getNode(1))[0], (*e->getNode(1))[1]);
        Eigen::Vector2d const dist1(n1 - anchor);
        Eigen::Vector2d const dist2(n2 - anchor);
        std::size_t const id = e->getID();
        // elements located at left or right side
        if ((dist1 - dist2).squaredNorm() < eps)
        {
            top_bound_idx.push_back(id);
            bottom_bound_idx.push_back(id);

            if (dist1.squaredNorm() < eps)
            {
                right_bound_idx.push_back(id);
            }
            else
            {
                left_bound_idx.push_back(id);
            }
            continue;
        }
        // elements located at top or bottom
        if (dist2.squaredNorm() < dist1.squaredNorm())
        {
            top_bound_idx.push_back(id);
            left_bound_idx.push_back(id);
            right_bound_idx.push_back(id);
        }
        else
        {
            bottom_bound_idx.push_back(id);
            left_bound_idx.push_back(id);
            right_bound_idx.push_back(id);
        }
    }

    writeBoundary(*boundary_mesh, left_bound_idx, output_name + "_left");
    writeBoundary(*boundary_mesh, right_bound_idx, output_name + "_right");
    writeBoundary(*boundary_mesh, top_bound_idx, output_name + "_top");
    writeBoundary(*boundary_mesh, bottom_bound_idx, output_name + "_bottom");
}

int main(int argc, char* argv[])
{
    QCoreApplication a(argc, argv);
    TCLAP::CmdLine cmd(
        "Creates a triangle-mesh of a vertical slice out of a list of input "
        "layer meshes. The slice is defined by a start- and end-point. In "
        "addition, the resolution for meshing the extracted slice (i.e. the "
        "maximum edge length of the domain discretisation) needs to be "
        "specified. The utility requires access to the meshing utility GMSH to "
        "work correctly.\n\n"
        "OpenGeoSys-6 software, version " +
            GitInfoLib::GitInfo::ogs_version +
            ".\n"
            "Copyright (c) 2012-2021, OpenGeoSys Community "
            "(http://www.opengeosys.org)",
        ' ', GitInfoLib::GitInfo::ogs_version);
    TCLAP::SwitchArg test_arg("t", "testdata", "keep test data", false);
    cmd.add(test_arg);
    TCLAP::SwitchArg bound_arg("b", "bounds", "save mesh boundaries", false);
    cmd.add(bound_arg);
    TCLAP::ValueArg<double> res_arg(
        "r", "resolution",
        "desired edge length of triangles in the resulting slice.", true, 0,
        "floating point number");
    cmd.add(res_arg);
    TCLAP::ValueArg<double> end_y_arg(
        "", "end-y", "y-coordinates of the end point defining the slice", true,
        0, "floating point number");
    cmd.add(end_y_arg);
    TCLAP::ValueArg<double> end_x_arg(
        "", "end-x", "x-coordinates of the end point defining the slice", true,
        0, "floating point number");
    cmd.add(end_x_arg);
    TCLAP::ValueArg<double> start_y_arg(
        "", "start-y", "y-coordinates of the start point defining the slice",
        true, 0, "floating point number");
    cmd.add(start_y_arg);
    TCLAP::ValueArg<double> start_x_arg(
        "", "start-x", "x-coordinates of the start point defining the slice",
        true, 0, "floating point number");
    cmd.add(start_x_arg);
    TCLAP::ValueArg<std::string> output_arg(
        "o", "output", "name of output mesh (*.vtu)", true, "", "string");
    cmd.add(output_arg);
    TCLAP::ValueArg<std::string> input_arg(
        "i", "input",
        "name of the input file list containing the paths the all input layers "
        "in correct order from top to bottom",
        true, "", "string");
    cmd.add(input_arg);
    cmd.parse(argc, argv);

    std::string const input_name = input_arg.getValue();
    std::string const output_name = output_arg.getValue();

    MathLib::Point3d const pnt_start{
        {start_x_arg.getValue(), start_y_arg.getValue(), 0.0}};
    MathLib::Point3d const pnt_end{
        {end_x_arg.getValue(), end_y_arg.getValue(), 0.0}};
    double const length = std::sqrt(MathLib::sqrDist(pnt_start, pnt_end));

    std::size_t const res = std::ceil(length / res_arg.getValue());
    double const interval_length = length / res;

    std::vector<std::string> const layer_names =
        BaseLib::IO::readStringListFromFile(input_name);
    if (layer_names.size() < 2)
    {
        ERR("At least two layers are required to extract a slice.");
        return EXIT_FAILURE;
    }

    GeoLib::GEOObjects geo;
    std::vector<std::string> const geo_name_list =
        createGeometries(geo, layer_names, pnt_start, pnt_end, res);

    if (geo_name_list.size() < 2)
    {
        ERR("Less than two geometries could be created from layers. Aborting "
            "extraction...");
        return EXIT_FAILURE;
    }

    std::string merged_geo_name = "merged_geometries";
    mergeGeometries(geo, geo_name_list, merged_geo_name);
    std::vector<GeoLib::Point*> points = *geo.getPointVec(merged_geo_name);
    auto const [rot_mat, z_shift] = rotateGeometryToXY(points);
    consolidateGeometry(geo, output_name, merged_geo_name, test_arg.getValue());

    std::unique_ptr<MeshLib::Mesh> mesh(
        generateMesh(geo, merged_geo_name, output_name, interval_length));
    if (!test_arg.getValue())
    {
        BaseLib::removeFile(output_name + ".geo");
        BaseLib::removeFile(output_name + ".msh");
    }
    if (mesh == nullptr)
    {
        ERR("Error generating mesh... (GMSH was unable to output mesh)");
        return EXIT_FAILURE;
    }
    rotateMesh(*mesh, rot_mat, z_shift);
    std::unique_ptr<MeshLib::Mesh> new_mesh(removeLineElements(*mesh));
    if (new_mesh == nullptr)
    {
        ERR("Error generating mesh... (GMSH created line mesh)");
        return EXIT_FAILURE;
    }

    // collapse all nodes that might have been created due to gmsh physically
    // separating layers
    MeshLib::MeshRevision rev(*new_mesh);
    std::unique_ptr<MeshLib::Mesh> revised_mesh(
        rev.simplifyMesh("RevisedMesh", new_mesh->getMinEdgeLength() / 100.0));

    if (bound_arg.getValue())
    {
        extractBoundaries(*revised_mesh, output_name, pnt_start);
    }
    MeshLib::IO::VtuInterface vtu(revised_mesh.get());
    vtu.writeToFile(output_name + ".vtu");
    return EXIT_SUCCESS;
}
