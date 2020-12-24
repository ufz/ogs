/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
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
#include "MathLib/LinAlg/Dense/DenseMatrix.h"
#include "MathLib/MathTools.h"
#include "MathLib/Point3d.h"
#include "MathLib/Vector3.h"
#include "MeshGeoToolsLib/GeoMapper.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/IO/VtkIO/VtuInterface.h"
#include "MeshLib/IO/readMeshFromFile.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshEditing/RemoveMeshComponents.h"
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

    std::size_t const n_lines = lines->size();
    geo.addPointVec(std::move(points), merged_geo_name);
    geo.addPolylineVec(std::move(lines), merged_geo_name);
}

/// rotates the merged geometry into the XY-plane
void rotateGeometryToXY(std::vector<GeoLib::Point*>& points,
                        MathLib::DenseMatrix<double>& rotation_matrix,
                        double& z_shift)
{
    // compute the plane normal
    auto const [plane_normal, d] =
        GeoLib::getNewellPlane(points.begin(), points.end());
    // rotate points into x-y-plane
    Eigen::Matrix3d const rotation_matrix_ =
        GeoLib::computeRotationMatrixToXY(plane_normal);
    GeoLib::rotatePoints(rotation_matrix_, points.begin(), points.end());
    // Todo (TF) Remove when rotateGeometryToXY uses Eigen for rot_mat
    for (int r = 0; r < 3; r++)
    {
        for (int c = 0; c < 3; c++)
        {
            rotation_matrix(r, c) = rotation_matrix_(r, c);
        }
    }

    GeoLib::AABB aabb(points.begin(), points.end());
    z_shift = (aabb.getMinPoint()[2] + aabb.getMaxPoint()[2]) / 2.0;
    std::for_each(points.begin(), points.end(),
                  [z_shift](GeoLib::Point* p) { (*p)[2] -= z_shift; });
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
    xml.setNameForExport(merged_geo_name);
    xml.writeToFile(filename);

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
    gmsh_io.setPrecision(std::numeric_limits<double>::digits10);
    bool const success = gmsh_io.writeToFile(gmsh_geo_name);

    std::string const gmsh_mesh_name = output_name + ".msh";
    std::string gmsh_command = "gmsh -2 -algo meshadapt " + gmsh_geo_name;
    gmsh_command += " -o " + gmsh_mesh_name + " -format msh22";
    int const return_value = std::system(gmsh_command.c_str());
    if (return_value != 0)
    {
        ERR("Execution of gmsh command returned non-zero "
            "status, %d",
            return_value);
    }
    return FileIO::GMSH::readGMSHMesh(gmsh_mesh_name);
}

/// inverse rotation of the mesh, back into original position
void rotateMesh(MeshLib::Mesh& mesh,
                MathLib::DenseMatrix<double> const& rot_mat,
                double const z_shift)
{
    std::vector<MeshLib::Node*> const& nodes = mesh.getNodes();
    std::for_each(nodes.begin(), nodes.end(),
                  [z_shift](MeshLib::Node* n) { (*n)[2] += z_shift; });
    Eigen::Matrix3d rot_mat_eigen;
    // Todo (TF) Remove when rotateMesh uses Eigen for rot_mat
    for (int r = 0; r < 3; r++)
    {
        for (int c = 0; c < 3; c++)
        {
            rot_mat_eigen(r, c) = rot_mat(r, c);
        }
    }
    GeoLib::rotatePoints(rot_mat_eigen.transpose(), nodes.begin(), nodes.end());
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
            "Copyright (c) 2012-2020, OpenGeoSys Community "
            "(http://www.opengeosys.org)",
        ' ', GitInfoLib::GitInfo::ogs_version);
    TCLAP::SwitchArg test_arg("t", "testdata", "keep test data", false);
    cmd.add(test_arg);
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
    MathLib::DenseMatrix<double> rot_mat(3, 3);
    double z_shift(0);
    rotateGeometryToXY(points, rot_mat, z_shift);
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

    MeshLib::IO::VtuInterface vtu(new_mesh.get());
    vtu.writeToFile(output_name + ".vtu");
    return EXIT_SUCCESS;
}
