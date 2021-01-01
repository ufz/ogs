/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <iterator>
#include <limits>
#include <memory>
#include <string>
#include <vector>

// ThirdParty
#include <tclap/CmdLine.h>

#include "GeoLib/GEOObjects.h"
#include "GeoLib/IO/XmlIO/Boost/BoostXmlGmlInterface.h"
#include "InfoLib/GitInfo.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Elements/Line.h"
#include "MeshLib/IO/VtkIO/VtuInterface.h"
#include "MeshLib/IO/readMeshFromFile.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshEditing/DuplicateMeshComponents.h"
#include "MeshLib/Node.h"

std::vector<std::size_t> getNodes(
    GeoLib::Point const& pnt, std::vector<MeshLib::Node*> const& nodes,
    MeshLib::PropertyVector<int> const& mat_ids,
    std::pair<int, int> const& mat_limits,
    std::pair<double, double> const& elevation_limits)
{
    std::vector<std::size_t> pnt_nodes;
    for (auto node : nodes)
    {
        double const eps = std::numeric_limits<double>::epsilon();
        if (std::abs((*node)[0] - pnt[0]) < eps &&
            std::abs((*node)[1] - pnt[1]) < eps)
        {
            auto const& elems = node->getElements();
            for (auto e : elems)
            {
                if (mat_ids[e->getID()] >= mat_limits.first &&
                    mat_ids[e->getID()] <= mat_limits.second &&
                    (*node)[2] >= elevation_limits.first &&
                    (*node)[2] <= elevation_limits.second)
                {
                    pnt_nodes.push_back(node->getID());
                    break;
                }
            }
        }
    }
    if (pnt_nodes.size() < 2)
    {
        WARN("No nodes found for point {:d}...", pnt.getID());
    }
    else
    {
        // sort found nodes from top to bottom (required for BHE simulations)
        std::sort(pnt_nodes.begin(), pnt_nodes.end(),
                  [nodes](std::size_t a, std::size_t b) {
                      return (*nodes[a])[2] > (*nodes[b])[2];
                  });
    }
    return pnt_nodes;
}

int main(int argc, char* argv[])
{
    TCLAP::CmdLine cmd(
        "Integrates line elements representing boreholes into a pre-existing "
        "3D mesh. Corresponding nodes matching the (x,y)-coordinates given in "
        "the gml-file are found in the mesh and connected from top to bottom "
        "via line elements. Each borehole (i.e. all points at a given "
        "(x,y)-location but at different depths) is assigned a unique material "
        "ID. Vertical limits of boreholes can be specified via Material IDs "
        "and/or elevation. Point not matchin any mesh nodes or located outside "
        "the mesh are ignored.\n\n"
        "OpenGeoSys-6 software, version " +
            GitInfoLib::GitInfo::ogs_version +
            ".\n"
            "Copyright (c) 2012-2021, OpenGeoSys Community "
            "(http://www.opengeosys.org)",
        ' ', GitInfoLib::GitInfo::ogs_version);

    double const dmax = std::numeric_limits<double>::max();
    TCLAP::ValueArg<double> max_elevation_arg(
        "", "max-elevation", "Maximum elevation for an integrated borehole",
        false, 0, "a number");
    cmd.add(max_elevation_arg);
    TCLAP::ValueArg<double> min_elevation_arg(
        "", "min-elevation", "Minimum elevation for an integrated borehole",
        false, 0, "a number");
    cmd.add(min_elevation_arg);
    TCLAP::ValueArg<int> max_id_arg(
        "", "max-id", "Maximum MaterialID for an integrated borehole", false,
        -1, "a number");
    cmd.add(max_id_arg);
    TCLAP::ValueArg<int> min_id_arg(
        "", "min-id", "Minimum MaterialID for an integrated borehole", false,
        -1, "a number");
    cmd.add(min_id_arg);
    TCLAP::ValueArg<std::string> geo_arg("g", "geo",
                                         "Name of the geometry file (*.gml)",
                                         true, "", "geometry file name");
    cmd.add(geo_arg);
    TCLAP::ValueArg<std::string> output_arg("o", "output",
                                            "Name of the output mesh (*.vtu)",
                                            true, "", "output file name");
    cmd.add(output_arg);
    TCLAP::ValueArg<std::string> input_arg("i", "input",
                                           "Name of the input mesh (*.vtu)",
                                           true, "", "input file name");
    cmd.add(input_arg);
    cmd.parse(argc, argv);

    std::pair<int, int> mat_limits(0, std::numeric_limits<int>::max());
    std::pair<double, double> elevation_limits(
        std::numeric_limits<double>::lowest(), dmax);

    if (min_id_arg.isSet() != max_id_arg.isSet())
    {
        ERR("If minimum MaterialID is set, maximum ID must be set, too (and "
            "vice versa).");
        return EXIT_FAILURE;
    }
    if (min_id_arg.isSet() && max_id_arg.isSet())
    {
        mat_limits =
            std::make_pair(min_id_arg.getValue(), max_id_arg.getValue());
    }
    if (mat_limits.first > mat_limits.second)
    {
        std::swap(mat_limits.first, mat_limits.second);
    }
    if (min_id_arg.isSet() && (mat_limits.first < 0 || mat_limits.second < 0))
    {
        ERR("Specified MaterialIDs must have non-negative values.");
        return EXIT_FAILURE;
    }
    if (min_elevation_arg.isSet() != max_elevation_arg.isSet())
    {
        ERR("If minimum elevation is set, maximum elevation must be set, too "
            "(and vice versa).");
        return EXIT_FAILURE;
    }
    if (min_elevation_arg.isSet() && max_elevation_arg.isSet())
    {
        elevation_limits = std::make_pair(min_elevation_arg.getValue(),
                                          max_elevation_arg.getValue());
    }
    if (elevation_limits.first > elevation_limits.second)
    {
        std::swap(elevation_limits.first, elevation_limits.second);
    }

    std::string const& mesh_name = input_arg.getValue();
    std::string const& output_name = output_arg.getValue();
    std::string const& geo_name = geo_arg.getValue();

    GeoLib::GEOObjects geo;
    GeoLib::IO::BoostXmlGmlInterface xml_io(geo);
    if (!xml_io.readFile(geo_name))
    {
        ERR("Failed to read geometry file `{:s}'.", geo_name);
        return EXIT_FAILURE;
    }
    std::vector<GeoLib::Point*> const& points =
        *geo.getPointVec(geo.getGeometryNames()[0]);

    std::unique_ptr<MeshLib::Mesh> const mesh(
        MeshLib::IO::readMeshFromFile(mesh_name));
    if (mesh == nullptr)
    {
        ERR("Failed to read input mesh file `{:s}'.", mesh_name);
        return EXIT_FAILURE;
    }
    if (mesh->getDimension() != 3)
    {
        ERR("Method can only be applied to 3D meshes.");
        return EXIT_FAILURE;
    }

    auto const& nodes = mesh->getNodes();
    auto const& mat_ids = MeshLib::materialIDs(*mesh);
    if (mat_ids == nullptr)
    {
        ERR("Mesh is required to have MaterialIDs");
        return EXIT_FAILURE;
    }

    auto const& elems = mesh->getElements();
    MeshLib::Properties props;
    auto new_mat_ids = props.createNewPropertyVector<int>(
        "MaterialIDs", MeshLib::MeshItemType::Cell);
    std::copy(mat_ids->cbegin(), mat_ids->cend(),
              std::back_inserter(*new_mat_ids));
    int const max_id = *std::max_element(mat_ids->begin(), mat_ids->end());
    std::vector<MeshLib::Node*> new_nodes = MeshLib::copyNodeVector(nodes);
    std::size_t const n_points = points.size();
    std::vector<MeshLib::Element*> new_elems =
        MeshLib::copyElementVector(elems, new_nodes);

    for (std::size_t i = 0; i < n_points; ++i)
    {
        std::vector<std::size_t> const& line_nodes =
            getNodes(*points[i], nodes, *mat_ids, mat_limits, elevation_limits);
        std::size_t const n_line_nodes = line_nodes.size();
        if (n_line_nodes < 2)
        {
            continue;
        }
        for (std::size_t j = 0; j < n_line_nodes - 1; ++j)
        {
            new_elems.push_back(new MeshLib::Line(
                {new_nodes[line_nodes[j]], new_nodes[line_nodes[j + 1]]},
                elems.size()));
            new_mat_ids->push_back(max_id + i + 1);
        }
    }

    MeshLib::Mesh const result("result", new_nodes, new_elems, props);
    MeshLib::IO::VtuInterface vtu(&result);
    vtu.writeToFile(output_name);
    return EXIT_SUCCESS;
}
