/**
 * @brief Computes the areas associated nodes of the surface mesh.
 *
 * @copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include <fstream>
#include <memory>
#include <numeric>
#include <string>
#include <vector>

#include <tclap/CmdLine.h>

#include "Applications/ApplicationsLib/LogogSetup.h"

#include "BaseLib/Error.h"
#include "BaseLib/FileTools.h"

#include "MeshLib/IO/readMeshFromFile.h"

#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshLib/MeshSurfaceExtraction.h"

static
void writeToFile(std::string const& id_area_fname, std::string const& csv_fname,
    std::vector<std::pair<std::size_t, double>> const& ids_and_areas,
    std::vector<MeshLib::Node*> const& mesh_nodes)
{
    std::ofstream ids_and_area_out(id_area_fname);
    if (!ids_and_area_out) {
        OGS_FATAL("Unable to open the file \"%s\" - aborting.", id_area_fname.c_str());
    }
    std::ofstream csv_out(csv_fname);
    if (!csv_out) {
        OGS_FATAL("Unable to open the file \"%s\" - aborting.", csv_fname.c_str());
    }

    ids_and_area_out.precision(std::numeric_limits<double>::digits10);
    csv_out.precision(std::numeric_limits<double>::digits10);

    csv_out << "ID x y z area node_id\n"; // CSV header
    for (std::size_t k(0); k<ids_and_areas.size(); k++) {
        ids_and_area_out << ids_and_areas[k].first << " "
                         << ids_and_areas[k].second << "\n";
        csv_out << k << " " << *(mesh_nodes[k]) << ids_and_areas[k].second
                << " " << ids_and_areas[k].first << "\n";
    }
    ids_and_area_out << "\n";
    csv_out << "\n";
}

int main (int argc, char* argv[])
{
    ApplicationsLib::LogogSetup logog_setup;

    TCLAP::CmdLine cmd("The tool computes the area per node of the surface mesh"
        " and writes the information as txt and csv data.", ' ', "0.1");
    TCLAP::ValueArg<std::string> mesh_in("i", "mesh-input-file",
        "the name of the file containing the input mesh", true,
        "", "file name of input mesh");
    cmd.add(mesh_in);
    TCLAP::ValueArg<std::string> id_prop_name("", "id-prop-name",
        "the name of the property containing the id information", false,
        "OriginalSubsurfaceNodeIDs", "property name");
    cmd.add(id_prop_name);
    TCLAP::ValueArg<std::string> out_base_fname("p", "output-base-name",
        "the path and base file name the output will be written to", false,
        "", "output path and base name as one string");
    cmd.add(out_base_fname);

    cmd.parse(argc, argv);

    std::unique_ptr<MeshLib::Mesh> surface_mesh(
        MeshLib::IO::readMeshFromFile(mesh_in.getValue()));
    INFO("Mesh read: %u nodes, %u elements.", surface_mesh->getNumberOfNodes(),
         surface_mesh->getNumberOfElements());
    // ToDo check if mesh is read correct and if the mesh is a surface mesh

    // check if a node property containing the subsurface ids is available
    auto* orig_node_ids =
        surface_mesh->getProperties().getPropertyVector<std::size_t>(
            id_prop_name.getValue());
    // if the node property is not available generate it
    if (!orig_node_ids)
    {
        orig_node_ids =
            surface_mesh->getProperties().createNewPropertyVector<std::size_t>(
                id_prop_name.getValue(), MeshLib::MeshItemType::Node, 1);
        if (!orig_node_ids)
        {
            ERR("Fatal error: could not create property.");
            return EXIT_FAILURE;
        }
        orig_node_ids->resize(surface_mesh->getNumberOfNodes());
        std::iota(orig_node_ids->begin(), orig_node_ids->end(), 0);
    }

    std::vector<double> areas(
        MeshLib::MeshSurfaceExtraction::getSurfaceAreaForNodes(*surface_mesh));

    // pack area and node id together
    std::vector<std::pair<std::size_t, double>> ids_and_areas;
    std::transform(orig_node_ids->cbegin(), orig_node_ids->cend(),
                   areas.cbegin(), std::back_inserter(ids_and_areas),
                   std::make_pair<std::size_t const&, double const&>);

    // generate file names for output
    std::string path(out_base_fname.getValue());
    if (path.empty())
        path = BaseLib::dropFileExtension(mesh_in.getValue());
    std::string const id_and_area_fname(path+".txt");
    std::string const csv_fname(path+".csv");

    writeToFile(id_and_area_fname, csv_fname, ids_and_areas,
                surface_mesh->getNodes());

    return EXIT_SUCCESS;
}
