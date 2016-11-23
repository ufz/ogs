/**
 * @copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include <map>
#include <memory>
#include <vector>

#include <tclap/CmdLine.h>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include "Applications/ApplicationsLib/LogogSetup.h"

#include "BaseLib/FileTools.h"

#include "MeshLib/IO/readMeshFromFile.h"
#include "MeshLib/IO/writeMeshToFile.h"

#include "MeshLib/Mesh.h"

#include "ProcessLib/LIE/Common/MeshUtils.h"
#include "ProcessLib/LIE/Common/PostUtils.h"


int main (int argc, char* argv[])
{
    ApplicationsLib::LogogSetup logog_setup;

    TCLAP::CmdLine cmd("Postp-process results of the LIE approach",
                       ' ', "0.1");
    TCLAP::ValueArg<std::string> arg_out_pvd("o", "output-file",
                                          "the name of the new PVD file", true,
                                          "", "path");
    cmd.add(arg_out_pvd);
    TCLAP::ValueArg<std::string> arg_in_pvd("i", "input-file",
                                         "the original PVD file name", true,
                                         "", "path");
    cmd.add(arg_in_pvd);

    cmd.parse(argc, argv);

    auto const in_pvd_filename = arg_in_pvd.getValue();
    auto const in_pvd_file_dir = BaseLib::extractPath(in_pvd_filename);
    auto const out_pvd_filename = arg_out_pvd.getValue();
    auto const out_pvd_file_dir = BaseLib::extractPath(out_pvd_filename);
    INFO("start reading the PVD file %s", in_pvd_filename.c_str());
    boost::property_tree::ptree pt;
    read_xml(arg_in_pvd.getValue(), pt,  boost::property_tree::xml_parser::trim_whitespace);

    for (auto& dataset : pt.get_child("VTKFile.Collection"))
    {
        if (dataset.first != "DataSet")
            continue;

        // read VTU with simulation results
        auto const org_vtu_filename = dataset.second.get<std::string>("<xmlattr>.file");
        auto const org_vtu_filebasename = BaseLib::extractBaseName(org_vtu_filename);
        auto org_vtu_dir = BaseLib::extractPath(org_vtu_filename);
        if (org_vtu_dir.empty())
            org_vtu_dir = in_pvd_file_dir;
        auto const org_vtu_filepath = BaseLib::joinPaths(org_vtu_dir, org_vtu_filebasename);
        INFO("processing %s...", org_vtu_filepath.c_str());

        std::unique_ptr<MeshLib::Mesh const> mesh(
            MeshLib::IO::readMeshFromFile(org_vtu_filepath));

        // post-process
        std::vector<MeshLib::Element*> vec_matrix_elements;
        std::vector<int> vec_fracture_mat_IDs;
        std::vector<std::vector<MeshLib::Element*>> vec_fracture_elements;
        std::vector<std::vector<MeshLib::Element*>> vec_fracture_matrix_elements;
        std::vector<std::vector<MeshLib::Node*>> vec_fracture_nodes;
        ProcessLib::LIE::getFractureMatrixDataInMesh(
            *mesh, vec_matrix_elements, vec_fracture_mat_IDs, vec_fracture_elements,
            vec_fracture_matrix_elements, vec_fracture_nodes);

        ProcessLib::LIE::PostProcessTool post(
            *mesh, vec_fracture_nodes, vec_fracture_matrix_elements);

        // create a new VTU file and update XML
        auto const dest_vtu_filename = "post_" + org_vtu_filebasename;
        auto const dest_vtu_filepath = BaseLib::joinPaths(out_pvd_file_dir, dest_vtu_filename);
        INFO("create %s", dest_vtu_filepath.c_str());
        MeshLib::IO::writeMeshToFile(post.getOutputMesh(), dest_vtu_filepath);

        dataset.second.put("<xmlattr>.file", dest_vtu_filename);
    }

    // save into the new PVD file
    INFO("save into the new PVD file %s", out_pvd_filename.c_str());
    boost::property_tree::xml_writer_settings<std::string> settings('\t', 1);
    write_xml(arg_out_pvd.getValue(), pt, std::locale(), settings);

    return EXIT_SUCCESS;
}
