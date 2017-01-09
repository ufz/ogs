/**
 * @copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
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
#include "MeshLib/MeshEditing/ConvertToLinearMesh.h"

#include "ProcessLib/LIE/Common/MeshUtils.h"
#include "ProcessLib/LIE/Common/PostUtils.h"

namespace
{

void postVTU(std::string const& int_vtu_filename, std::string const& out_vtu_filename)
{
    // read VTU with simulation results
    std::unique_ptr<MeshLib::Mesh const> mesh(
        MeshLib::IO::readMeshFromFile(int_vtu_filename));
    if (mesh->isNonlinear())
        mesh = MeshLib::convertToLinearMesh(*mesh, mesh->getName());

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

    // create a new VTU file
    INFO("create %s", out_vtu_filename.c_str());
    MeshLib::IO::writeMeshToFile(post.getOutputMesh(), out_vtu_filename);
}


void postPVD(std::string const& in_pvd_filename, std::string const& out_pvd_filename)
{
    auto const in_pvd_file_dir = BaseLib::extractPath(in_pvd_filename);
    auto const out_pvd_file_dir = BaseLib::extractPath(out_pvd_filename);
    INFO("start reading the PVD file %s", in_pvd_filename.c_str());
    boost::property_tree::ptree pt;
    read_xml(in_pvd_filename, pt,  boost::property_tree::xml_parser::trim_whitespace);

    for (auto& dataset : pt.get_child("VTKFile.Collection"))
    {
        if (dataset.first != "DataSet")
            continue;

        // get VTU file name
        auto const org_vtu_filename = dataset.second.get<std::string>("<xmlattr>.file");
        auto const org_vtu_filebasename = BaseLib::extractBaseName(org_vtu_filename);
        auto org_vtu_dir = BaseLib::extractPath(org_vtu_filename);
        if (org_vtu_dir.empty())
            org_vtu_dir = in_pvd_file_dir;
        auto const org_vtu_filepath = BaseLib::joinPaths(org_vtu_dir, org_vtu_filebasename);
        INFO("processing %s...", org_vtu_filepath.c_str());

        // post-process the VTU and save into the new file
        auto const dest_vtu_filename = "post_" + org_vtu_filebasename;
        auto const dest_vtu_filepath = BaseLib::joinPaths(out_pvd_file_dir, dest_vtu_filename);
        postVTU(org_vtu_filepath, dest_vtu_filepath);

        // create a new VTU file and update XML
        dataset.second.put("<xmlattr>.file", dest_vtu_filename);
    }

    // save into the new PVD file
    INFO("save into the new PVD file %s", out_pvd_filename.c_str());
    boost::property_tree::xml_writer_settings<std::string> settings('\t', 1);
    write_xml(out_pvd_filename, pt, std::locale(), settings);
}

} // unnamed namespace


int main (int argc, char* argv[])
{
    ApplicationsLib::LogogSetup logog_setup;

    TCLAP::CmdLine cmd("Post-process results of the LIE approach",
                       ' ', "0.1");
    TCLAP::ValueArg<std::string> arg_out_file("o", "output-file",
                                          "the name of the new PVD or VTU file", true,
                                          "", "path");
    cmd.add(arg_out_file);
    TCLAP::ValueArg<std::string> arg_in_file("i", "input-file",
                                         "the original PVD or VTU file name", true,
                                         "", "path");
    cmd.add(arg_in_file);

    cmd.parse(argc, argv);

    auto const in_file_ext = BaseLib::getFileExtension(arg_in_file.getValue());
    if (in_file_ext == "pvd")
        postPVD(arg_in_file.getValue(), arg_out_file.getValue());
    else if (in_file_ext == "vtu")
        postVTU(arg_in_file.getValue(), arg_out_file.getValue());
    else
        OGS_FATAL("The given file type (%s) is not supported.", in_file_ext.c_str());

    return EXIT_SUCCESS;
}
