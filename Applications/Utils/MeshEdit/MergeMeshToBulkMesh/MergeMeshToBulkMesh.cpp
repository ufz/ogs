/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on March 26, 2025, 12:58 PM
 */

#include "MeshToolsLib/MeshEditing/MergeMeshToBulkMesh.h"

#include <tclap/CmdLine.h>

#include <algorithm>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "BaseLib/FileTools.h"
#include "BaseLib/MPI.h"
#include "BaseLib/RunTime.h"
#include "BaseLib/StringTools.h"
#include "InfoLib/GitInfo.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/IO/VtkIO/VtuInterface.h"
#include "MeshLib/IO/readMeshFromFile.h"
#include "MeshLib/IO/writeMeshToFile.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"

int main(int argc, char* argv[])
{
    TCLAP::CmdLine cmd(
        "Tool merges one mesh to a bulk mesh.\n\n"
        "OpenGeoSys-6 software, version " +
            GitInfoLib::GitInfo::ogs_version +
            ".\n"
            "Copyright (c) 2012-2025, OpenGeoSys Community "
            "(http://www.opengeosys.org)",
        ' ', GitInfoLib::GitInfo::ogs_version);
    TCLAP::ValueArg<std::string> bulk_mesh_in(
        "b", "bulk-mesh-input-file",
        "Input (.vtk | .msh). The name of the file containing the input bulk "
        "mesh",
        true, "", "INPUT_FILE");
    cmd.add(bulk_mesh_in);
    TCLAP::ValueArg<std::string> mesh_in(
        "i", "mesh-input-file",
        "Input (.vtk | .msh). The name of the file containing the input mesh "
        "to be merged",
        true, "", "INPUT_FILE");
    cmd.add(mesh_in);
    TCLAP::ValueArg<std::string> mesh_out(
        "o", "mesh-output-file",
        "Output (.vtk | .msh). The name of the file the merged mesh should be "
        "written to",
        true, "", "OUTPUT_FILE");
    cmd.add(mesh_out);

    TCLAP::ValueArg<double> p("", "pressure",
                              "initial pressure value in the mesh to be merged",
                              false, 0.0, "PRESSURE");
    cmd.add(p);

    TCLAP::ValueArg<double> pg(
        "", "gas_pressure",
        "initial gas pressure value in the mesh to be merged", false, 0.0,
        "GAS_PRESSURE");
    cmd.add(pg);

    TCLAP::ValueArg<double> pc(
        "", "capillary_pressure",
        "initial capillary pressure value in the mesh to be merged", false, 0.0,
        "CAPILLARY_PRESSURE");
    cmd.add(pc);

    TCLAP::ValueArg<double> T(
        "", "temperature", "initial temperature value in the mesh to be merged",
        false, 290.0, "TEMPERATURE");
    cmd.add(T);

    TCLAP::ValueArg<double> sxx(
        "", "sigma_xx", "initial stress xx value in the mesh to be merged",
        false, 0.0, "SIGMA_XX");
    cmd.add(sxx);
    TCLAP::ValueArg<double> syy(
        "", "sigma_yy", "initial stress yy value in the mesh to be merged",
        false, 0.0, "SIGMA_YY");
    cmd.add(syy);
    TCLAP::ValueArg<double> szz(
        "", "sigma_zz", "initial stress zz value in the mesh to be merged",
        false, 0.0, "SIGMA_ZZ");
    cmd.add(szz);

    TCLAP::ValueArg<int> mat_id(
        "", "material_id", "Material ID of the mesh to be merged, (min = 0)",
        false, 0.0, "MATERIAL_ID");
    cmd.add(mat_id);

    cmd.parse(argc, argv);

    BaseLib::MPI::Setup mpi_setup(argc, argv);

    BaseLib::RunTime timer;
    timer.start();

    std::unordered_map<std::string, double> initial_value_dict;
    initial_value_dict.insert({"p", p.getValue()});
    initial_value_dict.insert({"pg", pg.getValue()});
    initial_value_dict.insert({"pc", pc.getValue()});
    initial_value_dict.insert({"T", T.getValue()});
    initial_value_dict.insert({"sxx", sxx.getValue()});
    initial_value_dict.insert({"syy", syy.getValue()});
    initial_value_dict.insert({"szz", szz.getValue()});
    initial_value_dict.insert({"mat_id", mat_id.getValue()});

    auto read_mesh = [](std::string const& mesh_file_name)
    {
        std::unique_ptr<MeshLib::Mesh> mesh(MeshLib::IO::readMeshFromFile(
            mesh_file_name, true /* compute_element_neighbors */));

        if (!mesh)
        {
            OGS_FATAL("Could not read the mesh {:s}", mesh_file_name);
        }

        INFO("Read {:s}: {:d} nodes, {:d} elements.", mesh_file_name,
             mesh->getNumberOfNodes(), mesh->getNumberOfElements());
        return mesh;
    };

    auto bulk_mesh = read_mesh(bulk_mesh_in.getValue());
    auto const mesh_to_be_merged = read_mesh(mesh_in.getValue());

    auto merged_mesh = MeshToolsLib::mergeMeshToBulkMesh(
        *bulk_mesh, *mesh_to_be_merged, initial_value_dict);

    MeshLib::IO::VtuInterface writer(merged_mesh.get());

    auto const result = writer.writeToFile(mesh_out.getValue());
    if (!result)
    {
        ERR("Could not write mesh to '{:s}'.", mesh_out.getValue());
        return EXIT_FAILURE;
    }
    INFO("It took {} s", timer.elapsed());

    // The merged nodes are moved to the new mesh, and the duplicated nodes on
    // the mesh interface have already been deleted in
    // MeshToolsLib::mergeMeshToBulkMesh. The elements are also moved to the
    // merged mesh. The merged mesh object handles the release of the memory
    // allocated for nodes and elements. Therefore, shallowClean is called.
    mesh_to_be_merged->shallowClean();
    bulk_mesh->shallowClean();

    return EXIT_SUCCESS;
}
