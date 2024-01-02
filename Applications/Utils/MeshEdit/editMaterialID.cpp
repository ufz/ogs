/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <tclap/CmdLine.h>

#ifdef USE_PETSC
#include <mpi.h>
#endif

#include <memory>

#include "InfoLib/GitInfo.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/IO/readMeshFromFile.h"
#include "MeshLib/IO/writeMeshToFile.h"
#include "MeshLib/Mesh.h"
#include "MeshToolsLib/MeshEditing/ElementValueModification.h"
#include "MeshToolsLib/MeshInformation.h"

int main(int argc, char* argv[])
{
    TCLAP::CmdLine cmd(
        "Edit material IDs of mesh elements.\n\n"
        "OpenGeoSys-6 software, version " +
            GitInfoLib::GitInfo::ogs_version +
            ".\n"
            "Copyright (c) 2012-2024, OpenGeoSys Community "
            "(http://www.opengeosys.org)",
        ' ', GitInfoLib::GitInfo::ogs_version);
    TCLAP::SwitchArg replaceArg("r", "replace", "replace material IDs", false);
    TCLAP::SwitchArg condenseArg("c", "condense", "condense material IDs",
                                 false);
    TCLAP::SwitchArg specifyArg(
        "s", "specify", "specify material IDs by element types (-e)", false);
    std::vector<TCLAP::Arg*> vec_xors;
    vec_xors.push_back(&replaceArg);
    vec_xors.push_back(&condenseArg);
    vec_xors.push_back(&specifyArg);
    cmd.xorAdd(vec_xors);
    TCLAP::ValueArg<std::string> mesh_in(
        "i", "mesh-input-file",
        "the name of the file containing the input mesh", true, "",
        "file name");
    cmd.add(mesh_in);
    TCLAP::ValueArg<std::string> mesh_out(
        "o", "mesh-output-file",
        "the name of the file the mesh will be written to", true, "",
        "file name");
    cmd.add(mesh_out);
    TCLAP::MultiArg<unsigned> matIDArg("m", "current-material-id",
                                       "current material id to be replaced",
                                       false, "number");
    cmd.add(matIDArg);
    TCLAP::ValueArg<unsigned> newIDArg("n", "new-material-id",
                                       "new material id", false, 0, "number");
    cmd.add(newIDArg);
    std::vector<std::string> eleList(MeshLib::getMeshElemTypeStringsShort());
    TCLAP::ValuesConstraint<std::string> allowedVals(eleList);
    TCLAP::ValueArg<std::string> eleTypeArg("e", "element-type", "element type",
                                            false, "", &allowedVals);
    cmd.add(eleTypeArg);

    cmd.parse(argc, argv);

#ifdef USE_PETSC
    MPI_Init(&argc, &argv);
#endif

    if (!replaceArg.isSet() && !condenseArg.isSet() && !specifyArg.isSet())
    {
        INFO("Please select editing mode: -r or -c or -s");
#ifdef USE_PETSC
        MPI_Finalize();
#endif
        return 0;
    }
    if (replaceArg.isSet() && condenseArg.isSet())
    {
        INFO("Please select only one editing mode: -r or -c or -s");
#ifdef USE_PETSC
        MPI_Finalize();
#endif
        return 0;
    }
    if (replaceArg.isSet())
    {
        if (!matIDArg.isSet() || !newIDArg.isSet())
        {
            INFO(
                "current and new material IDs must be provided for "
                "replacement");
#ifdef USE_PETSC
            MPI_Finalize();
#endif
            return 0;
        }
    }
    else if (specifyArg.isSet())
    {
        if (!eleTypeArg.isSet() || !newIDArg.isSet())
        {
            INFO(
                "element type and new material IDs must be provided to specify "
                "elements");
#ifdef USE_PETSC
            MPI_Finalize();
#endif
            return 0;
        }
    }

    std::unique_ptr<MeshLib::Mesh> mesh(
        MeshLib::IO::readMeshFromFile(mesh_in.getValue()));
    INFO("Mesh read: {:d} nodes, {:d} elements.", mesh->getNumberOfNodes(),
         mesh->getNumberOfElements());

    if (condenseArg.isSet())
    {
        INFO("Condensing material ID...");
        INFO("The MaterialIDs of the input file: [{}]",
             fmt::join(MeshToolsLib::MeshInformation::getMaterialIDs(*mesh),
                       ", "));
        MeshToolsLib::ElementValueModification::condense(*mesh);
        INFO("The MaterialIDs of the output file: [{}]",
             fmt::join(MeshToolsLib::MeshInformation::getMaterialIDs(*mesh),
                       ", "));
    }
    else if (replaceArg.isSet())
    {
        INFO("Replacing material ID...");
        INFO("The MaterialIDs of the input file: [{}]",
             fmt::join(MeshToolsLib::MeshInformation::getMaterialIDs(*mesh),
                       ", "));

        const auto vecOldID = matIDArg.getValue();
        const unsigned newID = newIDArg.getValue();
        for (auto oldID : vecOldID)
        {
            INFO("{:d} -> {:d}", oldID, newID);
            MeshToolsLib::ElementValueModification::replace(*mesh, oldID, newID,
                                                            true);
        }
        INFO("The MaterialIDs of the output file: [{}]",
             fmt::join(MeshToolsLib::MeshInformation::getMaterialIDs(*mesh),
                       ", "));
    }
    else if (specifyArg.isSet())
    {
        INFO("Specifying material ID...");
        INFO("The MaterialIDs of the input file: [{}]",
             fmt::join(MeshToolsLib::MeshInformation::getMaterialIDs(*mesh),
                       ", "));

        const std::string eleTypeName(eleTypeArg.getValue());
        const MeshLib::MeshElemType eleType =
            MeshLib::String2MeshElemType(eleTypeName);
        const unsigned newID = newIDArg.getValue();
        unsigned cnt = MeshToolsLib::ElementValueModification::setByElementType(
            *mesh, eleType, newID);
        INFO("updated {:d} elements", cnt);
        INFO("The MaterialIDs of the output file: [{}]",
             fmt::join(MeshToolsLib::MeshInformation::getMaterialIDs(*mesh),
                       ", "));
    }
    MeshLib::IO::writeMeshToFile(*mesh, mesh_out.getValue());

#ifdef USE_PETSC
    MPI_Finalize();
#endif
    return EXIT_SUCCESS;
}
