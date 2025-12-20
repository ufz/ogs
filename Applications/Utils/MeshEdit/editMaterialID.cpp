// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <spdlog/fmt/ranges.h>
#include <tclap/CmdLine.h>

#include <memory>

#include "BaseLib/Logging.h"
#include "BaseLib/MPI.h"
#include "BaseLib/TCLAPArguments.h"
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
            "Copyright (c) 2012-2025, OpenGeoSys Community "
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
        "Input (.vtu | .msh). The name of the file containing the input mesh",
        true, "", "INPUT_FILE");
    cmd.add(mesh_in);
    TCLAP::ValueArg<std::string> mesh_out(
        "o", "mesh-output-file",
        "Output (.vtu | .msh). The name of the file the mesh will be written "
        "to",
        true, "", "OUTPUT_FILE");
    cmd.add(mesh_out);
    TCLAP::MultiArg<unsigned> matIDArg("m", "current-material-id",
                                       "current material id to be replaced",
                                       false, "CURRENT_MATERIAL_ID");
    cmd.add(matIDArg);
    TCLAP::ValueArg<unsigned> newIDArg(
        "n", "new-material-id", "new material id", false, 0, "NEW_MATERIAL_ID");
    cmd.add(newIDArg);

    // TODO: FIND A BETTER SOLUTION FOR ALLOWED ELEM TYPES DEFINITION
    std::vector<std::string> eleList(MeshLib::getMeshElemTypeStringsShort());
    TCLAP::ValuesConstraint<std::string> allowedVals(eleList);
    std::vector<std::string> allowed_elems_vector{
        "point", "line", "quad", "hex", "tri", "tet", "pris", "pyra"};
    TCLAP::ValuesConstraint<std::string> allowed_elems(allowed_elems_vector);
    TCLAP::ValueArg<std::string> eleTypeArg("e", "element-type", "element type",
                                            false, "", &allowed_elems);
    cmd.add(eleTypeArg);
    auto log_level_arg = BaseLib::makeLogLevelArg();
    cmd.add(log_level_arg);

    cmd.parse(argc, argv);

    BaseLib::MPI::Setup mpi_setup(argc, argv);
    BaseLib::initOGSLogger(log_level_arg.getValue());

    if (!replaceArg.isSet() && !condenseArg.isSet() && !specifyArg.isSet())
    {
        INFO("Please select editing mode: -r or -c or -s");
        return 0;
    }
    if (replaceArg.isSet() && condenseArg.isSet())
    {
        INFO("Please select only one editing mode: -r or -c or -s");
        return 0;
    }
    if (replaceArg.isSet())
    {
        if (!matIDArg.isSet() || !newIDArg.isSet())
        {
            INFO(
                "current and new material IDs must be provided for "
                "replacement");
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

    return EXIT_SUCCESS;
}
