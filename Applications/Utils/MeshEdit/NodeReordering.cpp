// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <tclap/CmdLine.h>

#include <vector>

#include "BaseLib/Algorithm.h"
#include "BaseLib/Logging.h"
#include "BaseLib/MPI.h"
#include "BaseLib/TCLAPArguments.h"
#include "InfoLib/GitInfo.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/IO/readMeshFromFile.h"
#include "MeshLib/IO/writeMeshToFile.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshEnums.h"
#include "MeshLib/Node.h"
#include "MeshLib/Utils/GetSpaceDimension.h"
#include "MeshToolsLib/MeshEditing/ReverseNodeOrder.h"

/// Fixes inconsistencies between VTK's and OGS' node order for prism
/// elements.
/// In particular, this fixes issues between OGS6 meshes with and without
/// InSitu-Lib
void fixVtkInconsistencies(std::vector<MeshLib::Element*>& elements)
{
    for (auto* const element : elements)
    {
        const unsigned nElemNodes(element->getNumberOfBaseNodes());
        std::vector<MeshLib::Node*> nodes(element->getNodes(),
                                          element->getNodes() + nElemNodes);

        for (std::size_t j = 0; j < nElemNodes; ++j)
        {
            if (element->getGeomType() == MeshLib::MeshElemType::PRISM)
            {
                for (std::size_t k = 0; k < 3; ++k)
                {
                    element->setNode(k, nodes[k + 3]);
                    element->setNode(k + 3, nodes[k]);
                }
                break;
            }
        }
    }
}

/// Orders the base nodes of each elements before its non-linear nodes.
void reorderNonlinearNodes(MeshLib::Mesh& mesh)
{
    std::vector<MeshLib::Node*> base_nodes;
    std::vector<MeshLib::Node*> nonlinear_nodes;
    for (MeshLib::Element const* e : mesh.getElements())
    {
        for (unsigned i = 0; i < e->getNumberOfBaseNodes(); i++)
        {
            base_nodes.push_back(const_cast<MeshLib::Node*>(e->getNode(i)));
        }
        for (unsigned i = e->getNumberOfBaseNodes(); i < e->getNumberOfNodes();
             i++)
        {
            nonlinear_nodes.push_back(
                const_cast<MeshLib::Node*>(e->getNode(i)));
        }
    }

    BaseLib::makeVectorUnique(base_nodes,
                              MeshLib::idsComparator<MeshLib::Node*>);
    BaseLib::makeVectorUnique(nonlinear_nodes,
                              MeshLib::idsComparator<MeshLib::Node*>);

    std::vector<MeshLib::Node*>& allnodes =
        const_cast<std::vector<MeshLib::Node*>&>(mesh.getNodes());
    allnodes.clear();

    allnodes.insert(allnodes.end(), base_nodes.begin(), base_nodes.end());
    allnodes.insert(allnodes.end(), nonlinear_nodes.begin(),
                    nonlinear_nodes.end());

    mesh.resetNodeIDs();
}

int main(int argc, char* argv[])
{
    TCLAP::CmdLine cmd(
        "Reorders mesh nodes in elements to make old or incorrectly ordered "
        "meshes compatible with OGS6.\n"
        "Four options are available:\n"
        "Method 0: As Method 1, but additionally verifying the Jacobian "
        "determinant of the reordered elements.\n"
        "Method 1: Reversing order of nodes unless it's perceived correct by "
        "OGS6 standards. This is the default selection.\n"
        "Method 2: Fixing node ordering issues between VTK and OGS6 (only "
        "applies to prism-elements)\n"
        "Method 3: Re-ordering of mesh node vector such that all base nodes "
        "are sorted before all nonlinear nodes.\n\n"
        "OpenGeoSys-6 software, version " +
            GitInfoLib::GitInfo::ogs_version +
            ".\n"
            "Copyright (c) 2012-2026, OpenGeoSys Community "
            "(http://www.opengeosys.org)",
        ' ', GitInfoLib::GitInfo::ogs_version);

    TCLAP::SwitchArg no_volume_check_arg(
        "", "no_volume_check",
        "By default the volumes of original and reordered elements are "
        "compared if they are numerically equal, i.e., if the absolute "
        "difference or the relative volume difference is at most 100 machine "
        "epsilons. This switch disables the volume comparison.");
    cmd.add(no_volume_check_arg);

    std::vector<int> method_ids{0, 1, 2, 3};
    TCLAP::ValuesConstraint<int> allowed_values(method_ids);
    TCLAP::ValueArg<int> method_arg("m", "method",
                                    "reordering method selection", false, 1,
                                    &allowed_values);
    cmd.add(method_arg);
    TCLAP::ValueArg<std::string> output_mesh_arg(
        "o", "output_mesh", "Output (.vtu). The name of the output mesh file",
        true, "", "OUTPUT_FILE");
    cmd.add(output_mesh_arg);
    TCLAP::ValueArg<std::string> input_mesh_arg(
        "i", "input_mesh",
        "Input (.vtu | .vtk | .msh). The name of the input mesh file", true, "",
        "INPUT_FILE");
    cmd.add(input_mesh_arg);
    auto log_level_arg = BaseLib::makeLogLevelArg();
    cmd.add(log_level_arg);
    cmd.parse(argc, argv);

    BaseLib::MPI::Setup mpi_setup(argc, argv);
    BaseLib::initOGSLogger(log_level_arg.getValue());

    std::unique_ptr<MeshLib::Mesh> mesh(
        MeshLib::IO::readMeshFromFile(input_mesh_arg.getValue()));

    if (!mesh)
    {
        return EXIT_FAILURE;
    }

    INFO("Reordering nodes... ");
    if (!method_arg.isSet() || method_arg.getValue() < 2)
    {
        bool const forced = (method_arg.getValue() == 0);

        if (forced)
        {
            INFO(
                "Method 0: The Jacobian determinant of the reordered elements "
                "is checked again.");
        }
        INFO(
            "Method: Reversing order of nodes unless it is considered correct "
            "by the OGS6 standard, i.e. such that det(J) > 0, where J is the "
            "Jacobian of the global-to-local coordinate transformation.");
        int const mesh_space_dimension =
            MeshLib::getSpaceDimension(mesh->getNodes());
        INFO("Corrected {:d} elements.",
             MeshToolsLib::reverseNodeOrder(mesh->getElements(),
                                            mesh_space_dimension, forced,
                                            !no_volume_check_arg.isSet()));
    }
    else if (method_arg.getValue() == 2)
    {
        fixVtkInconsistencies(
            const_cast<std::vector<MeshLib::Element*>&>(mesh->getElements()));
    }
    else if (method_arg.getValue() == 3)
    {
        reorderNonlinearNodes(*mesh);
    }

    if (MeshLib::IO::writeMeshToFile(*mesh, output_mesh_arg.getValue()) != 0)
    {
        return EXIT_FAILURE;
    }
    INFO("VTU file written.");

    return EXIT_SUCCESS;
}
