/**
 * \file
 * 2013/13/06 KR Initial implementation
 *
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <tclap/CmdLine.h>

#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <memory>
#include <vector>

#include "BaseLib/Algorithm.h"
#include "BaseLib/Logging.h"
#include "BaseLib/MPI.h"
#include "BaseLib/TCLAPArguments.h"
#include "InfoLib/GitInfo.h"
#include "MeshLib/ElementCoordinatesMappingLocal.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/IO/readMeshFromFile.h"
#include "MeshLib/IO/writeMeshToFile.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshLib/Utils/GetSpaceDimension.h"
#include "NumLib/Fem/ShapeFunction/ShapeHex20.h"
#include "NumLib/Fem/ShapeFunction/ShapeHex8.h"
#include "NumLib/Fem/ShapeFunction/ShapeLine2.h"
#include "NumLib/Fem/ShapeFunction/ShapeLine3.h"
#include "NumLib/Fem/ShapeFunction/ShapePoint1.h"
#include "NumLib/Fem/ShapeFunction/ShapePrism15.h"
#include "NumLib/Fem/ShapeFunction/ShapePrism6.h"
#include "NumLib/Fem/ShapeFunction/ShapePyra13.h"
#include "NumLib/Fem/ShapeFunction/ShapePyra5.h"
#include "NumLib/Fem/ShapeFunction/ShapeQuad4.h"
#include "NumLib/Fem/ShapeFunction/ShapeQuad8.h"
#include "NumLib/Fem/ShapeFunction/ShapeQuad9.h"
#include "NumLib/Fem/ShapeFunction/ShapeTet10.h"
#include "NumLib/Fem/ShapeFunction/ShapeTet4.h"
#include "NumLib/Fem/ShapeFunction/ShapeTri3.h"
#include "NumLib/Fem/ShapeFunction/ShapeTri6.h"

template <typename GradShapeFunction, int Dim>
bool checkJacobianDeterminant(MeshLib::Element const& e,
                              int const mesh_space_dimension,
                              std::array<double, Dim> const& xi,
                              bool const check_reordered = false)
{
    Eigen::Matrix<double, GradShapeFunction::DIM, GradShapeFunction::NPOINTS,
                  Eigen::RowMajor>
        dNdxi;
    Eigen::Map<Eigen::VectorXd> dN_vec(dNdxi.data(), dNdxi.size());
    GradShapeFunction::computeGradShapeFunction(xi.data(), dN_vec);
    MeshLib::ElementCoordinatesMappingLocal ele_local_coord(
        e, mesh_space_dimension);
    // check the determinant of the Jacobian
    Eigen::MatrixXd J = Eigen::MatrixXd::Zero(Dim, Dim);

    for (unsigned k = 0; k < GradShapeFunction::NPOINTS; k++)
    {
        const MathLib::Point3d& mapped_pt =
            ele_local_coord.getMappedCoordinates(k);
        // outer product of dNdr and mapped_pt for a particular node
        J += dNdxi.col(k) * mapped_pt.asEigenVector3d().head(Dim).transpose();
    }

    if (!check_reordered)
    {
        return !(J.determinant() < 0);
    }

    if (J.determinant() < 0)
    {
        OGS_FATAL(
            "Element {:d} has negative Jacobian determinant {:g}. "
            "NodeReordering fails.",
            e.getID(), J.determinant());
    }

    return true;
}

/**
 * \brief Reverses order of nodes. In particular, this fixes issues between
 * (Gmsh or OGS5) and OGS6 meshes.
 *
 * \param elements  Mesh elements whose nodes should be reordered
 * \param forced    If true, nodes are reordered for all
 * elements, if false it is first checked if the node order is correct
 * according to OGS6 element definitions.
 */
void reverseNodeOrder(std::vector<MeshLib::Element*>& elements,
                      int const mesh_space_dimension, bool const forced)
{
    std::size_t n_corrected_elements = 0;

    auto reverse_nodes_non_3d =
        [](MeshLib::Element& element, const std::vector<MeshLib::Node*>& nodes)
    {
        const unsigned nElemNodes = nodes.size();
        for (std::size_t j = 0; j < nElemNodes; ++j)
        {
            element.setNode(j, nodes[nElemNodes - j - 1]);
        }
    };

    for (auto const element : elements)
    {
        const unsigned nElemNodes(element->getNumberOfBaseNodes());
        std::vector<MeshLib::Node*> nodes(element->getNodes(),
                                          element->getNodes() + nElemNodes);

        switch (element->getGeomType())
        {
            case MeshLib::MeshElemType::TETRAHEDRON:
            {
                std::array<double, 3> xi = {0.25, 0.25, 0.25};
                if (bool const correct_node_ordering =
                        checkJacobianDeterminant<NumLib::ShapeTet4, 3>(
                            *element, mesh_space_dimension, xi);
                    correct_node_ordering)
                {
                    continue;
                };
                for (std::size_t j = 0; j < 4; ++j)
                {
                    element->setNode(j, nodes[(j + 1) % 4]);
                }
                if (forced)
                {
                    checkJacobianDeterminant<NumLib::ShapeTet4, 3>(
                        *element, mesh_space_dimension, xi,
                        true /*check_reordered*/);
                }

                n_corrected_elements++;
            }
            break;
            case MeshLib::MeshElemType::PRISM:
            {
                std::array<double, 3> xi = {1.0 / 3.0, 1.0 / 3.0, 0.5};
                if (bool const correct_node_ordering =
                        checkJacobianDeterminant<NumLib::ShapePrism6, 3>(
                            *element, mesh_space_dimension, xi);
                    correct_node_ordering)
                {
                    continue;
                };
                for (std::size_t j = 0; j < 3; ++j)
                {
                    element->setNode(j, nodes[j + 3]);
                    element->setNode(j + 3, nodes[j]);
                }
                if (forced)
                {
                    checkJacobianDeterminant<NumLib::ShapePrism6, 3>(
                        *element, mesh_space_dimension, xi,
                        true /*check_reordered*/);
                }

                n_corrected_elements++;
            }
            break;
            case MeshLib::MeshElemType::HEXAHEDRON:
            {
                std::array<double, 3> xi = {0.5, 0.5, 0.5};
                if (bool const correct_node_ordering =
                        checkJacobianDeterminant<NumLib::ShapeHex8, 3>(
                            *element, mesh_space_dimension, xi);
                    correct_node_ordering)
                {
                    continue;
                };
                for (std::size_t j = 0; j < 4; ++j)
                {
                    element->setNode(j, nodes[j + 4]);
                    element->setNode(j + 4, nodes[j]);
                }
                if (forced)
                {
                    checkJacobianDeterminant<NumLib::ShapeHex8, 3>(
                        *element, mesh_space_dimension, xi,
                        true /*check_reordered*/);
                }
                n_corrected_elements++;
            }
            break;
            case MeshLib::MeshElemType::LINE:
            {
                std::array<double, 1> xi = {0.5};
                if (bool const correct_node_ordering =
                        checkJacobianDeterminant<NumLib::ShapeLine2, 1>(
                            *element, mesh_space_dimension, xi);
                    correct_node_ordering)
                {
                    continue;
                };

                reverse_nodes_non_3d(*element, nodes);

                if (forced)
                {
                    checkJacobianDeterminant<NumLib::ShapeLine2, 1>(
                        *element, mesh_space_dimension, xi,
                        true /*check_reordered*/);
                }

                n_corrected_elements++;
            }
            break;
            case MeshLib::MeshElemType::QUAD:
            {
                std::array<double, 2> xi = {0.0, 0.0};
                if (bool const correct_node_ordering =
                        checkJacobianDeterminant<NumLib::ShapeQuad4, 2>(
                            *element, mesh_space_dimension, xi);
                    correct_node_ordering)
                {
                    continue;
                };
                reverse_nodes_non_3d(*element, nodes);

                if (forced)
                {
                    checkJacobianDeterminant<NumLib::ShapeQuad4, 2>(
                        *element, mesh_space_dimension, xi,
                        true /*check_reordered*/);
                }
                n_corrected_elements++;
            }
            break;
            case MeshLib::MeshElemType::TRIANGLE:
            {
                std::array<double, 2> xi = {1.0 / 3.0, 1.0 / 3.0};
                if (bool const correct_node_ordering =
                        checkJacobianDeterminant<NumLib::ShapeTri3, 2>(
                            *element, mesh_space_dimension, xi);
                    correct_node_ordering)
                {
                    continue;
                };
                reverse_nodes_non_3d(*element, nodes);
                if (forced)
                {
                    checkJacobianDeterminant<NumLib::ShapeTri3, 2>(
                        *element, mesh_space_dimension, xi,
                        true /*check_reordered*/);
                }
                n_corrected_elements++;
            }
            break;
            case MeshLib::MeshElemType::POINT:
                break;
            default:
                OGS_FATAL(
                    "Element type {:d} not supported for node "
                    "reordering.",
                    static_cast<int>(element->getGeomType()));
        }
    }

    INFO("Corrected {:d} elements.", n_corrected_elements);
}

/// Fixes inconsistencies between VTK's and OGS' node order for prism elements.
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
        "Three options are available:\n"
        "Method 0: Reversing order of nodes and checking again for all "
        "elements.\n"
        "Method 1: Reversing order of nodes unless it's perceived correct by "
        "OGS6 standards. This is the default selection.\n"
        "Method 2: Fixing node ordering issues between VTK and OGS6 (only "
        "applies to prism-elements)\n"
        "Method 3: Re-ordering of mesh node vector such that all base nodes "
        "are sorted before all nonlinear nodes.\n\n"
        "OpenGeoSys-6 software, version " +
            GitInfoLib::GitInfo::ogs_version +
            ".\n"
            "Copyright (c) 2012-2025, OpenGeoSys Community "
            "(http://www.opengeosys.org)",
        ' ', GitInfoLib::GitInfo::ogs_version);

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
            INFO("Method 0: Reversing order of nodes will be checked again.");
        }
        INFO(
            "Method: Reversing order of nodes unless it is considered correct "
            "by the OGS6 standard, i.e. such that det(J) > 0, where J is the "
            "Jacobian of the global-to-local coordinate transformation.");
        int const mesh_space_dimension =
            MeshLib::getSpaceDimension(mesh->getNodes());
        reverseNodeOrder(
            const_cast<std::vector<MeshLib::Element*>&>(mesh->getElements()),
            mesh_space_dimension, forced);
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

    MeshLib::IO::writeMeshToFile(*mesh, output_mesh_arg.getValue());

    INFO("VTU file written.");

    return EXIT_SUCCESS;
}
