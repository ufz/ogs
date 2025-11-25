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
#include <concepts>
#include <functional>
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
#include "MeshLib/MeshEnums.h"
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

    Eigen::MatrixXd J = Eigen::MatrixXd::Zero(Dim, Dim);
    assert(e.getNumberOfNodes() == GradShapeFunction::NPOINTS);
    for (unsigned k = 0; k < GradShapeFunction::NPOINTS; k++)
    {
        const MathLib::Point3d& mapped_pt =
            ele_local_coord.getMappedCoordinates(k);
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

template <typename T>
concept ShapeFunction = requires(double* xi, Eigen::Map<Eigen::VectorXd> dN) {
    { T::DIM } -> std::convertible_to<int>;
    { T::NPOINTS } -> std::convertible_to<int>;
    T::computeGradShapeFunction(xi, dN);
};

struct ElementReorderConfigBase
{
    std::function<bool(MeshLib::Element&, int, bool)> is_node_ordering_correct;
    std::function<void(MeshLib::Element&, const std::vector<MeshLib::Node*>&)>
        reorder_element_nodes;
};

template <ShapeFunction ShapeFunc, int Dim>
ElementReorderConfigBase makeElementConfig(
    std::array<double, Dim> xi,
    std::function<void(MeshLib::Element&, const std::vector<MeshLib::Node*>&)>
        reorder)
{
    return ElementReorderConfigBase{
        [xi](MeshLib::Element& e, int mesh_space_dimension,
             bool check_reordered)
        {
            return checkJacobianDeterminant<ShapeFunc, Dim>(
                e, mesh_space_dimension, xi, check_reordered);
        },
        reorder};
}

static const std::array<ElementReorderConfigBase,
                        static_cast<int>(MeshLib::CellType::enum_length)>
    element_configs_array = []
{
    auto swap_nodes_i_j = [](MeshLib::Element& element,
                             const std::vector<MeshLib::Node*>& nodes,
                             unsigned const i, unsigned const j)
    {
        element.setNode(i, nodes[j]);
        element.setNode(j, nodes[i]);
    };

    auto order_nodes_quadratic_quad =
        [swap_nodes_i_j](MeshLib::Element& element,
                         const std::vector<MeshLib::Node*>& nodes)
    {
        swap_nodes_i_j(element, nodes, 1, 3);
        swap_nodes_i_j(element, nodes, 5, 6);
        swap_nodes_i_j(element, nodes, 4, 7);
    };

    std::array<ElementReorderConfigBase,
               static_cast<int>(MeshLib::CellType::enum_length)>
        arr{};

    arr[static_cast<int>(MeshLib::CellType::LINE2)] =
        makeElementConfig<NumLib::ShapeLine2, 1>(
            {0.5},
            [&swap_nodes_i_j](MeshLib::Element& element,
                              const std::vector<MeshLib::Node*>& nodes)
            { swap_nodes_i_j(element, nodes, 0, 1); });

    arr[static_cast<int>(MeshLib::CellType::LINE3)] =
        makeElementConfig<NumLib::ShapeLine3, 1>(
            {0.5},
            [&swap_nodes_i_j](MeshLib::Element& element,
                              const std::vector<MeshLib::Node*>& nodes)
            { swap_nodes_i_j(element, nodes, 0, 1); });

    arr[static_cast<int>(MeshLib::CellType::TRI3)] =
        makeElementConfig<NumLib::ShapeTri3, 2>(
            {1.0 / 3.0, 1.0 / 3.0},
            [&swap_nodes_i_j](MeshLib::Element& element,
                              const std::vector<MeshLib::Node*>& nodes)
            { swap_nodes_i_j(element, nodes, 1, 2); });

    arr[static_cast<int>(MeshLib::CellType::TRI6)] =
        makeElementConfig<NumLib::ShapeTri6, 2>(
            {1.0 / 3.0, 1.0 / 3.0},
            [&swap_nodes_i_j](MeshLib::Element& element,
                              const std::vector<MeshLib::Node*>& nodes)
            {
                swap_nodes_i_j(element, nodes, 1, 2);
                swap_nodes_i_j(element, nodes, 3, 5);
            });

    arr[static_cast<int>(MeshLib::CellType::QUAD4)] =
        makeElementConfig<NumLib::ShapeQuad4, 2>(
            {0.0, 0.0},
            [&swap_nodes_i_j](MeshLib::Element& element,
                              const std::vector<MeshLib::Node*>& nodes)
            { swap_nodes_i_j(element, nodes, 0, 2); });

    arr[static_cast<int>(MeshLib::CellType::QUAD8)] =
        makeElementConfig<NumLib::ShapeQuad8, 2>(
            {0.0, 0.0},
            [&order_nodes_quadratic_quad](
                MeshLib::Element& element,
                const std::vector<MeshLib::Node*>& nodes)
            { order_nodes_quadratic_quad(element, nodes); });
    arr[static_cast<int>(MeshLib::CellType::QUAD9)] =
        makeElementConfig<NumLib::ShapeQuad9, 2>(
            {0.0, 0.0},
            [&order_nodes_quadratic_quad](
                MeshLib::Element& element,
                const std::vector<MeshLib::Node*>& nodes)
            { order_nodes_quadratic_quad(element, nodes); });
    arr[static_cast<int>(MeshLib::CellType::TET4)] =
        makeElementConfig<NumLib::ShapeTet4, 3>(
            {0.25, 0.25, 0.25},
            [&swap_nodes_i_j](MeshLib::Element& element,
                              const std::vector<MeshLib::Node*>& nodes)
            { swap_nodes_i_j(element, nodes, 1, 2); });

    arr[static_cast<int>(MeshLib::CellType::TET10)] =
        makeElementConfig<NumLib::ShapeTet10, 3>(
            {0.25, 0.25, 0.25},
            [&swap_nodes_i_j](MeshLib::Element& element,
                              const std::vector<MeshLib::Node*>& nodes)
            {
                swap_nodes_i_j(element, nodes, 1, 2);
                swap_nodes_i_j(element, nodes, 4, 6);
                swap_nodes_i_j(element, nodes, 8, 9);
            });
    arr[static_cast<int>(MeshLib::CellType::PRISM6)] =
        makeElementConfig<NumLib::ShapePrism6, 3>(
            {1.0 / 3.0, 1.0 / 3.0, 0.5},
            [&swap_nodes_i_j](MeshLib::Element& element,
                              const std::vector<MeshLib::Node*>& nodes)
            {
                swap_nodes_i_j(element, nodes, 1, 2);
                swap_nodes_i_j(element, nodes, 4, 5);
            });

    arr[static_cast<int>(MeshLib::CellType::PRISM15)] =
        makeElementConfig<NumLib::ShapePrism15, 3>(
            {1.0 / 3.0, 1.0 / 3., 0.5},
            [&swap_nodes_i_j](MeshLib::Element& element,
                              const std::vector<MeshLib::Node*>& nodes)
            {
                swap_nodes_i_j(element, nodes, 0, 3);
                swap_nodes_i_j(element, nodes, 1, 4);
                swap_nodes_i_j(element, nodes, 2, 5);
                swap_nodes_i_j(element, nodes, 6, 9);
                swap_nodes_i_j(element, nodes, 7, 10);
                swap_nodes_i_j(element, nodes, 8, 11);
            });

    arr[static_cast<int>(MeshLib::CellType::PYRAMID5)] =
        makeElementConfig<NumLib::ShapePyra5, 3>(
            {0.25, 0.25, 0.5},
            [&swap_nodes_i_j](MeshLib::Element& element,
                              const std::vector<MeshLib::Node*>& nodes)
            { swap_nodes_i_j(element, nodes, 0, 2); });

    arr[static_cast<int>(MeshLib::CellType::PYRAMID13)] =
        makeElementConfig<NumLib::ShapePyra13, 3>(
            {0.25, 0.25, 0.5},
            [&swap_nodes_i_j](MeshLib::Element& element,
                              const std::vector<MeshLib::Node*>& nodes)
            {
                swap_nodes_i_j(element, nodes, 0, 2);
                swap_nodes_i_j(element, nodes, 9, 11);
                swap_nodes_i_j(element, nodes, 5, 6);
                swap_nodes_i_j(element, nodes, 7, 8);
            });
    arr[static_cast<int>(MeshLib::CellType::HEX8)] =
        makeElementConfig<NumLib::ShapeHex8, 3>(
            {0.5, 0.5, 0.5},
            [&swap_nodes_i_j](MeshLib::Element& element,
                              const std::vector<MeshLib::Node*>& nodes)
            {
                swap_nodes_i_j(element, nodes, 0, 2);
                swap_nodes_i_j(element, nodes, 4, 6);
            });

    arr[static_cast<int>(MeshLib::CellType::HEX20)] =
        makeElementConfig<NumLib::ShapeHex20, 3>(
            {0.5, 0.5, 0.5},
            [&swap_nodes_i_j](MeshLib::Element& element,
                              const std::vector<MeshLib::Node*>& nodes)
            {
                swap_nodes_i_j(element, nodes, 0, 2);
                swap_nodes_i_j(element, nodes, 4, 6);
                swap_nodes_i_j(element, nodes, 16, 18);
                swap_nodes_i_j(element, nodes, 8, 9);
                swap_nodes_i_j(element, nodes, 10, 11);
                swap_nodes_i_j(element, nodes, 12, 13);
                swap_nodes_i_j(element, nodes, 14, 15);
            });
    return arr;
}();

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

    for (auto* element : elements)
    {
        auto const cell_type = element->getCellType();

        if (cell_type == MeshLib::CellType::INVALID)
        {
            OGS_FATAL("Element type `{:s}' does not exist.",
                      MeshLib::CellType2String(cell_type));
        }
        // This check is already done in MeshLib::readMeshFromFile(). Therefore,
        // this is just a safeguard.
        if (cell_type == MeshLib::CellType::HEX27 ||
            cell_type == MeshLib::CellType::PRISM18)
        {
            OGS_FATAL("Element type `{:s}' is not supported in OGS.",
                      MeshLib::CellType2String(cell_type));
        }

        if (cell_type == MeshLib::CellType::POINT1)
        {
            continue;
        }

        const auto& element_config =
            element_configs_array[static_cast<int>(cell_type)];

        if (element_config.is_node_ordering_correct(
                *element, mesh_space_dimension, false /*check_reordered*/))
        {
            continue;
        }

        // Save nodes before reordering
        const unsigned nElemNodes = element->getNumberOfNodes();
        std::vector<MeshLib::Node*> nodes(element->getNodes(),
                                          element->getNodes() + nElemNodes);

        double const element_volume_origin = element->computeVolume();

        element_config.reorder_element_nodes(*element, nodes);

        // Ensure that the element volume did not change.
        double const element_volume = element->computeVolume();
        //  We use a fixed tolerance here because for very small elements the
        //  machine epsilon might be too small.
        if (std::fabs(element_volume - element_volume_origin) /
                element_volume_origin >
            10 * std::numeric_limits<double>::epsilon())
        {
            OGS_FATAL(
                "Reordering the nodes of element {:d} failed as its volume "
                "changed from {:.20g} to {:.20g}, the difference is {:.20g} "
                "and the threshold is {:.20g}.",
                element->getID(), element_volume_origin, element_volume,
                std::fabs(element_volume_origin - element_volume),
                10 * std::numeric_limits<double>::epsilon());
        }
        // Element::computeVolume() uses the element vertecies to compute
        // the element volume, so the change of edge nodes are not
        // considered. Therefore, we need to additionally check the Jacobian
        // determinant here if forced is true.
        if (forced)
        {
            element_config.is_node_ordering_correct(
                *element, mesh_space_dimension, true /*check_reordered*/);
        }

        ++n_corrected_elements;
    }

    INFO("Corrected {:d} elements.", n_corrected_elements);
}

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
    enum class ExpectedCellType
    {
        INVALID = 0,
        POINT1 = 1,
        LINE2 = 2,
        LINE3 = 3,
        TRI3 = 4,
        TRI6 = 5,
        QUAD4 = 6,
        QUAD8 = 7,
        QUAD9 = 8,
        TET4 = 9,
        TET10 = 10,
        HEX8 = 11,
        HEX20 = 12,
        HEX27 = 13,
        PRISM6 = 14,
        PRISM15 = 15,
        PRISM18 = 16,
        PYRAMID5 = 17,
        PYRAMID13 = 18,
        enum_length
    };

    constexpr bool are_expected_cell_types =
        static_cast<int>(ExpectedCellType::enum_length) ==
            static_cast<int>(MeshLib::CellType::enum_length) &&
        static_cast<int>(ExpectedCellType::INVALID) ==
            static_cast<int>(MeshLib::CellType::INVALID) &&
        static_cast<int>(ExpectedCellType::POINT1) ==
            static_cast<int>(MeshLib::CellType::POINT1) &&
        static_cast<int>(ExpectedCellType::LINE2) ==
            static_cast<int>(MeshLib::CellType::LINE2) &&
        static_cast<int>(ExpectedCellType::LINE3) ==
            static_cast<int>(MeshLib::CellType::LINE3) &&
        static_cast<int>(ExpectedCellType::TRI3) ==
            static_cast<int>(MeshLib::CellType::TRI3) &&
        static_cast<int>(ExpectedCellType::TRI6) ==
            static_cast<int>(MeshLib::CellType::TRI6) &&
        static_cast<int>(ExpectedCellType::QUAD4) ==
            static_cast<int>(MeshLib::CellType::QUAD4) &&
        static_cast<int>(ExpectedCellType::QUAD8) ==
            static_cast<int>(MeshLib::CellType::QUAD8) &&
        static_cast<int>(ExpectedCellType::QUAD9) ==
            static_cast<int>(MeshLib::CellType::QUAD9) &&
        static_cast<int>(ExpectedCellType::TET4) ==
            static_cast<int>(MeshLib::CellType::TET4) &&
        static_cast<int>(ExpectedCellType::TET10) ==
            static_cast<int>(MeshLib::CellType::TET10) &&
        static_cast<int>(ExpectedCellType::HEX8) ==
            static_cast<int>(MeshLib::CellType::HEX8) &&
        static_cast<int>(ExpectedCellType::HEX20) ==
            static_cast<int>(MeshLib::CellType::HEX20) &&
        static_cast<int>(ExpectedCellType::HEX27) ==
            static_cast<int>(MeshLib::CellType::HEX27) &&
        static_cast<int>(ExpectedCellType::PRISM6) ==
            static_cast<int>(MeshLib::CellType::PRISM6) &&
        static_cast<int>(ExpectedCellType::PRISM15) ==
            static_cast<int>(MeshLib::CellType::PRISM15) &&
        static_cast<int>(ExpectedCellType::PRISM18) ==
            static_cast<int>(MeshLib::CellType::PRISM18) &&
        static_cast<int>(ExpectedCellType::PYRAMID5) ==
            static_cast<int>(MeshLib::CellType::PYRAMID5) &&
        static_cast<int>(ExpectedCellType::PYRAMID13) ==
            static_cast<int>(MeshLib::CellType::PYRAMID13);
    // This error should never occur unless someone has changed the
    // MeshLib::CellType enum class. These assertions ensure that the
    // array 'element_configs_array' is up to date.
    if (!are_expected_cell_types)
    {
        OGS_FATAL(
            "The enum class MeshLib::CellType has changed. Please adapt the "
            "array 'element_configs_array' in NodeReordering.cpp accordingly.");
    }

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
