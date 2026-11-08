// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "ReverseNodeOrder.h"

#include <Eigen/LU>
#include <array>
#include <cassert>
#include <cmath>
#include <concepts>
#include <functional>
#include <limits>
#include <vector>

#include "BaseLib/Error.h"
#include "MeshLib/ElementCoordinatesMappingLocal.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/MeshEnums.h"
#include "MeshLib/Node.h"
#include "NumLib/Fem/ShapeFunction/ShapeHex20.h"
#include "NumLib/Fem/ShapeFunction/ShapeHex8.h"
#include "NumLib/Fem/ShapeFunction/ShapeLine2.h"
#include "NumLib/Fem/ShapeFunction/ShapeLine3.h"
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

namespace
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
static_assert(are_expected_cell_types,
              "The enum class MeshLib::CellType has changed. Please adapt the "
              "array 'element_configs_array' in ReverseNodeOrder.cpp "
              "accordingly.");

template <typename GradShapeFunction, int Dim>
bool checkJacobianDeterminant(MeshLib::Element const& e,
                              int const mesh_space_dimension,
                              std::array<double, Dim> const& xi,
                              bool const check_reordered)
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
            "reverseNodeOrder() failed: element {:d} has negative Jacobian "
            "determinant {:g}.",
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

const std::array<ElementReorderConfigBase,
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
            [swap_nodes_i_j](MeshLib::Element& element,
                             const std::vector<MeshLib::Node*>& nodes)
            { swap_nodes_i_j(element, nodes, 0, 1); });

    arr[static_cast<int>(MeshLib::CellType::LINE3)] =
        makeElementConfig<NumLib::ShapeLine3, 1>(
            {0.5},
            [swap_nodes_i_j](MeshLib::Element& element,
                             const std::vector<MeshLib::Node*>& nodes)
            { swap_nodes_i_j(element, nodes, 0, 1); });

    arr[static_cast<int>(MeshLib::CellType::TRI3)] =
        makeElementConfig<NumLib::ShapeTri3, 2>(
            {1.0 / 3.0, 1.0 / 3.0},
            [swap_nodes_i_j](MeshLib::Element& element,
                             const std::vector<MeshLib::Node*>& nodes)
            { swap_nodes_i_j(element, nodes, 1, 2); });

    arr[static_cast<int>(MeshLib::CellType::TRI6)] =
        makeElementConfig<NumLib::ShapeTri6, 2>(
            {1.0 / 3.0, 1.0 / 3.0},
            [swap_nodes_i_j](MeshLib::Element& element,
                             const std::vector<MeshLib::Node*>& nodes)
            {
                swap_nodes_i_j(element, nodes, 1, 2);
                swap_nodes_i_j(element, nodes, 3, 5);
            });

    arr[static_cast<int>(MeshLib::CellType::QUAD4)] =
        makeElementConfig<NumLib::ShapeQuad4, 2>(
            {0.0, 0.0},
            [swap_nodes_i_j](MeshLib::Element& element,
                             const std::vector<MeshLib::Node*>& nodes)
            { swap_nodes_i_j(element, nodes, 0, 2); });

    arr[static_cast<int>(MeshLib::CellType::QUAD8)] =
        makeElementConfig<NumLib::ShapeQuad8, 2>({0.0, 0.0},
                                                 order_nodes_quadratic_quad);
    arr[static_cast<int>(MeshLib::CellType::QUAD9)] =
        makeElementConfig<NumLib::ShapeQuad9, 2>({0.0, 0.0},
                                                 order_nodes_quadratic_quad);
    arr[static_cast<int>(MeshLib::CellType::TET4)] =
        makeElementConfig<NumLib::ShapeTet4, 3>(
            {0.25, 0.25, 0.25},
            [swap_nodes_i_j](MeshLib::Element& element,
                             const std::vector<MeshLib::Node*>& nodes)
            { swap_nodes_i_j(element, nodes, 1, 2); });

    arr[static_cast<int>(MeshLib::CellType::TET10)] =
        makeElementConfig<NumLib::ShapeTet10, 3>(
            {0.25, 0.25, 0.25},
            [swap_nodes_i_j](MeshLib::Element& element,
                             const std::vector<MeshLib::Node*>& nodes)
            {
                swap_nodes_i_j(element, nodes, 1, 2);
                swap_nodes_i_j(element, nodes, 4, 6);
                swap_nodes_i_j(element, nodes, 8, 9);
            });
    arr[static_cast<int>(MeshLib::CellType::PRISM6)] =
        makeElementConfig<NumLib::ShapePrism6, 3>(
            {1.0 / 3.0, 1.0 / 3.0, 0.5},
            [swap_nodes_i_j](MeshLib::Element& element,
                             const std::vector<MeshLib::Node*>& nodes)
            {
                swap_nodes_i_j(element, nodes, 1, 2);
                swap_nodes_i_j(element, nodes, 4, 5);
            });

    arr[static_cast<int>(MeshLib::CellType::PRISM15)] =
        makeElementConfig<NumLib::ShapePrism15, 3>(
            {1.0 / 3.0, 1.0 / 3., 0.5},
            [swap_nodes_i_j](MeshLib::Element& element,
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
            [swap_nodes_i_j](MeshLib::Element& element,
                             const std::vector<MeshLib::Node*>& nodes)
            { swap_nodes_i_j(element, nodes, 0, 2); });

    arr[static_cast<int>(MeshLib::CellType::PYRAMID13)] =
        makeElementConfig<NumLib::ShapePyra13, 3>(
            {0.25, 0.25, 0.5},
            [swap_nodes_i_j](MeshLib::Element& element,
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
            [swap_nodes_i_j](MeshLib::Element& element,
                             const std::vector<MeshLib::Node*>& nodes)
            {
                swap_nodes_i_j(element, nodes, 0, 2);
                swap_nodes_i_j(element, nodes, 4, 6);
            });

    arr[static_cast<int>(MeshLib::CellType::HEX20)] =
        makeElementConfig<NumLib::ShapeHex20, 3>(
            {0.5, 0.5, 0.5},
            [swap_nodes_i_j](MeshLib::Element& element,
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

void checkElementVolume(int const element_id, double const v_old,
                        double const v_new)
{
    //  We use a fixed tolerance here because for very small elements
    //  the machine epsilon might be too small.
    double const eps = 100.0 * std::numeric_limits<double>::epsilon();

    double const diff = std::abs(v_new - v_old);

    if (diff <= eps || diff / v_old <= eps)
    {
        return;
    }
    OGS_FATAL(
        "reverseNodeOrder() failed: reordering the nodes of element {:d} "
        "changed its volume from {:.20g} to {:.20g}. The absolute difference "
        "is {:.20g} and the relative difference is {:.20g} (thresholds are "
        "{:.20g} and {:.20g}).",
        element_id, v_old, v_new, diff, diff / v_old, eps, eps);
}

}  // namespace

namespace MeshToolsLib
{
std::size_t reverseNodeOrder(std::vector<MeshLib::Element*> const& elements,
                             int const mesh_space_dimension, bool const forced,
                             bool const volume_check)
{
    std::size_t n_corrected_elements = 0;

    for (auto* element : elements)
    {
        auto const cell_type = element->getCellType();

        if (cell_type == MeshLib::CellType::INVALID)
        {
            OGS_FATAL(
                "reverseNodeOrder() failed: element type `{:s}' does not "
                "exist.",
                MeshLib::CellType2String(cell_type));
        }
        // These types have no entry in `element_configs_array`, so reordering
        // them would call an empty `std::function`.
        if (cell_type == MeshLib::CellType::HEX27 ||
            cell_type == MeshLib::CellType::PRISM18)
        {
            OGS_FATAL(
                "reverseNodeOrder() failed: element type `{:s}' is not "
                "supported in OGS.",
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

        if (volume_check)
        {
            double const element_volume_origin = element->computeVolume();

            element_config.reorder_element_nodes(*element, nodes);

            // Ensure that the element volume did not change.
            double const element_volume = element->computeVolume();

            checkElementVolume(element->getID(), element_volume_origin,
                               element_volume);
        }
        else
        {
            element_config.reorder_element_nodes(*element, nodes);
        }

        // Element::computeVolume() uses the element vertices to compute
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

    return n_corrected_elements;
}

}  // namespace MeshToolsLib
