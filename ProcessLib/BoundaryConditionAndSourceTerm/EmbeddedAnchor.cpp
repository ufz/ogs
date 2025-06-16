/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "EmbeddedAnchor.h"

#include <Eigen/Core>
#include <cassert>
#include <numeric>
#include <range/v3/view/enumerate.hpp>
#include <unsupported/Eigen/KroneckerProduct>

#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "NumLib/DOF/DOFTableUtil.h"
#include "NumLib/Fem/FiniteElement/ElementTraitsLagrange.h"
#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/InitShapeMatrices.h"
#include "NumLib/Fem/Interpolation.h"
#include "NumLib/Fem/ShapeFunction/ShapeQuad4.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "ProcessLib/BoundaryConditionAndSourceTerm/Utils/CreateLocalAssemblers.h"

namespace ProcessLib
{
// The local indices are, due to the nature of the DOF table, all even
// for the first node and odd for the second node.
template <int GlobalDim>
auto nodeLocalIndices(std::size_t const i)
{
    return Eigen::seqN(i, Eigen::fix<GlobalDim>, Eigen::fix<2>);
}

template <int GlobalDim>
Eigen::RowVectorXd computeShapeMatrix(
    MeshLib::Element const& bulk_element,
    std::span<double const, 3> const nat_coords)
{
    using ShapeMatricesType =
        EigenDynamicShapeMatrixPolicy<void /* dummy */, -1 /* dynamic order */>;

    typename ShapeMatricesType::ShapeMatrices shape_matrix(
        bulk_element.getDimension(), GlobalDim,
        bulk_element.getNumberOfNodes());

    bool matched = false;

    BaseLib::TMP::foreach<NumLib::AllElementTraitsLagrange>(
        [&](auto* tag)
        {
            using ElementTrait = std::remove_pointer_t<decltype(tag)>;
            using ShapeFunction = typename ElementTrait::ShapeFunction;
            using MeshElement = typename ElementTrait::Element;

            if (!matched &&
                dynamic_cast<MeshElement const*>(&bulk_element) != nullptr)
            {
                auto const fe = NumLib::createIsoparametricFiniteElement<
                    ShapeFunction, ShapeMatricesType>(bulk_element);

                fe.template computeShapeFunctions<NumLib::ShapeMatrixType::N>(
                    nat_coords.data(), shape_matrix, GlobalDim,
                    false /* axisymmetric */);

                matched = true;
            }
        });

    if (!matched)
    {
        OGS_FATAL("Element type '{:s}' is not supported as anchor element.",
                  MeshLib::CellType2String(bulk_element.getCellType()));
    }

    return shape_matrix.N;
}

template <int GlobalDim>
std::tuple<Eigen::VectorXd, Eigen::MatrixXd> assembleLocalBJac(
    Eigen::Vector<double, GlobalDim> const& f,
    Eigen::Matrix<double, GlobalDim, GlobalDim> const& Df,
    std::vector<Eigen::RowVectorXd> const& shape_matrices,
    std::size_t const num_dof,
    std::array<std::size_t, 2> const& nodes_per_element)
{
    // Signs for the two nodes alternate.
    constexpr auto even_odd_sign = [](std::size_t const n)
    { return (n % 2 == 0) ? 1.0 : -1.0; };

    Eigen::VectorXd local_rhs(num_dof);
    Eigen::MatrixXd local_Jac(num_dof, num_dof);
    local_rhs.setZero();
    local_Jac.setZero();

    for (std::size_t node_idx = 0; node_idx < 2; ++node_idx)
    {
        // Select the correct block according to the node index.
        // Replicate shape matrices of the corresponding element to GlobalDim
        // dimensions. Multiply component-wise with the force vector f.
        local_rhs.segment(node_idx * GlobalDim * nodes_per_element[0],
                          GlobalDim * nodes_per_element[node_idx]) =
            even_odd_sign(node_idx) *
            shape_matrices[node_idx]
                .transpose()
                .replicate<GlobalDim, 1>()
                .cwiseProduct(f.transpose()
                                  .replicate(nodes_per_element[node_idx], 1)
                                  .reshaped());
        for (std::size_t node_idx_inner = 0; node_idx_inner < 2;
             ++node_idx_inner)
        {
            Eigen::MatrixXd const& ones = Eigen::MatrixXd::Ones(
                nodes_per_element[node_idx], nodes_per_element[node_idx_inner]);

            // Select the correct block according to the node index.
            // Replicate shape matrices of the corresponding element to
            // GlobalDim dimensions. Multiply component-wise with the derivative
            // Df of force vector f.
            local_Jac.block(node_idx * GlobalDim * nodes_per_element[0],
                            node_idx_inner * GlobalDim * nodes_per_element[0],
                            GlobalDim * nodes_per_element[node_idx],
                            GlobalDim * nodes_per_element[node_idx_inner]) =
                (even_odd_sign(node_idx) *
                 shape_matrices[node_idx]
                     .transpose()
                     .replicate<GlobalDim, 1>() *
                 even_odd_sign(node_idx_inner) *
                 shape_matrices[node_idx_inner].replicate<1, GlobalDim>())
                    .cwiseProduct(kroneckerProduct(Df, ones));
        }
    }
    return {local_rhs, local_Jac};
}

template <int GlobalDim>
EmbeddedAnchor<GlobalDim>::EmbeddedAnchor(
    MeshLib::Mesh const& bulk_mesh,
    NumLib::LocalToGlobalIndexMap const& dof_table_bulk,
    std::size_t const source_term_mesh_id,
    MeshLib::Mesh const& st_mesh,
    const int variable_id,
    ParameterLib::Parameter<double> const& parameter)
    : dof_table_bulk_(dof_table_bulk),
      bulk_mesh_(bulk_mesh),
      source_term_mesh_id_(source_term_mesh_id),
      st_mesh_(st_mesh),
      variable_id_(variable_id),
      component_ids_(
          []
          {
              std::array<int, GlobalDim> arr{};
              std::iota(arr.begin(), arr.end(), 0);
              return arr;
          }()),
      parameter_(parameter)
{
    DBUG("Create EmbeddedAnchor.");
    std::string_view const bulk_element_ids_string = "bulk_element_ids";
    bulk_element_ids_ =
        st_mesh_.getProperties().template getPropertyVector<std::size_t>(
            bulk_element_ids_string);
    assert(bulk_element_ids_ != nullptr);

    std::string_view const natural_coordinates_string = "natural_coordinates";
    natural_coordinates_ =
        st_mesh_.getProperties().template getPropertyVector<double>(
            natural_coordinates_string);
    assert(natural_coordinates_ != nullptr);
}

template <int GlobalDim>
void EmbeddedAnchor<GlobalDim>::
    getShapeMatricesAndGlobalIndicesAndDisplacements(
        MeshLib::Element const* const anchor_element,
        std::array<std::size_t, 2>& nodes_per_element,
        std::vector<Eigen::RowVectorXd>& shape_matrices,
        std::vector<GlobalIndexType>& global_indices,
        Eigen::Vector<double, 2 * GlobalDim>& local_x,
        GlobalVector const& x,
        ParameterLib::SpatialPosition& pos) const
{
    for (auto&& [node_idx, anchor_node_id] : ranges::views::enumerate(
             anchor_element->nodes() | MeshLib::views::ids))
    {
        // used pos in integrate function will be of the last node in anchor
        // element, i,e, only one uniform stiffness is supported
        auto const bulk_element_id = (*bulk_element_ids_)[anchor_node_id];
        pos.setElementID(bulk_element_id);

        auto const& bulk_element = *bulk_mesh_.getElement(bulk_element_id);
        nodes_per_element[node_idx] = bulk_element.nodes().size();

        std::span<double const, 3> const nat_coords(
            &(*natural_coordinates_)[3 * anchor_node_id], 3);
        auto N = computeShapeMatrix<GlobalDim>(bulk_element, nat_coords);
        shape_matrices.push_back(std::move(N));

        for (int component = 0; component < GlobalDim; ++component)
        {
            auto const& global_indices_percomponent = getIndices(
                variable_id_, component, bulk_element_id, dof_table_bulk_);

            Eigen::VectorXd const local_x_element_percomponent =
                MathLib::toVector(x.get(global_indices_percomponent));

            global_indices.insert(global_indices.cend(),
                                  global_indices_percomponent.begin(),
                                  global_indices_percomponent.end());

            // this maps the displacement of each component of all nodes
            // of the element to the desired anchor location in that element
            // where all even local_x entries are for the first node and all odd
            // local_x entries are for the second node
            NumLib::shapeFunctionInterpolate(
                local_x_element_percomponent,
                shape_matrices[node_idx],
                local_x[nodeLocalIndices<GlobalDim>(node_idx)[component]]);
        }
    }
}

template <int GlobalDim>
void EmbeddedAnchor<GlobalDim>::integrate(const double t, GlobalVector const& x,
                                          GlobalVector& b,
                                          GlobalMatrix* jac) const
{
    DBUG("Assemble EmbeddedAnchor.");

    using GlobalDimVector = Eigen::Vector<double, GlobalDim>;
    using GlobalDimMatrix = Eigen::Matrix<double, GlobalDim, GlobalDim>;

    for (MeshLib::Element const* const anchor_element : st_mesh_.getElements())
    {
        std::vector<GlobalIndexType> global_indices;
        Eigen::Vector<double, 2 * GlobalDim> local_x;
        std::vector<Eigen::RowVectorXd> shape_matrices;
        std::array<std::size_t, 2> nodes_per_element;
        ParameterLib::SpatialPosition pos;
        getShapeMatricesAndGlobalIndicesAndDisplacements(
            anchor_element, nodes_per_element, shape_matrices, global_indices,
            local_x, x, pos);

        auto node_coords = [anchor_element](int const i)
        { return anchor_element->getNode(i)->asEigenVector3d(); };
        GlobalDimVector const l_original =
            (node_coords(1) - node_coords(0)).template head<GlobalDim>();
        double const l_original_norm = l_original.norm();

        // Displacement in the two nodes.
        auto u = [&local_x](int const i)
        { return local_x(nodeLocalIndices<GlobalDim>(i)); };
        GlobalDimVector const l = l_original + u(1) - u(0);

        double const K = parameter_(t, pos)[0];
        GlobalDimVector const f = l_original / l_original_norm * K *
                                  (l.norm() - l_original_norm) /
                                  l_original_norm;

        GlobalDimMatrix const Df = l_original / l_original_norm * K *
                                   l.transpose() / l.norm() / l_original_norm;

        auto const& [local_rhs, local_Jac] = assembleLocalBJac<GlobalDim>(
            f, Df, shape_matrices, global_indices.size(), nodes_per_element);

        b.add(global_indices, local_rhs);
        if (jac)
        {
            jac->add({global_indices, global_indices}, local_Jac);
        }
    }
}

template class EmbeddedAnchor<2>;
template class EmbeddedAnchor<3>;

}  // namespace ProcessLib
