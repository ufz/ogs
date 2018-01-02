/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "MeshLib/PropertyVector.h"
#include "NumLib/DOF/DOFTableUtil.h"
#include "NumLib/Function/Interpolation.h"

#include "GenericNonuniformNaturalBoundaryConditionLocalAssembler.h"

namespace ProcessLib
{
struct NonuniformNeumannBoundaryConditionData
{
    /* TODO use Parameter */
    MeshLib::PropertyVector<double> const& values;

    // Used for mapping boundary nodes to bulk nodes.
    std::size_t bulk_mesh_id;
    MeshLib::PropertyVector<std::size_t> const& mapping_to_bulk_nodes;
    NumLib::LocalToGlobalIndexMap const& dof_table_bulk;
    int const variable_id_bulk;
    int const component_id_bulk;
};

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
class NonuniformNeumannBoundaryConditionLocalAssembler final
    : public GenericNonuniformNaturalBoundaryConditionLocalAssembler<
          ShapeFunction, IntegrationMethod, GlobalDim>
{
    using Base = GenericNonuniformNaturalBoundaryConditionLocalAssembler<
        ShapeFunction, IntegrationMethod, GlobalDim>;

public:
    /// The neumann_bc_value factor is directly integrated into the local
    /// element matrix.
    NonuniformNeumannBoundaryConditionLocalAssembler(
        MeshLib::Element const& e,
        std::size_t const local_matrix_size,
        bool is_axially_symmetric,
        unsigned const integration_order,
        NonuniformNeumannBoundaryConditionData const& data)
        : Base(e, is_axially_symmetric, integration_order),
          _data(data),
          _local_rhs(local_matrix_size)
    {
    }

    void assemble(std::size_t const id,
                  NumLib::LocalToGlobalIndexMap const& dof_table_boundary,
                  double const /*t*/, const GlobalVector& /*x*/,
                  GlobalMatrix& /*K*/, GlobalVector& b) override
    {
        _local_rhs.setZero();

        auto indices = NumLib::getIndices(id, dof_table_boundary);

        // TODO lots of heap allocations.
        std::vector<double> neumann_param_nodal_values_local;
        neumann_param_nodal_values_local.reserve(indices.size());
        for (auto i : indices)
        {
            neumann_param_nodal_values_local.push_back(
                _data.values.getComponent(i, 0));
        }

        auto const n_integration_points = Base::_ns_and_weights.size();
        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            auto const& n_and_weight = Base::_ns_and_weights[ip];
            double flux;
            NumLib::shapeFunctionInterpolate(neumann_param_nodal_values_local,
                                             n_and_weight.N, flux);
            _local_rhs.noalias() +=
                n_and_weight.N * (n_and_weight.weight * flux);
        }

        // map boundary dof indices to bulk dof indices
        for (auto& i : indices)
        {
            auto const bulk_node_id =
                _data.mapping_to_bulk_nodes.getComponent(i, 0);

            MeshLib::Location const l{
                _data.bulk_mesh_id, MeshLib::MeshItemType::Node, bulk_node_id};

            i = _data.dof_table_bulk.getGlobalIndex(l, _data.variable_id_bulk,
                                                    _data.component_id_bulk);
            assert(i != NumLib::MeshComponentMap::nop);
        }
        b.add(indices, _local_rhs);
    }

private:
    NonuniformNeumannBoundaryConditionData const& _data;
    typename Base::NodalVectorType _local_rhs;

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

}  // namespace ProcessLib
