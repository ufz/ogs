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
#include "ProcessLib/Parameter/MeshNodeParameter.h"

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
    using NodalVectorType = typename Base::NodalVectorType;

public:
    /// The neumann_bc_value factor is directly integrated into the local
    /// element matrix.
    NonuniformNeumannBoundaryConditionLocalAssembler(
        MeshLib::Element const& e,
        std::size_t const local_matrix_size,
        bool const is_axially_symmetric,
        unsigned const integration_order,
        NonuniformNeumannBoundaryConditionData const& data)
        : Base(e, is_axially_symmetric, integration_order),
          _data(data),
          _local_rhs(local_matrix_size)
    {
    }

    void assemble(std::size_t const id,
                  NumLib::LocalToGlobalIndexMap const& dof_table_boundary,
                  double const t, const GlobalVector& /*x*/,
                  GlobalMatrix& /*K*/, GlobalVector& b,
                  GlobalMatrix* /*Jac*/) override
    {
        _local_rhs.setZero();

        MeshNodeParameter<double> neumann_values{"NeumannValues", _data.values};
        // Get element nodes for the interpolation from nodes to
        // integration point.
        NodalVectorType parameter_node_values =
            neumann_values.getNodalValuesOnElement(Base::_element, t);

        unsigned const n_integration_points =
            Base::_integration_method.getNumberOfPoints();
        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            auto const& n_and_weight = Base::_ns_and_weights[ip];
            auto const& N = n_and_weight.N;
            auto const& w = n_and_weight.weight;
            _local_rhs.noalias() += N * parameter_node_values.dot(N) * w;
        }

        auto const indices = NumLib::getIndices(id, dof_table_boundary);
        b.add(indices, _local_rhs);
    }

private:
    NonuniformNeumannBoundaryConditionData const& _data;
    NodalVectorType _local_rhs;

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

}  // namespace ProcessLib
