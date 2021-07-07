/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "MeshLib/PropertyVector.h"
#include "NumLib/DOF/DOFTableUtil.h"
#include "NumLib/Function/Interpolation.h"
#include "ParameterLib/MeshNodeParameter.h"

#include "GenericNaturalBoundaryConditionLocalAssembler.h"

namespace ProcessLib
{
struct VariableDependentNeumannBoundaryConditionData
{
    ParameterLib::Parameter<double> const& constant;
    ParameterLib::Parameter<double> const& coefficient_current_variable;
    ParameterLib::Parameter<double> const& coefficient_other_variable;
    ParameterLib::Parameter<double> const& coefficient_mixed_variables;
    // Used for mapping boundary nodes to bulk nodes.
    NumLib::LocalToGlobalIndexMap const& dof_table_boundary_other_variable;
};

template <typename ShapeFunction, typename IntegrationMethod, int GlobalDim>
class VariableDependentNeumannBoundaryConditionLocalAssembler final
    : public GenericNaturalBoundaryConditionLocalAssembler<
          ShapeFunction, IntegrationMethod, GlobalDim>
{
    using Base = GenericNaturalBoundaryConditionLocalAssembler<
        ShapeFunction, IntegrationMethod, GlobalDim>;
    using NodalVectorType = typename Base::NodalVectorType;
    using NodalMatrixType = typename Base::NodalMatrixType;

public:
    /// The neumann_bc_term factor is directly integrated into the local
    /// element matrix.
    VariableDependentNeumannBoundaryConditionLocalAssembler(
        MeshLib::Element const& e,
        std::size_t const local_matrix_size,
        bool const is_axially_symmetric,
        unsigned const integration_order,
        VariableDependentNeumannBoundaryConditionData const& data)
        : Base(e, is_axially_symmetric, integration_order),
          _data(data),
          _local_matrix_size(local_matrix_size)
    {
    }

    void assemble(std::size_t const mesh_item_id,
                  NumLib::LocalToGlobalIndexMap const& dof_table_boundary,
                  double const t, std::vector<GlobalVector*> const& x,
                  int const process_id, GlobalMatrix& /*K*/, GlobalVector& b,
                  GlobalMatrix* /*Jac*/) override
    {
        NodalVectorType _local_rhs(_local_matrix_size);
        _local_rhs.setZero();
        // Get element nodes for the interpolation from nodes to
        // integration point.
        NodalVectorType const constant_node_values =
            _data.constant.getNodalValuesOnElement(Base::_element, t)
                .template topRows<ShapeFunction::MeshElement::n_all_nodes>();
        NodalVectorType const coefficient_current_variable_node_values =
            _data.coefficient_current_variable
                .getNodalValuesOnElement(Base::_element, t)
                .template topRows<ShapeFunction::MeshElement::n_all_nodes>();
        NodalVectorType const coefficient_other_variable_node_values =
            _data.coefficient_other_variable
                .getNodalValuesOnElement(Base::_element, t)
                .template topRows<ShapeFunction::MeshElement::n_all_nodes>();
        NodalVectorType const coefficient_mixed_variables_node_values =
            _data.coefficient_mixed_variables
                .getNodalValuesOnElement(Base::_element, t)
                .template topRows<ShapeFunction::MeshElement::n_all_nodes>();
        unsigned const n_integration_points =
            Base::_integration_method.getNumberOfPoints();

        auto const indices_current_variable =
            NumLib::getIndices(mesh_item_id, dof_table_boundary);
        auto const indices_other_variable = NumLib::getIndices(
            mesh_item_id, _data.dof_table_boundary_other_variable);
        std::vector<double> const local_current_variable =
            x[process_id]->get(indices_current_variable);
        std::vector<double> const local_other_variable =
            x[process_id]->get(indices_other_variable);

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            auto const& n_and_weight = Base::_ns_and_weights[ip];
            auto const& N = n_and_weight.N;
            auto const& w = n_and_weight.weight;

            double current_variable_int_pt = 0.0;
            double other_variable_int_pt = 0.0;

            NumLib::shapeFunctionInterpolate(local_current_variable, N,
                                             current_variable_int_pt);
            NumLib::shapeFunctionInterpolate(local_other_variable, N,
                                             other_variable_int_pt);
            NodalVectorType const neumann_node_values =
                constant_node_values +
                coefficient_current_variable_node_values *
                    current_variable_int_pt +
                coefficient_other_variable_node_values * other_variable_int_pt +
                coefficient_mixed_variables_node_values *
                    current_variable_int_pt * other_variable_int_pt;
            _local_rhs.noalias() += N * neumann_node_values.dot(N) * w;
        }

        b.add(indices_current_variable, _local_rhs);
    }

private:
    VariableDependentNeumannBoundaryConditionData const& _data;
    std::size_t const _local_matrix_size;
};

}  // namespace ProcessLib
