/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
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

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
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
          data_(data),
          local_matrix_size_(local_matrix_size)
    {
    }

    void assemble(std::size_t const mesh_item_id,
                  NumLib::LocalToGlobalIndexMap const& dof_table_boundary,
                  double const t, std::vector<GlobalVector*> const& x,
                  int const process_id, GlobalMatrix& /*K*/, GlobalVector& b,
                  GlobalMatrix* /*Jac*/) override
    {
        NodalVectorType local_rhs_(local_matrix_size_);
        local_rhs_.setZero();
        // Get element nodes for the interpolation from nodes to
        // integration point.
        NodalVectorType const constant_node_values =
            data_.constant.getNodalValuesOnElement(Base::element_, t)
                .template topRows<ShapeFunction::MeshElement::n_all_nodes>();
        NodalVectorType const coefficient_current_variable_node_values =
            data_.coefficient_current_variable
                .getNodalValuesOnElement(Base::element_, t)
                .template topRows<ShapeFunction::MeshElement::n_all_nodes>();
        NodalVectorType const coefficient_other_variable_node_values =
            data_.coefficient_other_variable
                .getNodalValuesOnElement(Base::element_, t)
                .template topRows<ShapeFunction::MeshElement::n_all_nodes>();
        NodalVectorType const coefficient_mixed_variables_node_values =
            data_.coefficient_mixed_variables
                .getNodalValuesOnElement(Base::element_, t)
                .template topRows<ShapeFunction::MeshElement::n_all_nodes>();
        unsigned const n_integration_points =
            Base::integration_method_.getNumberOfPoints();

        auto const indices_current_variable =
            NumLib::getIndices(mesh_item_id, dof_table_boundary);
        auto const indices_other_variable = NumLib::getIndices(
            mesh_item_id, data_.dof_table_boundary_other_variable);
        std::vector<double> const local_current_variable =
            x[process_id]->get(indices_current_variable);
        std::vector<double> const local_other_variable =
            x[process_id]->get(indices_other_variable);

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            auto const& n_and_weight = Base::ns_and_weights_[ip];
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
            local_rhs_.noalias() += N * neumann_node_values.dot(N) * w;
        }

        b.add(indices_current_variable, local_rhs_);
    }

private:
    VariableDependentNeumannBoundaryConditionData const& data_;
    std::size_t const local_matrix_size_;
};

}  // namespace ProcessLib
