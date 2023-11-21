/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <typeinfo>

#include "GenericNaturalBoundaryConditionLocalAssembler.h"
#include "NumLib/DOF/DOFTableUtil.h"
#include "NumLib/Fem/InitShapeMatrices.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "ParameterLib/MeshNodeParameter.h"
#include "ParameterLib/Parameter.h"

namespace ProcessLib
{
struct NeumannBoundaryConditionData final
{
    ParameterLib::Parameter<double> const& neumann_bc_parameter;
    ParameterLib::Parameter<double> const* const integral_measure;
};

template <typename ShapeFunction, int GlobalDim>
class NeumannBoundaryConditionLocalAssembler final
    : public GenericNaturalBoundaryConditionLocalAssembler<ShapeFunction,
                                                           GlobalDim>
{
    using Base =
        GenericNaturalBoundaryConditionLocalAssembler<ShapeFunction, GlobalDim>;
    using ShapeMatricesType = ShapeMatrixPolicyType<ShapeFunction, GlobalDim>;
    using NodalVectorType = typename Base::NodalVectorType;

public:
    /// The neumann_bc_value factor is directly integrated into the local
    /// element matrix.
    NeumannBoundaryConditionLocalAssembler(
        MeshLib::Element const& e,
        std::size_t const local_matrix_size,
        NumLib::GenericIntegrationMethod const& integration_method,
        bool const is_axially_symmetric,
        NeumannBoundaryConditionData const& data)
        : Base(e, is_axially_symmetric, integration_method),
          _data(data),
          _local_rhs(local_matrix_size)
    {
    }

    void assemble(std::size_t const id,
                  NumLib::LocalToGlobalIndexMap const& dof_table_boundary,
                  double const t, std::vector<GlobalVector*> const& /*x*/,
                  int const /*process_id*/, GlobalMatrix& /*K*/,
                  GlobalVector& b, GlobalMatrix* /*Jac*/) override
    {
        _local_rhs.setZero();

        unsigned const n_integration_points =
            Base::_integration_method.getNumberOfPoints();

        NodalVectorType parameter_node_values;
        if (typeid(ParameterLib::MeshNodeParameter<double> const) ==
            typeid(_data.neumann_bc_parameter))
        {
            // Get element nodes for the interpolation from nodes to integration
            // point.
            parameter_node_values =
                _data.neumann_bc_parameter
                    .getNodalValuesOnElement(Base::_element, t)
                    .template topRows<
                        ShapeFunction::MeshElement::n_all_nodes>();
        }

        double integral_measure = 1.0;
        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            auto const& ip_data = Base::_ns_and_weights[ip];
            auto const& N = ip_data.N;
            auto const& w = ip_data.weight;

            ParameterLib::SpatialPosition const position{
                std::nullopt, Base::_element.getID(), ip,
                MathLib::Point3d(
                    NumLib::interpolateCoordinates<ShapeFunction,
                                                   ShapeMatricesType>(
                        Base::_element, N))};

            if (_data.integral_measure)
            {
                integral_measure = (*_data.integral_measure)(t, position)[0];
            }
            if (typeid(ParameterLib::MeshNodeParameter<double> const) ==
                typeid(_data.neumann_bc_parameter))
            {
                _local_rhs.noalias() +=
                    N * parameter_node_values.dot(N) * w * integral_measure;
            }
            else
            {
                auto const value = _data.neumann_bc_parameter(t, position)[0];
                _local_rhs.noalias() += N * value * w * integral_measure;
            }
        }

        auto const indices = NumLib::getIndices(id, dof_table_boundary);
        b.add(indices, _local_rhs);
    }

private:
    NeumannBoundaryConditionData const& _data;

    NodalVectorType _local_rhs;

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

}  // namespace ProcessLib
