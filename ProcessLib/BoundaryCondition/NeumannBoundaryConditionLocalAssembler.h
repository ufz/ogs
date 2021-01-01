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

#include "GenericNaturalBoundaryConditionLocalAssembler.h"
#include "NumLib/DOF/DOFTableUtil.h"
#include "NumLib/Fem/InitShapeMatrices.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "ParameterLib/Parameter.h"

namespace ProcessLib
{
struct NeumannBoundaryConditionData final
{
    ParameterLib::Parameter<double> const& neumann_bc_parameter;
    ParameterLib::Parameter<double> const* const integral_measure;
};

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
class NeumannBoundaryConditionLocalAssembler final
    : public GenericNaturalBoundaryConditionLocalAssembler<
          ShapeFunction, IntegrationMethod, GlobalDim>
{
    using Base = GenericNaturalBoundaryConditionLocalAssembler<
        ShapeFunction, IntegrationMethod, GlobalDim>;
    using ShapeMatricesType = ShapeMatrixPolicyType<ShapeFunction, GlobalDim>;
    using NodalVectorType = typename Base::NodalVectorType;

public:
    /// The neumann_bc_value factor is directly integrated into the local
    /// element matrix.
    NeumannBoundaryConditionLocalAssembler(
        MeshLib::Element const& e,
        std::size_t const local_matrix_size,
        bool const is_axially_symmetric,
        unsigned const integration_order,
        NeumannBoundaryConditionData const& data)
        : Base(e, is_axially_symmetric, integration_order),
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

        // Get element nodes for the interpolation from nodes to integration
        // point.
        NodalVectorType parameter_node_values =
            _data.neumann_bc_parameter
                .getNodalValuesOnElement(Base::_element, t)
                .template topRows<ShapeFunction::MeshElement::n_all_nodes>();

        double integral_measure = 1.0;
        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            auto const& ip_data = Base::_ns_and_weights[ip];
            auto const& N = ip_data.N;
            auto const& w = ip_data.weight;

            ParameterLib::SpatialPosition const position{
                boost::none, Base::_element.getID(), ip,
                MathLib::Point3d(
                    NumLib::interpolateCoordinates<ShapeFunction,
                                                   ShapeMatricesType>(
                        Base::_element, N))};

            if (_data.integral_measure)
            {
                integral_measure = (*_data.integral_measure)(t, position)[0];
            }
            _local_rhs.noalias() +=
                N * parameter_node_values.dot(N) * w * integral_measure;
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
