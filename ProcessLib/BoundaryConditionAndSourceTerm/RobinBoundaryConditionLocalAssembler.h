/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "GenericNaturalBoundaryConditionLocalAssembler.h"
#include "NumLib/DOF/DOFTableUtil.h"
#include "ParameterLib/Parameter.h"

namespace ProcessLib
{
struct RobinBoundaryConditionData final
{
    ParameterLib::Parameter<double> const& alpha;
    ParameterLib::Parameter<double> const& u_0;
    ParameterLib::Parameter<double> const* const integral_measure;
};

template <typename ShapeFunction, typename IntegrationMethod, int GlobalDim>
class RobinBoundaryConditionLocalAssembler final
    : public GenericNaturalBoundaryConditionLocalAssembler<
          ShapeFunction, IntegrationMethod, GlobalDim>
{
    using Base = GenericNaturalBoundaryConditionLocalAssembler<
        ShapeFunction, IntegrationMethod, GlobalDim>;
    using ShapeMatricesType = ShapeMatrixPolicyType<ShapeFunction, GlobalDim>;

public:
    RobinBoundaryConditionLocalAssembler(MeshLib::Element const& e,
                                         std::size_t const local_matrix_size,
                                         bool is_axially_symmetric,
                                         unsigned const integration_order,
                                         RobinBoundaryConditionData const& data)
        : Base(e, is_axially_symmetric, integration_order),
          _data(data),
          _local_K(local_matrix_size, local_matrix_size),
          _local_rhs(local_matrix_size)
    {
    }

    // TODO also implement derivative for Jacobian in Newton scheme.
    void assemble(std::size_t const id,
                  NumLib::LocalToGlobalIndexMap const& dof_table_boundary,
                  double const t, std::vector<GlobalVector*> const& /*x*/,
                  int const /*process_id*/, GlobalMatrix& K, GlobalVector& b,
                  GlobalMatrix* /*Jac*/) override
    {
        _local_K.setZero();
        _local_rhs.setZero();

        unsigned const n_integration_points =
            Base::_integration_method.getNumberOfPoints();

        typename Base::NodalVectorType const alpha =
            _data.alpha.getNodalValuesOnElement(Base::_element, t)
                .template topRows<ShapeFunction::MeshElement::n_all_nodes>();
        typename Base::NodalVectorType const u_0 =
            _data.u_0.getNodalValuesOnElement(Base::_element, t)
                .template topRows<ShapeFunction::MeshElement::n_all_nodes>();

        for (unsigned ip = 0; ip < n_integration_points; ++ip)
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

            double integral_measure = 1.0;
            if (_data.integral_measure)
            {
                integral_measure = (*_data.integral_measure)(t, position)[0];
            }

            // flux = alpha * ( u_0 - u )
            // adding a alpha term to the diagonal of the stiffness matrix
            // and a alpha * u_0 term to the rhs vector
            _local_K.diagonal().noalias() +=
                N * alpha.dot(N) * w * integral_measure;
            _local_rhs.noalias() +=
                N * alpha.dot(N) * u_0.dot(N) * w * integral_measure;
        }

        auto const indices = NumLib::getIndices(id, dof_table_boundary);
        K.add(NumLib::LocalToGlobalIndexMap::RowColumnIndices(indices, indices),
              _local_K);
        b.add(indices, _local_rhs);
    }

private:
    RobinBoundaryConditionData const& _data;

    typename Base::NodalMatrixType _local_K;
    typename Base::NodalVectorType _local_rhs;

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

}  // namespace ProcessLib
