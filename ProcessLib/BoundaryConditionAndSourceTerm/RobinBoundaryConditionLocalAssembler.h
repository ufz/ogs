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

#include "GenericNaturalBoundaryConditionLocalAssembler.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
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

template <typename ShapeFunction, int GlobalDim>
class RobinBoundaryConditionLocalAssembler final
    : public GenericNaturalBoundaryConditionLocalAssembler<ShapeFunction,
                                                           GlobalDim>
{
    using Base =
        GenericNaturalBoundaryConditionLocalAssembler<ShapeFunction, GlobalDim>;
    using ShapeMatricesType = ShapeMatrixPolicyType<ShapeFunction, GlobalDim>;

public:
    RobinBoundaryConditionLocalAssembler(
        MeshLib::Element const& e,
        std::size_t const local_matrix_size,
        NumLib::GenericIntegrationMethod const& integration_method,
        bool is_axially_symmetric,
        RobinBoundaryConditionData const& data)
        : Base(e, is_axially_symmetric, integration_method),
          _data(data),
          _local_K(local_matrix_size, local_matrix_size),
          _local_rhs(local_matrix_size)
    {
    }

    void assemble(std::size_t const id,
                  NumLib::LocalToGlobalIndexMap const& dof_table_boundary,
                  double const t, std::vector<GlobalVector*> const& xs,
                  int const process_id, GlobalMatrix& K, GlobalVector& b,
                  GlobalMatrix* Jac) override
    {
        _local_K.setZero();
        _local_rhs.setZero();

        auto& x = *xs[process_id];

        auto const indices = NumLib::getIndices(id, dof_table_boundary);
        auto const local_x = x.get(indices);
        auto const u =
            MathLib::toVector<Eigen::Matrix<double, ShapeFunction::NPOINTS, 1>>(
                local_x, ShapeFunction::NPOINTS);

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

            double const a = alpha.dot(N) * w * integral_measure;

            // The local K matrix is used for both, the Newton and Picard
            // methods.
            _local_K.noalias() += N.transpose() * N * a;

            if (Jac != nullptr)
            {
                _local_rhs.noalias() -= N.transpose() * (u - u_0).dot(N) * a;
            }
            else
            {
                _local_rhs.noalias() += N.transpose() * (u_0.dot(N) * a);
            }
        }

        b.add(indices, _local_rhs);
        if (Jac != nullptr)
        {
            Jac->add(NumLib::LocalToGlobalIndexMap::RowColumnIndices(indices,
                                                                     indices),
                     _local_K);
        }
        else
        {
            K.add(NumLib::LocalToGlobalIndexMap::RowColumnIndices(indices,
                                                                  indices),
                  _local_K);
        }
    }

private:
    RobinBoundaryConditionData const& _data;

    typename Base::NodalMatrixType _local_K;
    typename Base::NodalVectorType _local_rhs;

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

}  // namespace ProcessLib
