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

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
class RobinBoundaryConditionLocalAssembler final
    : public GenericNaturalBoundaryConditionLocalAssembler<
          ShapeFunction, IntegrationMethod, GlobalDim>
{
    using Base = GenericNaturalBoundaryConditionLocalAssembler<
        ShapeFunction, IntegrationMethod, GlobalDim>;

public:
    RobinBoundaryConditionLocalAssembler(MeshLib::Element const& e,
                                         std::size_t const local_matrix_size,
                                         bool is_axially_symmetric,
                                         unsigned const integration_order,
                                         RobinBoundaryConditionData const& data)
        : Base(e, is_axially_symmetric, integration_order),
          data_(data),
          local_K_(local_matrix_size, local_matrix_size),
          local_rhs_(local_matrix_size)
    {
    }

    // TODO also implement derivative for Jacobian in Newton scheme.
    void assemble(std::size_t const id,
                  NumLib::LocalToGlobalIndexMap const& dof_table_boundary,
                  double const t, std::vector<GlobalVector*> const& /*x*/,
                  int const /*process_id*/, GlobalMatrix& K, GlobalVector& b,
                  GlobalMatrix* /*Jac*/) override
    {
        local_K_.setZero();
        local_rhs_.setZero();

        unsigned const n_integration_points =
            Base::integration_method_.getNumberOfPoints();

        typename Base::NodalVectorType const alpha =
            data_.alpha.getNodalValuesOnElement(Base::element_, t)
                .template topRows<ShapeFunction::MeshElement::n_all_nodes>();
        typename Base::NodalVectorType const u_0 =
            data_.u_0.getNodalValuesOnElement(Base::element_, t)
                .template topRows<ShapeFunction::MeshElement::n_all_nodes>();

        ParameterLib::SpatialPosition position;
        position.setElementID(Base::element_.getID());

        for (unsigned ip = 0; ip < n_integration_points; ++ip)
        {
            position.setIntegrationPoint(ip);
            auto const& ip_data = Base::ns_and_weights_[ip];
            auto const& N = ip_data.N;
            auto const& w = ip_data.weight;

            double integral_measure = 1.0;
            if (data_.integral_measure)
            {
                integral_measure = (*data_.integral_measure)(t, position)[0];
            }

            // flux = alpha * ( u_0 - u )
            // adding a alpha term to the diagonal of the stiffness matrix
            // and a alpha * u_0 term to the rhs vector
            local_K_.diagonal().noalias() +=
                N * alpha.dot(N) * w * integral_measure;
            local_rhs_.noalias() +=
                N * alpha.dot(N) * u_0.dot(N) * w * integral_measure;
        }

        auto const indices = NumLib::getIndices(id, dof_table_boundary);
        K.add(NumLib::LocalToGlobalIndexMap::RowColumnIndices(indices, indices),
              local_K_);
        b.add(indices, local_rhs_);
    }

private:
    RobinBoundaryConditionData const& data_;

    typename Base::NodalMatrixType local_K_;
    typename Base::NodalVectorType local_rhs_;

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

}  // namespace ProcessLib
