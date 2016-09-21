/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_ROBINBOUNDARYCONDITIONLOCALASSEMBLER_H
#define PROCESSLIB_ROBINBOUNDARYCONDITIONLOCALASSEMBLER_H

#include "NumLib/DOF/DOFTableUtil.h"
#include "ProcessLib/Parameter/Parameter.h"
#include "GenericNaturalBoundaryConditionLocalAssembler.h"

namespace ProcessLib
{
struct RobinBoundaryConditionData final {
    Parameter<double> const& alpha;
    Parameter<double> const& u_0;
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
          _data(data),
          _local_K(local_matrix_size, local_matrix_size),
          _local_rhs(local_matrix_size)
    {
    }

    // TODO also implement derivative for Jacobian in Newton scheme.
    void assemble(std::size_t const id,
                  NumLib::LocalToGlobalIndexMap const& dof_table_boundary,
                  double const t, const GlobalVector& /*x*/, GlobalMatrix& K,
                  GlobalVector& b) override
    {
        _local_K.setZero();
        _local_rhs.setZero();

        unsigned const n_integration_points =
            Base::_integration_method.getNumberOfPoints();

        SpatialPosition pos;
        pos.setElementID(id);

        for (unsigned ip = 0; ip < n_integration_points; ++ip)
        {
            pos.setIntegrationPoint(ip);
            auto const& sm = Base::_shape_matrices[ip];
            auto const& wp = Base::_integration_method.getWeightedPoint(ip);

            double const alpha = _data.alpha(t, pos)[0];
            double const u_0 = _data.u_0(t, pos)[0];

            // flux = alpha * ( u_0 - u )
            // adding a alpha term to the diagonal of the stiffness matrix
            // and a alpha * u_0 term to the rhs vector
            _local_K.diagonal().noalias() +=
                sm.N * alpha * sm.detJ * wp.getWeight() * sm.integralMeasure;
            _local_rhs.noalias() += sm.N * alpha * u_0 * sm.detJ *
                                    wp.getWeight() * sm.integralMeasure;
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
};

}  // ProcessLib

#endif  // PROCESSLIB_ROBINBOUNDARYCONDITIONLOCALASSEMBLER_H
