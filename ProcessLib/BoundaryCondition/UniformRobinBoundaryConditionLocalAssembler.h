/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_UNIFORMROBINBOUNDARYCONDITIONLOCALASSEMBLER_H
#define PROCESSLIB_UNIFORMROBINBOUNDARYCONDITIONLOCALASSEMBLER_H

#include "NumLib/DOF/DOFTableUtil.h"

namespace ProcessLib
{
struct UniformRobinBoundaryConditionData final {
    double const alpha;
    double const f_0;
};

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
class UniformRobinBoundaryConditionLocalAssembler final
    : public GenericNaturalBoundaryConditionLocalAssembler<
          ShapeFunction, IntegrationMethod, GlobalDim>
{
    using Base = GenericNaturalBoundaryConditionLocalAssembler<
        ShapeFunction, IntegrationMethod, GlobalDim>;

public:
    UniformRobinBoundaryConditionLocalAssembler(
        MeshLib::Element const& e, std::size_t const local_matrix_size,
        unsigned const integration_order,
        UniformRobinBoundaryConditionData const& data)
        : Base(e, integration_order),
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
        (void)t;  // TODO time-dependent Robin BCs

        _local_K.setZero();
        _local_rhs.setZero();

        IntegrationMethod integration_method(Base::_integration_order);
        std::size_t const n_integration_points =
            integration_method.getNumberOfPoints();

        for (std::size_t ip = 0; ip < n_integration_points; ++ip) {
            auto const& sm = Base::_shape_matrices[ip];
            auto const& wp = integration_method.getWeightedPoint(ip);

            // df/dn = alpha ( f - f_0 )
            // adding a -alpha term to the diagonal of the stiffness matrix
            // and a -alpha * f_0 term to the rhs vector
            _local_K.diagonal().noalias() -=
                sm.N * _data.alpha * sm.detJ * wp.getWeight();
            _local_rhs.noalias() -=
                sm.N * _data.alpha * _data.f_0 * sm.detJ * wp.getWeight();
        }

        auto const indices = NumLib::getIndices(id, dof_table_boundary);
        K.add(NumLib::LocalToGlobalIndexMap::RowColumnIndices(indices, indices),
              _local_K);
        b.add(indices, _local_rhs);
    }

private:
    UniformRobinBoundaryConditionData const& _data;

    typename Base::NodalMatrixType _local_K;
    typename Base::NodalVectorType _local_rhs;
};

}  // ProcessLib

#endif  // PROCESSLIB_UNIFORMROBINBOUNDARYCONDITIONLOCALASSEMBLER_H
