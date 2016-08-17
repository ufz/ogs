/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_UNIFORMNEUMANNBOUNDARYCONDITIONLOCALASSEMBLER_H
#define PROCESSLIB_UNIFORMNEUMANNBOUNDARYCONDITIONLOCALASSEMBLER_H

#include "NumLib/DOF/DOFTableUtil.h"
#include "ProcessLib/Parameter/Parameter.h"
#include "GenericNaturalBoundaryConditionLocalAssembler.h"

namespace ProcessLib
{
template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
class UniformNeumannBoundaryConditionLocalAssembler final
    : public GenericNaturalBoundaryConditionLocalAssembler<
          ShapeFunction, IntegrationMethod, GlobalDim>
{
    using Base = GenericNaturalBoundaryConditionLocalAssembler<
        ShapeFunction, IntegrationMethod, GlobalDim>;

public:
    /// The neumann_bc_value factor is directly integrated into the local
    /// element matrix.
    UniformNeumannBoundaryConditionLocalAssembler(
        MeshLib::Element const& e,
        std::size_t const local_matrix_size,
        unsigned const integration_order,
        Parameter<double> const& neumann_bc_value)
        : Base(e, integration_order),
          _neumann_bc_value(neumann_bc_value),
          _local_rhs(local_matrix_size)
    {
    }

    void assemble(std::size_t const id,
                  NumLib::LocalToGlobalIndexMap const& dof_table_boundary,
                  double const t, const GlobalVector& /*x*/,
                  GlobalMatrix& /*K*/, GlobalVector& b) override
    {
        _local_rhs.setZero();

        IntegrationMethod integration_method(Base::_integration_order);
        std::size_t const n_integration_points =
            integration_method.getNumberOfPoints();

        SpatialPosition pos;
        pos.setElementID(id);

        for (std::size_t ip(0); ip < n_integration_points; ip++)
        {
            pos.setIntegrationPoint(ip);
            auto const& sm = Base::_shape_matrices[ip];
            auto const& wp = integration_method.getWeightedPoint(ip);
            _local_rhs.noalias() += sm.N *
                                    _neumann_bc_value.getTuple(t, pos).front() *
                                    sm.detJ * wp.getWeight();
        }

        auto const indices = NumLib::getIndices(id, dof_table_boundary);
        b.add(indices, _local_rhs);
    }

private:
    Parameter<double> const& _neumann_bc_value;
    typename Base::NodalVectorType _local_rhs;
};

}   // namespace ProcessLib

#endif  // PROCESSLIB_UNIFORMNEUMANNBOUNDARYCONDITIONLOCALASSEMBLER_H
