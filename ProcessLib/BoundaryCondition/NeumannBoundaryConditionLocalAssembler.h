/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_NEUMANNBOUNDARYCONDITIONLOCALASSEMBLER_H
#define PROCESSLIB_NEUMANNBOUNDARYCONDITIONLOCALASSEMBLER_H

#include "NumLib/DOF/DOFTableUtil.h"
#include "ProcessLib/Parameter/Parameter.h"
#include "GenericNaturalBoundaryConditionLocalAssembler.h"

namespace ProcessLib
{
template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
class NeumannBoundaryConditionLocalAssembler final
    : public GenericNaturalBoundaryConditionLocalAssembler<
          ShapeFunction, IntegrationMethod, GlobalDim>
{
    using Base = GenericNaturalBoundaryConditionLocalAssembler<
        ShapeFunction, IntegrationMethod, GlobalDim>;

public:
    /// The neumann_bc_value factor is directly integrated into the local
    /// element matrix.
    NeumannBoundaryConditionLocalAssembler(
        MeshLib::Element const& e,
        std::size_t const local_matrix_size,
        bool is_axially_symmetric,
        unsigned const integration_order,
        Parameter<double> const& neumann_bc_parameter)
        : Base(e, is_axially_symmetric, integration_order),
          _neumann_bc_parameter(neumann_bc_parameter),
          _local_rhs(local_matrix_size)
    {
    }

    void assemble(std::size_t const id,
                  NumLib::LocalToGlobalIndexMap const& dof_table_boundary,
                  double const t, const GlobalVector& /*x*/,
                  GlobalMatrix& /*K*/, GlobalVector& b) override
    {
        _local_rhs.setZero();

        unsigned const n_integration_points =
            Base::_integration_method.getNumberOfPoints();

        SpatialPosition pos;
        pos.setElementID(id);

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            pos.setIntegrationPoint(ip);
            auto const& sm = Base::_shape_matrices[ip];
            auto const& wp = Base::_integration_method.getWeightedPoint(ip);
            _local_rhs.noalias() += sm.N * _neumann_bc_parameter(t, pos)[0] *
                                    sm.detJ * wp.getWeight() *
                                    sm.integralMeasure;
        }

        auto const indices = NumLib::getIndices(id, dof_table_boundary);
        b.add(indices, _local_rhs);
    }

private:
    Parameter<double> const& _neumann_bc_parameter;
    typename Base::NodalVectorType _local_rhs;
};

}   // namespace ProcessLib

#endif  // PROCESSLIB_NEUMANNBOUNDARYCONDITIONLOCALASSEMBLER_H
