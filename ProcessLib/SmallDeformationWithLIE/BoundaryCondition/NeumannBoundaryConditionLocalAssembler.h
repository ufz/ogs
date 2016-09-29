/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "NumLib/DOF/DOFTableUtil.h"

#include "ProcessLib/Parameter/Parameter.h"
#include "ProcessLib/BoundaryCondition/GenericNaturalBoundaryConditionLocalAssembler.h"

#include "ProcessLib/SmallDeformationWithLIE/Common/LevelSetFunction.h"
#include "ProcessLib/SmallDeformationWithLIE/Common/Utils.h"

namespace ProcessLib
{
namespace SmallDeformationWithLIE
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
        Parameter<double> const& neumann_bc_parameter,
        FractureProperty const& fracture_prop,
        int variable_id)
        : Base(e, is_axially_symmetric, integration_order),
          _neumann_bc_parameter(neumann_bc_parameter),
          _local_rhs(local_matrix_size),
          _element(e),
          _fracture_prop(fracture_prop),
          _variable_id(variable_id)
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

            // levelset functions
            auto const ip_physical_coords = computePhysicalCoordinates(_element, sm.N);
            double const levelsets = calculateLevelSetFunction(_fracture_prop, ip_physical_coords.getCoords());

            _local_rhs.noalias() += sm.N * levelsets * _neumann_bc_parameter(t, pos)[0] *
                                    sm.detJ * wp.getWeight() *
                                    sm.integralMeasure;
        }

        auto const indices = NumLib::getIndices(id, dof_table_boundary);
        b.add(indices, _local_rhs);
    }

private:
    Parameter<double> const& _neumann_bc_parameter;
    typename Base::NodalVectorType _local_rhs;
    MeshLib::Element const& _element;
    FractureProperty const& _fracture_prop;
    int const _variable_id;
};

}   // namespace SmallDeformationWithLIE
}   // namespace ProcessLib

