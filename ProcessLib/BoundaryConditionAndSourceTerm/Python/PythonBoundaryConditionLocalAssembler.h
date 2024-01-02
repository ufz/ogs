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

#include "NumLib/Fem/CoordinatesMapping/NaturalNodeCoordinates.h"
#include "PythonBoundaryCondition.h"
#include "Utils/BcAndStLocalAssemblerImpl.h"

namespace ProcessLib
{
template <typename ShapeFunction, typename LowerOrderShapeFunction,
          int GlobalDim>
class PythonBoundaryConditionLocalAssembler final
    : public PythonBoundaryConditionLocalAssemblerInterface
{
    using LocAsmImpl = ProcessLib::BoundaryConditionAndSourceTerm::Python::
        BcAndStLocalAssemblerImpl<PythonBcData, ShapeFunction,
                                  LowerOrderShapeFunction, GlobalDim>;
    using Traits = typename LocAsmImpl::Traits;

public:
    PythonBoundaryConditionLocalAssembler(
        MeshLib::Element const& e,
        std::size_t const /*local_matrix_size*/,
        NumLib::GenericIntegrationMethod const& integration_method,
        bool is_axially_symmetric,
        PythonBcData const& data)
        : impl_{e, integration_method, is_axially_symmetric, data}
    {
    }

    void assemble(std::size_t const boundary_element_id,
                  NumLib::LocalToGlobalIndexMap const& dof_table_boundary,
                  double const t, std::vector<GlobalVector*> const& xs,
                  int const process_id, GlobalMatrix& /*K*/, GlobalVector& b,
                  GlobalMatrix* const Jac) override
    {
        impl_.assemble(boundary_element_id, dof_table_boundary, t,
                       *xs[process_id], b, Jac);
    }

    double interpolate(unsigned const local_node_id,
                       NumLib::LocalToGlobalIndexMap const& dof_table_boundary,
                       GlobalVector const& x, int const var,
                       int const comp) const override
    {
        if constexpr (ShapeFunction::ORDER < 2 ||
                      LowerOrderShapeFunction::ORDER > 1)
        {
            return std::numeric_limits<double>::quiet_NaN();
        }
        else
        {
            auto const N = computeLowerOrderShapeMatrix(local_node_id);

            auto const nodal_values_base_node =
                ProcessLib::BoundaryConditionAndSourceTerm::Python::
                    collectDofsToMatrixOnBaseNodesSingleComponent(
                        impl_.element,
                        impl_.bc_or_st_data.bc_or_st_mesh.getID(),
                        dof_table_boundary, x, var, comp);

            return N * nodal_values_base_node;
        }
    }

private:
    typename Traits::LowerOrderShapeMatrix computeLowerOrderShapeMatrix(
        unsigned const local_node_id) const
    {
        using HigherOrderMeshElement = typename ShapeFunction::MeshElement;

        assert(local_node_id < impl_.element.getNumberOfNodes());

        std::array natural_coordss{MathLib::Point3d{NumLib::NaturalCoordinates<
            HigherOrderMeshElement>::coordinates[local_node_id]}};

        bool const is_axially_symmetric = false;  // does not matter for N
        auto const shape_matrices = NumLib::computeShapeMatrices<
            LowerOrderShapeFunction,
            typename Traits::LowerOrderShapeMatrixPolicy, GlobalDim,
            NumLib::ShapeMatrixType::N>(impl_.element, is_axially_symmetric,
                                        natural_coordss);
        return shape_matrices.front().N;
    }

    LocAsmImpl const impl_;
};

}  // namespace ProcessLib
