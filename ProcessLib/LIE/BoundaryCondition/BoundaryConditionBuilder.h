/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "ProcessLib/BoundaryCondition/BoundaryCondition.h"

namespace MeshLib
{
class Mesh;
}

namespace NumLib
{
class LocalToGlobalIndexMap;
}

namespace ProcessLib
{
namespace LIE
{
struct FractureProperty;
}
}  // namespace ProcessLib

namespace ProcessLib
{
namespace LIE
{
/// A boundary condition builder for displacement jumps. Boundary
/// integration, e.g. for Neumann BC, should take into account the leveset
/// function.
class BoundaryConditionBuilder : public ProcessLib::BoundaryConditionBuilder
{
public:
    explicit BoundaryConditionBuilder(FractureProperty const& fracture_prop)
        : _fracture_prop(fracture_prop)
    {
    }

private:
    std::unique_ptr<BoundaryCondition> createNeumannBoundaryCondition(
        const BoundaryConditionConfig& config,
        const NumLib::LocalToGlobalIndexMap& dof_table,
        const MeshLib::Mesh& mesh, const int variable_id,
        const unsigned integration_order, const unsigned shapefunction_order,
        const std::vector<std::unique_ptr<ProcessLib::ParameterBase>>&
            parameters) override;

    FractureProperty const& _fracture_prop;
};

}  // namespace LIE
}  // namespace ProcessLib
