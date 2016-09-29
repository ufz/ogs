/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "NumLib/NumericsConfig.h"

#include "ProcessLib/BoundaryCondition/BoundaryCondition.h"

#include "ProcessLib/SmallDeformationWithLIE/Common/FractureProperty.h"

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
namespace SmallDeformationWithLIE
{

class BoundaryConditionBuilder : public ProcessLib::BoundaryConditionBuilder
{
public:
    explicit BoundaryConditionBuilder(FractureProperty const& fracture_prop)
    : _fracture_prop(fracture_prop) {}

    std::unique_ptr<BoundaryCondition> createBoundaryCondition(
        const BoundaryConditionConfig& config,
        const NumLib::LocalToGlobalIndexMap& dof_table, const MeshLib::Mesh& mesh,
        const int variable_id, const unsigned integration_order,
        const std::vector<std::unique_ptr<ProcessLib::ParameterBase>>& parameters) override;

private:
    FractureProperty const& _fracture_prop;
};

}  // SmallDeformationWithLIE
}  // ProcessLib
