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

#include <memory>
#include <vector>

#include "MaterialLib/MPL/CreateMaterialSpatialDistributionMap.h"
namespace MeshLib
{
class Mesh;
}

namespace NumLib
{
class LocalToGlobalIndexMap;
}
namespace ParameterLib
{
struct ParameterBase;
}
namespace ProcessLib
{
class BoundaryCondition;
struct BoundaryConditionConfig;
class Process;
class ProcessVariable;

std::unique_ptr<BoundaryCondition> createBoundaryCondition(
    const BoundaryConditionConfig& config,
    const NumLib::LocalToGlobalIndexMap& dof_table,
    const MeshLib::Mesh& bulk_mesh, const int variable_id,
    const unsigned integration_order, const unsigned shapefunction_order,
    const std::vector<std::unique_ptr<ParameterLib::ParameterBase>>& parameters,
    const Process& process,
    std::vector<std::reference_wrapper<ProcessVariable>> const&
        all_process_variables_for_this_process,
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media);

}  // namespace ProcessLib
