/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_BOUNDARY_CONDITION_H_
#define PROCESS_LIB_BOUNDARY_CONDITION_H_

#include "BaseLib/ConfigTree.h"
#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "NumLib/NumericsConfig.h" // for GlobalIndexType

#include "DirichletBoundaryCondition.h"

namespace GeoLib
{
class GeoObject;
}

namespace MeshGeoToolsLib
{
class MeshNodeSearcher;
}

namespace ProcessLib
{
struct BoundaryConditionConfig;

/// The UniformDirichletBoundaryCondition class describes a constant in space
/// and time Dirichlet boundary condition.
/// The expected parameter in the passed configuration is "value" which, when
/// not present defaults to zero.
class UniformDirichletBoundaryCondition : public DirichletBoundaryCondition
{
public:
    UniformDirichletBoundaryCondition(
        NumLib::IndexValueVector<GlobalIndexType>&& bc)
        : _bc(std::move(bc))
    {
    }

    NumLib::IndexValueVector<GlobalIndexType> getBCValues()
    {
        return std::move(_bc);
    }

private:
    NumLib::IndexValueVector<GlobalIndexType> _bc;
};

std::unique_ptr<UniformDirichletBoundaryCondition>
createUniformDirichletBoundaryCondition(
    BoundaryConditionConfig const& config,
    MeshGeoToolsLib::MeshNodeSearcher& searcher,
    NumLib::LocalToGlobalIndexMap const& dof_table,
    int const variable_id);

}  // namespace ProcessLib

#endif  // PROCESS_LIB_BOUNDARY_CONDITION_H_
