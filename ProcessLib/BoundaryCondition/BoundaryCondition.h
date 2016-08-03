/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_BOUNDARYCONDITION_H
#define PROCESSLIB_BOUNDARYCONDITION_H

#include "NumLib/NumericsConfig.h"

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
struct BoundaryConditionConfig;

class BoundaryCondition
{
public:
    virtual void apply(const double t, GlobalVector const& x, GlobalMatrix& K,
                       GlobalVector& b) = 0;
    virtual ~BoundaryCondition() = default;
};

std::unique_ptr<BoundaryCondition> createBoundaryCondition(
    const BoundaryConditionConfig& config,
    const NumLib::LocalToGlobalIndexMap& dof_table,
    const MeshLib::Mesh& mesh,
    const int variable_id,
    const unsigned integration_order);

}  // ProcessLib

#endif  // PROCESSLIB_BOUNDARYCONDITION_H
