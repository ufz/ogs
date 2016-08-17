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
template <typename>
struct IndexValueVector;
}

namespace ProcessLib
{
struct BoundaryConditionConfig;
struct ParameterBase;

class BoundaryCondition
{
public:
    //! Applies natural BCs (i.e. non-Dirichlet BCs) to the stiffness matrix
    //! \c K and the vector \c b.
    virtual void applyNaturalBC(const double /*t*/, GlobalVector const& /*x*/,
                                GlobalMatrix& /*K*/, GlobalVector& /*b*/)
    {
        // By default it is assumed that the BC is not a natural BC. Therefore
        // there is nothing to do here.
    }

    //! Writes the values of Dirichlet BCs to \c bc_values.
    virtual void getDirichletBCValues(
        const double /*t*/,
        NumLib::IndexValueVector<GlobalIndexType>& /*bc_values*/) const
    {
        // By default it is assumed that the BC is not a Dirichlet BC. Therefore
        // there is nothing to do here.
    }

    virtual void preTimestep(const double /*t*/) {}

    virtual ~BoundaryCondition() = default;
};

std::unique_ptr<BoundaryCondition> createBoundaryCondition(
    const BoundaryConditionConfig& config,
    const NumLib::LocalToGlobalIndexMap& dof_table, const MeshLib::Mesh& mesh,
    const int variable_id, const unsigned integration_order,
    const std::vector<std::unique_ptr<ProcessLib::ParameterBase>>& parameters);

}  // ProcessLib

#endif  // PROCESSLIB_BOUNDARYCONDITION_H
