/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

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
}  // namespace NumLib
namespace ParameterLib
{
struct ParameterBase;
}
namespace ProcessLib
{
struct BoundaryConditionConfig;
class Process;

class BoundaryCondition
{
public:
    //! Applies natural BCs (i.e. non-Dirichlet BCs) to the stiffness matrix
    //! \c K and the vector \c b.
    virtual void applyNaturalBC(const double /*t*/,
                                std::vector<GlobalVector*> const& /*x*/,
                                int const /*process_id*/, GlobalMatrix* /*K*/,
                                GlobalVector& /*b*/, GlobalMatrix* /*Jac*/)
    {
        // By default it is assumed that the BC is not a natural BC. Therefore
        // there is nothing to do here.
    }

    //! Writes the values of essential BCs to \c bc_values.
    virtual void getEssentialBCValues(
        const double /*t*/, GlobalVector const& /*x*/,
        NumLib::IndexValueVector<GlobalIndexType>& /*bc_values*/) const
    {
        // By default it is assumed that the BC is not an essential BC.
        // Therefore there is nothing to do here.
    }

    virtual void preTimestep(const double /*t*/,
                             std::vector<GlobalVector*> const& /*x*/,
                             int const /*process_id*/)
    {
        // A hook added for solution dependent dirichlet
    }

    virtual void postTimestep(const double /*t*/,
                              std::vector<GlobalVector*> const& /*x*/,
                              int const /*process_id*/)
    {
        // A hook added for solution dependent dirichlet
    }

    virtual ~BoundaryCondition() = default;
};

}  // namespace ProcessLib
