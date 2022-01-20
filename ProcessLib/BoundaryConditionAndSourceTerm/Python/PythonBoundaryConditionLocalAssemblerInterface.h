/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "ProcessLib/BoundaryConditionAndSourceTerm/GenericNaturalBoundaryCondition.h"

namespace ProcessLib
{

class PythonBoundaryConditionLocalAssemblerInterface
    : public GenericNaturalBoundaryConditionLocalAssemblerInterface
{
public:
    //! Interpolates the given component of the given variable to the given \c
    //! local_node_id.
    //!
    //! The \c local_node_id is the number of the node within the current
    //! boundary element.
    virtual double interpolate(
        unsigned const local_node_id,
        NumLib::LocalToGlobalIndexMap const& dof_table_boundary,
        GlobalVector const& x, int const var, int const comp) const = 0;
};

}  // namespace ProcessLib
