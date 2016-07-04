/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef NUMLIB_DOF_DOFTABLEUTIL_H
#define NUMLIB_DOF_DOFTABLEUTIL_H

#include "NumLib/DOF/LocalToGlobalIndexMap.h"

namespace NumLib
{
//! Returns nodal indices for the item identified by \c id from the given
//!  \c dof_table.
std::vector<GlobalIndexType> getIndices(
    std::size_t const id, NumLib::LocalToGlobalIndexMap const& dof_table);

//! Returns the values for the given \c dof_indices from the vector \c x.
std::vector<double> getLocalNodalDOFs(
    GlobalVector const& x, std::vector<GlobalIndexType> const& dof_indices);
}  // namespace NumLib

#endif  // NUMLIB_DOF_DOFTABLEUTIL_H
