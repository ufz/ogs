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
//! Returns nodal indices for the item identified by \c mesh_item_id from the
//! given \c dof_table.
std::vector<GlobalIndexType> getIndices(
    std::size_t const mesh_item_id,
    NumLib::LocalToGlobalIndexMap const& dof_table);
}  // namespace NumLib

#endif  // NUMLIB_DOF_DOFTABLEUTIL_H
