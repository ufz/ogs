/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_UTILS_VECTORUTIL_H
#define PROCESSLIB_UTILS_VECTORUTIL_H

#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "NumLib/NumericsConfig.h"

namespace ProcessLib
{
//! Returns the value for the given \c node_id and \c global_component_id from
//! the vector \c x.
double getNodalValue(GlobalVector const& x, MeshLib::Mesh const& mesh,
                     NumLib::LocalToGlobalIndexMap const& dof_table,
                     std::size_t const node_id,
                     std::size_t const global_component_id);

} // namespace ProcessLib


#endif // PROCESSLIB_UTILS_VECTORUTIL_H
