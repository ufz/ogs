/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>
#include "ProcessLib/Process.h"

namespace ProcessLib
{
std::unique_ptr<NodalSourceTerm> createNodalSourceTerm(
    BaseLib::ConfigTree const& config,
    const NumLib::LocalToGlobalIndexMap& dof_table, std::size_t mesh_id,
    std::size_t const node_id, const int variable_id, const int component_id);

}   // namespace ProcessLib
