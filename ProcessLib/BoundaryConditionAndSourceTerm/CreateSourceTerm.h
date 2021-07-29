/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>
#include <vector>

#include "ParameterLib/Parameter.h"

namespace MeshLib
{
class Mesh;
}

namespace NumLib
{
class LocalToGlobalIndexMap;
}  // namespace NumLib

namespace ProcessLib
{
class SourceTerm;
class NodalSourceTerm;
struct SourceTermConfig;
}  // namespace ProcessLib

namespace ProcessLib
{
std::unique_ptr<SourceTerm> createSourceTerm(
    const SourceTermConfig& config,
    const NumLib::LocalToGlobalIndexMap& dof_table, const MeshLib::Mesh& mesh,
    const int variable_id, const unsigned integration_order,
    const unsigned shapefunction_order,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
        parameters);

}  // namespace ProcessLib
