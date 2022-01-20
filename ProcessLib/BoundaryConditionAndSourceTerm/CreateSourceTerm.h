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

#include <memory>
#include <vector>

namespace MeshLib
{
class Mesh;
}  // namespace MeshLib

namespace NumLib
{
class LocalToGlobalIndexMap;
}  // namespace NumLib

namespace ParameterLib
{
struct ParameterBase;
}  // namespace ParameterLib

namespace ProcessLib
{
class ProcessVariable;
class SourceTerm;
struct SourceTermConfig;

std::unique_ptr<SourceTerm> createSourceTerm(
    const SourceTermConfig& config,
    const NumLib::LocalToGlobalIndexMap& dof_table_bulk, const MeshLib::Mesh& source_term_mesh,
    const int variable_id, const unsigned integration_order,
    const unsigned shapefunction_order,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    std::vector<std::reference_wrapper<ProcessVariable>> const&
        all_process_variables_for_this_process);

}  // namespace ProcessLib
