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

#include <vector>

#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "Output.h"  // TODO avoid that

namespace ProcessLib
{

/// Prepare the output data, i.e. add the solution to vtu data structure.
void addProcessDataToMesh(
    const double t, std::vector<GlobalVector*> const& x, int const process_id,
    MeshLib::Mesh& mesh,
    std::vector<NumLib::LocalToGlobalIndexMap const*> const& mesh_dof_tables,
    Process const& process, bool const output_secondary_variable,
    OutputDataSpecification const& process_output);

/// Prepare the output data, i.e. add the solution to vtu data structure.
void addProcessDataToMesh(const double t, std::vector<GlobalVector*> const& x,
                          int const process_id, MeshLib::Mesh& mesh,
                          Process const& process,
                          bool const output_secondary_variable,
                          OutputDataSpecification const& process_output);
}  // namespace ProcessLib
