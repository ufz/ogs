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
#include "OutputDataSpecification.h"
#include "ProcessLib/Process.h"

namespace ProcessLib
{

/// Adds all data (e.g., primary & secondary variables, integration point data)
/// from the given \c process to the given \c mesh.
///
/// The given \c mesh can be the full simulation domain or a subset thereof.
///
/// \sa ProcessLib::addProcessDataToMesh()
void addProcessDataToSubMesh(
    const double t, std::vector<GlobalVector*> const& x, int const process_id,
    MeshLib::Mesh& mesh,
    std::vector<NumLib::LocalToGlobalIndexMap const*> const& mesh_dof_tables,
    Process const& process, bool const output_secondary_variable,
    OutputDataSpecification const& process_output);

/// Adds all data (e.g., primary & secondary variables, integration point data)
/// from the given \c process to the given \c mesh.
///
/// \note The given \c mesh must be the full simulation domain, i.e., the full
/// domain on which the given \c process has been solved.
void addProcessDataToMesh(const double t, std::vector<GlobalVector*> const& x,
                          int const process_id, MeshLib::Mesh& mesh,
                          Process const& process,
                          bool const output_secondary_variable,
                          OutputDataSpecification const& process_output);
}  // namespace ProcessLib
