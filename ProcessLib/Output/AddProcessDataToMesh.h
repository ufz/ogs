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

#include "OutputDataSpecification.h"
#include "ProcessOutputData.h"

namespace ProcessLib
{
/// Adds output data to a mesh.
///
/// \note The mesh is passed via \c process_output_data.
void addProcessDataToMesh(const double t, std::vector<GlobalVector*> const& xs,
                          int const process_id,
                          ProcessOutputData const& process_output_data,
                          bool const output_secondary_variables,
                          OutputDataSpecification const& process_output);
}  // namespace ProcessLib
