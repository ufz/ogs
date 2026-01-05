// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "OutputDataSpecification.h"
#include "ProcessOutputData.h"

namespace ProcessLib
{
/// Adds output data to a mesh.
///
/// \note The mesh is passed via \c process_output_data.
void addProcessDataToMesh(NumLib::Time const& t,
                          std::vector<GlobalVector*> const& xs,
                          int const process_id,
                          ProcessOutputData const& process_output_data,
                          bool const output_secondary_variables,
                          OutputDataSpecification const& process_output);
}  // namespace ProcessLib
