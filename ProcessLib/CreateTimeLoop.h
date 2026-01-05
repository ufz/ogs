// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <map>
#include <memory>
#include <string>
#include <vector>

namespace BaseLib
{
class ConfigTree;
}

namespace MeshLib
{
class Mesh;
}

namespace NumLib
{
class NonlinearSolverBase;
}

namespace ProcessLib
{
class TimeLoop;
class Process;
}  // namespace ProcessLib

namespace ProcessLib
{
//! Builds a TimeLoop from the given configuration.
std::unique_ptr<TimeLoop> createTimeLoop(
    BaseLib::ConfigTree const& config, std::string const& output_directory,
    std::vector<std::unique_ptr<Process>> const& processes,
    std::map<std::string, std::unique_ptr<NumLib::NonlinearSolverBase>> const&
        nonlinear_solvers,
    std::vector<std::unique_ptr<MeshLib::Mesh>>& meshes,
    bool const compensate_non_equilibrium_initial_residuum);

}  // namespace ProcessLib
