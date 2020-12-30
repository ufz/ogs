/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

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
}

namespace ProcessLib
{
//! Builds a TimeLoop from the given configuration.
std::unique_ptr<TimeLoop> createTimeLoop(
    BaseLib::ConfigTree const& config, std::string const& output_directory,
    std::vector<std::unique_ptr<Process>> const& processes,
    std::map<std::string, std::unique_ptr<NumLib::NonlinearSolverBase>> const&
        nonlinear_solvers,
    std::vector<std::unique_ptr<MeshLib::Mesh>> const& meshes);

}  // namespace ProcessLib
