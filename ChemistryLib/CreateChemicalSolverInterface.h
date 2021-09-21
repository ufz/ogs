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

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "ChemicalSolverType.h"
#include "MathLib/LinAlg/GlobalMatrixVectorTypes.h"

namespace BaseLib
{
class ConfigTree;
}

namespace MeshLib
{
class Mesh;
}

namespace ChemistryLib
{
class ChemicalSolverInterface;

template <ChemicalSolver chemical_solver>
std::unique_ptr<ChemicalSolverInterface> createChemicalSolverInterface(
    std::vector<std::unique_ptr<MeshLib::Mesh>> const& meshes,
    std::map<std::string, std::unique_ptr<GlobalLinearSolver>> const&
        linear_solvers,
    BaseLib::ConfigTree const& config, std::string const& output_directory);
}  // namespace ChemistryLib
