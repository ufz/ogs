// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

/**
 * \file
 * \brief Chemical-solver interface used in OGS operator-split reactive
 * transport.
 *
 * OpenGeoSys advances reactive transport by splitting it into two subproblems:
 *  (1) transport (advection / diffusion) in OpenGeoSys, and
 *  (2) local chemical reaction / speciation in an external geochemical solver.
 *
 * The domain is partitioned into independent chemical systems. Each chemical
 * system is identified by an integer \c chemical_system_id, typically mapped
 * from a reactive mesh node or control volume. During a chemistry step, each
 * chemical system is treated as a closed, well-mixed reactor. No mass is
 * exchanged between systems during this step.
 *
 * For each \c chemical_system_id, OpenGeoSys provides the external solver
 * (e.g. PHREEQC) with the current chemical state, executes the solver, and
 * reads back the reacted state (updated component totals, pH, mineral /
 * porosity effects). These updated values then serve as the chemical state
 * for the next transport solve in OpenGeoSys.
 *
 * This header declares the factory for creating a concrete chemistry backend
 * (e.g. a PHREEQC-based implementation).
 */
#pragma once

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "ChemicalSolverType.h"
#include "MathLib/LinAlg/GlobalLinearSolverType.h"

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
