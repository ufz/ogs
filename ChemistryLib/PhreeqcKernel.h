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

#include <phreeqcpp/Phreeqc.h>

#include <map>
#include <vector>

#include "ChemicalSolverInterface.h"
#include "PhreeqcKernelData/EquilibriumReactants.h"
#include "PhreeqcKernelData/KineticReactant.h"

class cxxSolution;
class cxxISolution;

namespace ChemistryLib
{
namespace PhreeqcKernelData
{
class AqueousSolution;
class ReactionRate;

class PhreeqcKernel final : public ChemicalSolverInterface, private Phreeqc
{
public:
    PhreeqcKernel(GlobalLinearSolver& linear_solver,
                  std::size_t const num_chemical_systems,
                  std::vector<std::pair<int, std::string>> const&
                      process_id_to_component_name_map,
                  std::string const& database, AqueousSolution aqueous_solution,
                  std::unique_ptr<EquilibriumReactants>&& equilibrium_reactants,
                  std::unique_ptr<Kinetics>&& kinetic_reactants,
                  std::vector<ReactionRate>&& reaction_rates);

    void executeSpeciationCalculation(double const dt) override;

    void setAqueousSolutions(
        std::vector<GlobalVector*> const& process_solutions);

    void callPhreeqc(std::vector<GlobalVector*>& process_solutions);

    std::vector<GlobalVector*> getIntPtProcessSolutions() const override
    {
        return {};
    }

    void updateNodalProcessSolutions(
        std::vector<GlobalVector*> const& process_solutions,
        std::size_t const node_id);

private:
    void initializePhreeqcGeneralSettings() { do_initialize(); }

    void tidyEquilibriumReactants(
        EquilibriumReactants const equilibrium_reactants);

    void loadDatabase(std::string const& database);

    void reinitializeRates();

    void setConvergenceTolerance() { convergence_tolerance = 1e-12; }

    void configureOutputSettings() { pr.all = false; }

    cxxISolution* getOrCreateInitialAqueousSolution(
        cxxSolution& aqueous_solution);

    bool isHydrogen(char const* element) const
    {
        return strcmp(element, "H") == 0;
    }

    void setTimeStepSize(double const dt);

    void reset(std::size_t const chemical_system_id);

    std::map<int, struct master*> _process_id_to_master_map;
    std::unique_ptr<cxxISolution const> _initial_aqueous_solution;
    std::unique_ptr<cxxSolution const> _aqueous_solution;
    std::vector<ReactionRate> const _reaction_rates;
};
}  // namespace PhreeqcKernelData
}  // namespace ChemistryLib
