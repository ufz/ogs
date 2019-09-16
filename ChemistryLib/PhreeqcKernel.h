/**
 * \file
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <map>
#include <vector>

#include "ChemicalSolverInterface.h"
#include "PhreeqcKernelData/InitialAqueousSolution.h"
#include "PhreeqcKernelData/KineticReactant.h"

#include "ThirdParty/iphreeqc/src/src/phreeqcpp/Phreeqc.h"

class cxxSolution;

namespace ChemistryLib
{
namespace PhreeqcKernelData
{
class AqueousSolution;
class ReactionRate;

class PhreeqcKernel final : public ChemicalSolverInterface, private Phreeqc
{
public:
    PhreeqcKernel(std::size_t const num_chemical_systems,
                  std::vector<std::pair<int, std::string>> const&
                      process_id_to_component_name_map,
                  std::string const& database,
                  AqueousSolution& aqueous_solution,
                  std::unique_ptr<Kinetics>
                      kinetic_reactants,
                  std::vector<ReactionRate>&& reaction_rates);

    void doWaterChemistryCalculation(
        std::vector<GlobalVector*>& process_solutions,
        double const dt) override;

    void setAqueousSolutions(
        std::vector<GlobalVector*> const& process_solutions);

    void execute(std::vector<GlobalVector*>& process_solutions);

    void updateNodalProcessSolutions(
        std::vector<GlobalVector*> const& process_solutions,
        std::size_t const node_id);

private:
    void reinitializeRates();

    void setConvergenceTolerance() { convergence_tolerance = 1e-12; }

    void configureOutputSettings() { pr.all = false; }

    cxxISolution* getOrCreateInitialAqueousSolution(
        cxxSolution& aqueous_solution);

    bool isHydrogen(char const* element) const
    {
        return strcmp(element, "H") == 0 ? true : false;
    }

    void setTimeStepSize(double const dt);

    std::map<int, struct master*> _process_id_to_master_map;
    cxxISolution _initial_aqueous_solution;
    std::vector<ReactionRate> const _reaction_rates;
};
}  // namespace PhreeqcKernelData
}  // namespace ChemistryLib
