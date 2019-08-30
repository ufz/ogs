/**
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
#include "PhreeqcKernelData/AqueousSolution.h"

#include "ThirdParty/iphreeqc/src/src/phreeqcpp/Phreeqc.h"

class cxxISolution;

namespace ChemistryLib
{
namespace PhreeqcKernelData
{

class PhreeqcKernel final : public ChemicalSolverInterface, private Phreeqc
{
public:
    PhreeqcKernel(std::size_t const num_chemical_systems,
                  std::vector<std::pair<int, std::string>> const&
                      process_id_to_component_name_map,
                  std::string const& database,
                  AqueousSolution& aqueous_solution,
                  cxxKinetics& kinetic_reactants,
                  std::tuple<rate*, int> const& reaction_rates);

    void doWaterChemistryCalculation(
        std::vector<GlobalVector*>& process_solutions,
        double const dt) override;

    void setAqueousSolutions(
        std::vector<GlobalVector*> const& process_solutions);

    void setTimeStep(double const dt) override;

    void execute(std::vector<GlobalVector*>& process_solutions);

    void updateNodalProcessSolutions(
        std::vector<GlobalVector*> const& process_solutions,
        std::size_t const node_id);

private:
    std::map<int, struct master*> _process_id_to_master_map;
    cxxISolution _templated_initial_aqueous_solution;
};
}  // namespace PhreeqcKernelData
}  // namespace ChemistryLib
