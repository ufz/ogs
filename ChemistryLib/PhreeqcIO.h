/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>

#include "MathLib/LinAlg/GlobalMatrixVectorTypes.h"
#include "Output.h"
#include "PhreeqcIOData/AqueousSolution.h"
#include "PhreeqcIOData/EquilibriumPhase.h"
#include "PhreeqcIOData/KineticReactant.h"
#include "PhreeqcIOData/ReactionRate.h"

namespace ChemistryLib
{
enum class Status
{
    SettingAqueousSolutions,
    UpdatingProcessSolutions
};

class PhreeqcIO
{
public:
    PhreeqcIO(std::string const& project_file_name,
              std::string&& database,
              std::vector<AqueousSolution>&& aqueous_solutions,
              std::vector<std::vector<EquilibriumPhase>>&& equilibrium_phases,
              std::vector<std::vector<KineticReactant>>&& kinetic_reactants,
              std::vector<ReactionRate>&& reaction_rates,
              std::unique_ptr<Output>&& output,
              std::vector<std::pair<int, std::string>> const&
                  process_id_to_component_name_map)
        : _phreeqc_input_file(project_file_name + "_phreeqc.inp"),
          _database(std::move(database)),
          _aqueous_solutions(std::move(aqueous_solutions)),
          _equilibrium_phases(std::move(equilibrium_phases)),
          _kinetic_reactants(std::move(kinetic_reactants)),
          _reaction_rates(std::move(reaction_rates)),
          _output(std::move(output)),
          _process_id_to_component_name_map(process_id_to_component_name_map)
    {
    }

    void setAqueousSolutionsOrUpdateProcessSolutions(
        std::vector<GlobalVector*> const& process_solutions,
        Status const status);

    void setTimeStep(double const dt) { _dt = dt; }

    std::string const _phreeqc_input_file;

private:
    std::string const _database;
    std::vector<AqueousSolution> _aqueous_solutions;
    std::vector<std::vector<EquilibriumPhase>> _equilibrium_phases;
    std::vector<std::vector<KineticReactant>> _kinetic_reactants;
    std::vector<ReactionRate> const _reaction_rates;
    std::unique_ptr<Output> const _output;
    std::vector<std::pair<int, std::string>> const&
        _process_id_to_component_name_map;
    double _dt;
};
}  // namespace ChemistryLib
