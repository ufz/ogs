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

#include <memory>

#include "ChemicalSolverInterface.h"

namespace ChemistryLib
{
namespace PhreeqcIOData
{
struct AqueousSolution;
struct EquilibriumPhase;
struct KineticReactant;
struct ReactionRate;
struct Output;
struct SurfaceSite;
struct Dump;
struct Knobs;
struct UserPunch;

enum class Status
{
    SettingAqueousSolutions,
    UpdatingProcessSolutions
};

class PhreeqcIO final : public ChemicalSolverInterface
{
public:
    PhreeqcIO(std::string const project_file_name,
              std::string&& database,
              std::vector<AqueousSolution>&& aqueous_solutions,
              std::vector<EquilibriumPhase>&& equilibrium_phases,
              std::vector<KineticReactant>&& kinetic_reactants,
              std::vector<ReactionRate>&& reaction_rates,
              std::vector<SurfaceSite>&& surface,
              std::unique_ptr<UserPunch>&& user_punch,
              std::unique_ptr<Output>&& output,
              std::unique_ptr<Dump>&& dump,
              std::unique_ptr<Knobs>&& knobs,
              std::vector<std::pair<int, std::string>> const&
                  process_id_to_component_name_map);

    void executeInitialCalculation(
        std::vector<GlobalVector*>& process_solutions) override;

    void doWaterChemistryCalculation(
        std::vector<GlobalVector*>& process_solutions,
        double const dt) override;

    void setAqueousSolutionsOrUpdateProcessSolutions(
        std::vector<GlobalVector*> const& process_solutions,
        Status const status);

    void writeInputsToFile(double const dt = 0);

    void execute();

    void readOutputsFromFile();

    friend std::ostream& operator<<(std::ostream& os,
                                    PhreeqcIO const& phreeqc_io);

    friend std::istream& operator>>(std::istream& in, PhreeqcIO& phreeqc_io);

    std::string const _phreeqc_input_file;

private:
    PhreeqcIO& operator<<(double const dt)
    {
        _dt = dt;
        return *this;
    }

    std::string const _database;
    std::vector<AqueousSolution> _aqueous_solutions;
    std::vector<EquilibriumPhase> _equilibrium_phases;
    std::vector<KineticReactant> _kinetic_reactants;
    std::vector<ReactionRate> const _reaction_rates;
    std::vector<SurfaceSite> const _surface;
    std::unique_ptr<UserPunch> _user_punch;
    std::unique_ptr<Output> const _output;
    std::unique_ptr<Dump> const _dump;
    std::unique_ptr<Knobs> const _knobs;
    std::vector<std::pair<int, std::string>> const&
        _process_id_to_component_name_map;
    double _dt = std::numeric_limits<double>::quiet_NaN();
    const int phreeqc_instance_id = 0;
};
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
