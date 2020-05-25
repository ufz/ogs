/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>

#include "ChemicalSolverInterface.h"
#include "PhreeqcIOData/Knobs.h"

namespace MeshLib
{
class Mesh;
}

namespace ChemistryLib
{
namespace PhreeqcIOData
{
struct AqueousSolution;
struct EquilibriumReactant;
struct KineticReactant;
struct ReactionRate;
struct Output;
struct SurfaceSite;
struct Dump;
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
              MeshLib::Mesh const& mesh,
              std::string&& database,
              std::vector<AqueousSolution>&& aqueous_solutions,
              std::vector<EquilibriumReactant>&& equilibrium_reactants,
              std::vector<KineticReactant>&& kinetic_reactants,
              std::vector<ReactionRate>&& reaction_rates,
              std::vector<SurfaceSite>&& surface,
              std::unique_ptr<UserPunch>&& user_punch,
              std::unique_ptr<Output>&& output,
              std::unique_ptr<Dump>&& dump,
              Knobs&& knobs);

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

    std::vector<std::string> const getComponentList() const override;

    std::string const phreeqc_input_file_;

private:
    PhreeqcIO& operator<<(double const dt)
    {
        dt_ = dt;
        return *this;
    }

    void setAqueousSolutionsPrevFromDumpFile();

    MeshLib::Mesh const& mesh_;
    std::string const database_;
    std::vector<AqueousSolution> aqueous_solutions_;
    std::vector<EquilibriumReactant> equilibrium_reactants_;
    std::vector<KineticReactant> kinetic_reactants_;
    std::vector<ReactionRate> const reaction_rates_;
    std::vector<SurfaceSite> const surface_;
    std::unique_ptr<UserPunch> user_punch_;
    std::unique_ptr<Output> const output_;
    std::unique_ptr<Dump> const dump_;
    Knobs const knobs_;
    double dt_ = std::numeric_limits<double>::quiet_NaN();
    const int phreeqc_instance_id = 0;
};
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
