// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

/**
 * \file
 * \brief PHREEQC-backed ChemicalSolverInterface implementation.
 *
 * PhreeqcIO is the concrete coupling layer between OpenGeoSys and PHREEQC.
 * For each \c chemical_system_id (one reactive control volume / local
 * chemical system), it:
 *
 *  - prepares PHREEQC input blocks (SOLUTION, EQUILIBRIUM_PHASES, KINETICS,
 *    EXCHANGE, SURFACE, etc.) from the transported component totals
 *    \f$c_{T\alpha}\f$, the porous-medium state (porosity, mineral fractions,
 *    reactive surface area), and the process variables (e.g. temperature,
 *    pressure);
 *
 *  - calls PHREEQC through the IPhreeqc API, advancing each
 *    \c chemical_system_id as a closed, well-mixed batch reactor over the
 *    current timestep \f$\Delta t\f$ (i.e. no mass exchange between different
 *    \c chemical_system_id inside PHREEQC during that solve);
 *
 *  - reads PHREEQC output (updated totals \f$c_{T\alpha}\f$, pH, pe/redox,
 *    mineral amounts, etc.) and writes these reacted values back into
 *    OpenGeoSys;
 *
 *  - updates porosity and solid volume fractions to reflect precipitation
 *    and dissolution, so that flow / transport / mechanics in the next
 *    timestep see the chemically modified medium.
 *
 * Temperature and pressure:
 *  The interface allows passing temperature \f$T\f$ [K] and pressure \f$p\f$
 *  [Pa] for each \c chemical_system_id via setChemicalSystemConcrete().
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
struct ChemicalSystem;
struct ReactionRate;
struct Output;
struct Dump;
struct UserPunch;

/**
 * \class PhreeqcIO
 * \brief Drives the chemistry step in the operator-split reactive transport
 *        loop.
 *
 * Per timestep:
 *
 *  1. initializeChemicalSystemConcrete(...) /
 *     setChemicalSystemConcrete(...):
 *        For each \c chemical_system_id, store the transported totals
 *        \f$c_{T\alpha}\f$, the current porous-medium state (porosity,
 *        mineral fractions, reactive surface area), and the process
 *        variables (e.g. temperature \f$T\f$, pressure \f$p\f$, time \f$t\f$,
 *        and timestep \f$\Delta t\f$) in the corresponding ChemicalSystem.
 *        Each \c chemical_system_id is treated as one local batch reactor.
 *
 *  2. executeSpeciationCalculation(dt):
 *        Build PHREEQC input for all \c chemical_system_id, run PHREEQC via
 *        IPhreeqc, and obtain the reacted state after advancing chemistry
 *        over \f$\Delta t\f$:
 *          - updated aqueous composition,
 *          - reaction/source terms \f$R_{\alpha}\f$ for each transported
 *            component,
 *          - updated mineral amounts.
 *
 *  3. updateVolumeFractionPostReaction(...),
 *     updatePorosityPostReaction(...):
 *        Apply precipitation / dissolution back into OpenGeoSys by updating
 *        mineral / solid volume fractions and porosity for each local system.
 *        These updates feed into flow / transport / mechanics in the next
 *        global step.
 *
 *  4. computeSecondaryVariable(...):
 *        Optionally expose derived quantities (e.g. pH fields, mineral
 *        fractions, surface loading) for output.
 *
 * Accessors:
 *  - getConcentration(component_id, chemical_system_id):
 *      Reacted total \f$c_{T\alpha}\f$ of a transported component for that
 *      local system after the last chemistry step. Used as the starting
 *      composition for the next transport solve.
 *
 *  - getComponentList():
 *      List of transported components in the order expected by the transport
 *      process.
 *
 * Internal data:
 *  - _chemical_system holds the per-\c chemical_system_id definition
 *    (aqueous solution, kinetic and equilibrium reactants, exchangers,
 *    surface sites).
 *  - writeInputsToFile(), callPhreeqc(), readOutputsFromFile() implement
 *    PHREEQC I/O.
 */
class PhreeqcIO final : public ChemicalSolverInterface
{
public:
    PhreeqcIO(MeshLib::Mesh const& mesh,
              GlobalLinearSolver& linear_solver,
              std::string const& project_file_name,
              std::string&& database,
              std::unique_ptr<ChemicalSystem>&& chemical_system,
              std::vector<ReactionRate>&& reaction_rates,
              std::unique_ptr<UserPunch>&& user_punch,
              std::unique_ptr<Output>&& output,
              std::unique_ptr<Dump>&& dump,
              Knobs&& knobs,
              bool use_stream_mode);

    ~PhreeqcIO();

    void initialize() override;

    void initializeChemicalSystemConcrete(
        std::vector<double> const& concentrations,
        GlobalIndexType const& chemical_system_id,
        MaterialPropertyLib::Medium const& medium,
        ParameterLib::SpatialPosition const& pos,
        double const t) override;

    void setChemicalSystemConcrete(
        std::vector<double> const& concentrations,
        GlobalIndexType const& chemical_system_id,
        MaterialPropertyLib::Medium const* medium,
        MaterialPropertyLib::VariableArray const& vars,
        ParameterLib::SpatialPosition const& pos, double const t,
        double const dt) override;

    void setAqueousSolutionsPrevFromDumpFile() override;

    void executeSpeciationCalculation(double const dt) override;

    double getConcentration(
        int const component_id,
        GlobalIndexType const chemical_system_id) const override;

    friend std::ostream& operator<<(std::ostream& os,
                                    PhreeqcIO const& phreeqc_io);

    friend std::istream& operator>>(std::istream& in, PhreeqcIO& phreeqc_io);

    void updateVolumeFractionPostReaction(
        GlobalIndexType const& chemical_system_id,
        MaterialPropertyLib::Medium const& medium,
        ParameterLib::SpatialPosition const& pos, double const porosity,
        double const t, double const dt) override;

    void updatePorosityPostReaction(GlobalIndexType const& chemical_system_id,
                                    MaterialPropertyLib::Medium const& medium,
                                    double& porosity) override;

    void computeSecondaryVariable(
        std::size_t const ele_id,
        std::vector<GlobalIndexType> const& chemical_system_indices) override;

    std::vector<std::string> const getComponentList() const override;

    std::string const _phreeqc_input_file;

private:
    void writeInputsToFile(double const dt);

    std::stringstream writeInputsToStringStream(double const dt);

    void setAqueousSolutionsPrevFromDumpString(std::string const& dump_content);

    void callPhreeqc() const;

    void callPhreeqcWithString(std::string const& input_content) const;

    std::string retrieveSelectedOutputString() const;

    std::string retrieveDumpString() const;

    void readOutputsFromFile();

    void readOutputsFromString(std::string const& output_content);

    PhreeqcIO& operator<<(double const dt)
    {
        _dt = dt;
        return *this;
    }

    std::string const _database;
    std::unique_ptr<ChemicalSystem> _chemical_system;
    std::vector<ReactionRate> const _reaction_rates;
    std::unique_ptr<UserPunch> _user_punch;
    std::unique_ptr<Output> const _output;
    std::unique_ptr<Dump> const _dump;
    Knobs const _knobs;
    double _dt = std::numeric_limits<double>::quiet_NaN();
    const int phreeqc_instance_id = 0;
    std::size_t _num_chemical_systems = -1;
    bool _use_stream_mode = false;
    mutable std::string _last_input_content;   // For diagnostics
    mutable std::string _last_output_content;  // For diagnostics
    mutable std::string _last_dump_content;    // For diagnostics
};
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
