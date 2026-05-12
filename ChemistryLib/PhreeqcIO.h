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
#include "PhreeqcIOData/PhreeqcInstancePool.h"

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
 *  - In stream mode, executeSpeciationCalculationParallel() handles both
 *    serial (pool size 1) and parallel cases via the PhreeqcInstancePool.
 *  - In file mode, writeInputsToFile(), callPhreeqc(), readOutputsFromFile()
 *    implement PHREEQC I/O.
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
              bool use_stream_mode,
              int num_chemistry_threads,
              double negative_concentration_tolerance = 1e-12);

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

    /// Write the global input header (PHASES, KNOBS, SELECTED_OUTPUT,
    /// USER_PUNCH, RATES) shared by both file and stream mode.
    void writeInputHeader(std::ostream& os) const;

    /// Write the per-system input block (SOLUTION, EQUILIBRIUM_PHASES,
    /// KINETICS, SURFACE, EXCHANGE).
    /// @param os                Output stream to write the PHREEQC input to.
    /// @param chemical_system_id  Zero-based index of the chemical system
    ///                            (mesh node/cell) being written.
    /// @param solution_id       PHREEQC solution number for this system
    ///                          (1 in stream mode, chemical_system_id+1 in
    ///                          file mode).
    /// @param prev_solution_id  Solution number of the previous-timestep
    ///                          aqueous solution used for SURFACE/EXCHANGE
    ///                          equilibration.
    /// @param dt                Time step size in seconds.
    void writeSystemBlock(std::ostream& os,
                          std::size_t chemical_system_id,
                          std::size_t solution_id,
                          std::size_t prev_solution_id,
                          double dt) const;

    /// Parse one line of PHREEQC selected output and update the chemical
    /// system state.  Used by both output-reading code paths (file-mode
    /// operator>> and stream-mode parseOutputForSystem()).
    void updateSystemFromOutputLine(std::string_view line,
                                    std::size_t chemical_system_id);

    void setAqueousSolutionsPrevFromDumpString(std::string_view dump_content);

    void callPhreeqc() const;

    void readOutputsFromFile();

    void executeSpeciationCalculationParallel(double const dt);

    std::string generateInputForSystem(std::size_t chemical_system_id,
                                       double const dt) const;

    void parseOutputForSystem(std::string_view output_content,
                              std::size_t chemical_system_id);

    void updateChemicalSystemFromOutput(
        std::vector<double> const& accepted_items,
        std::size_t chemical_system_id);

    PhreeqcIO& operator<<(double const dt)
    {
        _dt = dt;
        return *this;
    }

    // Member variables are ordered by size (largest first) to minimize padding.
    std::string const _database;
    Knobs const _knobs;
    std::vector<ReactionRate> const _reaction_rates;
    std::unique_ptr<ChemicalSystem> _chemical_system;
    std::unique_ptr<UserPunch> _user_punch;
    std::unique_ptr<Output> const _output;
    std::unique_ptr<Dump> const _dump;
    std::unique_ptr<PhreeqcInstancePool> instance_pool_;
    double _dt = std::numeric_limits<double>::quiet_NaN();
    std::size_t _num_chemical_systems;
    /// Concentrations between this magnitude (negative) and zero are
    /// silently clamped; more-negative values trigger a per-call WARN.
    double _negative_concentration_tolerance;
    /// ID of the standalone file-mode PHREEQC instance (independent of the
    /// parallel pool). Created by PhreeqcInstancePool::createInstance() in the
    /// ctor only in file mode; in stream mode the pool is used instead and
    /// this stays -1. It must not be hard-coded to 0 because the IPhreeqc
    /// library hands out monotonically increasing IDs and another consumer
    /// (e.g. a unit test that constructed an IPhreeqc instance earlier) may
    /// have used 0 first.
    int phreeqc_instance_id = -1;
    int num_chemistry_threads_;
    bool _use_stream_mode;
};
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
