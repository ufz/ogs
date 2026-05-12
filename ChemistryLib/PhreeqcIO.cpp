// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "PhreeqcIO.h"

#include <IPhreeqc.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include <algorithm>
#include <boost/algorithm/string.hpp>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <mutex>
#include <numeric>
#include <sstream>

#include "BaseLib/Algorithm.h"
#include "BaseLib/ConfigTreeUtil.h"
#include "BaseLib/MPI.h"
#include "MaterialLib/MPL/Medium.h"
#include "MathLib/LinAlg/Eigen/EigenVector.h"
#include "MeshLib/Mesh.h"
#include "PhreeqcIOData/AqueousSolution.h"
#include "PhreeqcIOData/ChemicalSystem.h"
#include "PhreeqcIOData/Dump.h"
#include "PhreeqcIOData/EquilibriumReactant.h"
#include "PhreeqcIOData/Exchange.h"
#include "PhreeqcIOData/KineticReactant.h"
#include "PhreeqcIOData/Knobs.h"
#include "PhreeqcIOData/Output.h"
#include "PhreeqcIOData/ReactionRate.h"
#include "PhreeqcIOData/Surface.h"
#include "PhreeqcIOData/UserPunch.h"

namespace ChemistryLib
{
namespace PhreeqcIOData
{
namespace
{
template <class... Ts>
struct overloaded : Ts...
{
    using Ts::operator()...;
};
template <class... Ts>
overloaded(Ts...) -> overloaded<Ts...>;

// Helper to iterate lines in a string_view without copying
class StringViewLineIterator
{
public:
    StringViewLineIterator(std::string_view const text, std::size_t pos = 0)
        : text_(text), pos_(pos)
    {
    }

    bool getline(std::string_view& line)
    {
        if (pos_ >= text_.size())
        {
            return false;
        }

        if (const auto newline_pos = text_.find('\n', pos_);
            newline_pos == std::string_view::npos)
        {
            line = text_.substr(pos_);
            pos_ = text_.size();
        }
        else
        {
            line = text_.substr(pos_, newline_pos - pos_);
            pos_ = newline_pos + 1;
        }

        // Remove trailing \r if present (Windows line endings)
        if (!line.empty() && line.back() == '\r')
        {
            line.remove_suffix(1);
        }

        return true;
    }

    void skip(int const num_lines)
    {
        for (int i = 0; i < num_lines && pos_ < text_.size(); ++i)
        {
            const auto newline_pos = text_.find('\n', pos_);
            pos_ = (newline_pos == std::string_view::npos) ? text_.size()
                                                           : newline_pos + 1;
        }
    }

private:
    std::string_view text_;
    std::size_t pos_;
};

std::vector<std::string> extractItemsFromLine(std::string_view const line)
{
    std::vector<std::string> items;
    std::string line_str(line);
    boost::trim_if(line_str, boost::is_any_of("\t "));
    boost::algorithm::split(items, line_str, boost::is_any_of("\t "),
                            boost::token_compress_on);
    return items;
}

std::vector<double> parseAndFilterChemicalData(
    std::string_view const line,
    std::vector<int> const& dropped_item_ids,
    std::size_t const chemical_system_id)
{
    std::vector<double> accepted_items;
    std::vector<std::string> const items = extractItemsFromLine(line);
    for (int item_id = 0; item_id < static_cast<int>(items.size()); ++item_id)
    {
        if (std::find(dropped_item_ids.begin(), dropped_item_ids.end(),
                      item_id) != dropped_item_ids.end())
        {
            continue;
        }
        double value;
        try
        {
            value = std::stod(items[item_id]);
        }
        catch (const std::invalid_argument& e)
        {
            OGS_FATAL(
                "Invalid argument. Could not convert string '{:s}' to "
                "double for chemical system {:d}, column {:d}. "
                "Exception '{:s}' was thrown.",
                items[item_id], chemical_system_id + 1, item_id, e.what());
        }
        catch (const std::out_of_range& e)
        {
            OGS_FATAL(
                "Out of range error. Could not convert string "
                "'{:s}' to double for chemical system {:d}, column "
                "{:d}. Exception '{:s}' was thrown.",
                items[item_id], chemical_system_id + 1, item_id, e.what());
        }
        accepted_items.push_back(value);
    }
    return accepted_items;
}

template <typename DataBlock>
std::ostream& operator<<(std::ostream& os,
                         std::vector<DataBlock> const& data_blocks)
{
    std::copy(data_blocks.begin(), data_blocks.end(),
              std::ostream_iterator<DataBlock>(os));
    return os;
}

struct ClampingStats
{
    std::size_t n_clamped = 0;
    std::size_t n_severe_clamped = 0;
    double total_clamped_amount = 0.0;
    double worst_negative_value = 0.0;
    std::string_view worst_component_name;
};

// Concentrations between -tolerance and 0 are floating-point noise from MPI
// partitioning and are silently clamped to zero. Values more negative than
// -tolerance trigger an aggregated WARN.
ClampingStats setAqueousSolution(std::vector<double> const& concentrations,
                                 GlobalIndexType const& chemical_system_id,
                                 AqueousSolution& aqueous_solution,
                                 double const negative_tolerance)
{
    ClampingStats stats;
    auto& components = aqueous_solution.components;
    for (unsigned component_id = 0; component_id < components.size();
         ++component_id)
    {
        double c = concentrations[component_id];
        if (c < 0.0)
        {
            if (c < -negative_tolerance)
            {
                stats.n_severe_clamped++;
                if (c < stats.worst_negative_value)
                {
                    stats.worst_negative_value = c;
                    stats.worst_component_name = components[component_id].name;
                }
            }
            stats.n_clamped++;
            stats.total_clamped_amount += -c;
            c = 0.0;
        }
        components[component_id].amount[chemical_system_id] = c;
    }

    if (stats.n_severe_clamped > 0)
    {
        WARN(
            "{:d} component(s) at chemical system {:d} had concentrations "
            "more negative than the tolerance ({:g} mol/L); worst: '{:s}' "
            "= {:g} mol/L. Clamping to zero.",
            stats.n_severe_clamped, chemical_system_id, -negative_tolerance,
            stats.worst_component_name, stats.worst_negative_value);
    }

    // The transport process carries the H+ activity 10^-pH as the last
    // entry of the concentrations vector for the pH "component".
    aqueous_solution.H_plus_activity[chemical_system_id] =
        concentrations.back();
    return stats;
}

template <typename Reactant>
void initializeReactantMolality(Reactant& reactant,
                                GlobalIndexType const& chemical_system_id,
                                MaterialPropertyLib::Phase const& solid_phase,
                                MaterialPropertyLib::Phase const& liquid_phase,
                                MaterialPropertyLib::Medium const& medium,
                                ParameterLib::SpatialPosition const& pos,
                                double const t)
{
    auto const& solid_constituent = solid_phase.component(reactant.name);

    if (solid_constituent.hasProperty(
            MaterialPropertyLib::PropertyType::molality))
    {
        auto const molality =
            solid_constituent[MaterialPropertyLib::PropertyType::molality]
                .template initialValue<double>(pos, t);

        (*reactant.molality)[chemical_system_id] = molality;
        (*reactant.molality_prev)[chemical_system_id] = molality;
    }
    else
    {
        auto const volume_fraction =
            solid_constituent
                [MaterialPropertyLib::PropertyType::volume_fraction]
                    .template initialValue<double>(pos, t);

        (*reactant.volume_fraction)[chemical_system_id] = volume_fraction;

        (*reactant.volume_fraction_prev)[chemical_system_id] = volume_fraction;

        auto const fluid_density =
            liquid_phase[MaterialPropertyLib::PropertyType::density]
                .template initialValue<double>(pos, t);

        auto const porosity =
            medium[MaterialPropertyLib::PropertyType::porosity]
                .template initialValue<double>(pos, t);

        auto const molar_volume =
            solid_constituent[MaterialPropertyLib::PropertyType::molar_volume]
                .template initialValue<double>(pos, t);

        (*reactant.molality)[chemical_system_id] =
            volume_fraction / fluid_density / porosity / molar_volume;

        (*reactant.molality_prev)[chemical_system_id] =
            (*reactant.molality)[chemical_system_id];
    }
}

template <typename Reactant>
void setReactantMolality(Reactant& reactant,
                         GlobalIndexType const& chemical_system_id,
                         MaterialPropertyLib::Phase const& solid_phase,
                         MaterialPropertyLib::Phase const& liquid_phase,
                         MaterialPropertyLib::VariableArray const& vars,
                         ParameterLib::SpatialPosition const& pos,
                         double const t, double const dt)
{
    auto const& solid_constituent = solid_phase.component(reactant.name);

    if (solid_constituent.hasProperty(
            MaterialPropertyLib::PropertyType::molality))
    {
        (*reactant.molality_prev)[chemical_system_id] =
            (*reactant.molality)[chemical_system_id];

        return;
    }

    auto const volume_fraction =
        (*reactant.volume_fraction)[chemical_system_id];

    (*reactant.volume_fraction_prev)[chemical_system_id] =
        (*reactant.volume_fraction)[chemical_system_id];

    auto const fluid_density =
        liquid_phase[MaterialPropertyLib::PropertyType::density]
            .template value<double>(vars, pos, t, dt);

    auto const molar_volume =
        solid_constituent[MaterialPropertyLib::PropertyType::molar_volume]
            .template value<double>(vars, pos, t, dt);

    (*reactant.molality)[chemical_system_id] =
        volume_fraction / fluid_density / vars.porosity / molar_volume;

    (*reactant.molality_prev)[chemical_system_id] =
        (*reactant.molality)[chemical_system_id];
}

template <typename Site>
void initializeSiteMolality(Site& site,
                            GlobalIndexType const& chemical_system_id,
                            MaterialPropertyLib::Phase const& solid_phase,
                            ParameterLib::SpatialPosition const& pos,
                            double const t)
{
    auto const& solid_constituent = solid_phase.component(site.name);

    auto const molality =
        solid_constituent[MaterialPropertyLib::PropertyType::molality]
            .template initialValue<double>(pos, t);

    (*site.molality)[chemical_system_id] = molality;
}

template <typename Reactant>
void updateReactantVolumeFraction(Reactant& reactant,
                                  GlobalIndexType const& chemical_system_id,
                                  MaterialPropertyLib::Medium const& medium,
                                  ParameterLib::SpatialPosition const& pos,
                                  double const porosity, double const t,
                                  double const dt)
{
    auto const& solid_phase =
        medium.phase(MaterialPropertyLib::PhaseName::Solid);
    auto const& liquid_phase =
        medium.phase(MaterialPropertyLib::PhaseName::AqueousLiquid);

    MaterialPropertyLib::VariableArray vars;

    auto const liquid_density =
        liquid_phase[MaterialPropertyLib::PropertyType::density]
            .template value<double>(vars, pos, t, dt);

    auto const& solid_constituent = solid_phase.component(reactant.name);

    if (solid_constituent.hasProperty(
            MaterialPropertyLib::PropertyType::molality))
    {
        return;
    }

    auto const molar_volume =
        solid_constituent[MaterialPropertyLib::PropertyType::molar_volume]
            .template value<double>(vars, pos, t, dt);

    (*reactant.volume_fraction)[chemical_system_id] +=
        ((*reactant.molality)[chemical_system_id] -
         (*reactant.molality_prev)[chemical_system_id]) *
        liquid_density * porosity * molar_volume;
}

template <typename Reactant>
void setPorosityPostReaction(Reactant& reactant,
                             GlobalIndexType const& chemical_system_id,
                             MaterialPropertyLib::Medium const& medium,
                             double& porosity)
{
    auto const& solid_phase =
        medium.phase(MaterialPropertyLib::PhaseName::Solid);

    auto const& solid_constituent = solid_phase.component(reactant.name);

    if (solid_constituent.hasProperty(
            MaterialPropertyLib::PropertyType::molality))
    {
        return;
    }

    porosity -= ((*reactant.volume_fraction)[chemical_system_id] -
                 (*reactant.volume_fraction_prev)[chemical_system_id]);
}

template <typename Reactant>
static double averageReactantMolality(
    Reactant const& reactant,
    std::vector<GlobalIndexType> const& chemical_system_indices)
{
    double const sum = std::accumulate(
        chemical_system_indices.begin(), chemical_system_indices.end(), 0.0,
        [&](double const s, GlobalIndexType const id)
        { return s + (*reactant.molality)[id]; });
    return sum / chemical_system_indices.size();
}
}  // namespace

extern std::string specifyFileName(std::string const& project_file_name,
                                   std::string const& file_extension);

PhreeqcIO::PhreeqcIO(MeshLib::Mesh const& mesh,
                     GlobalLinearSolver& linear_solver,
                     std::string const& project_file_name,
                     std::string&& database,
                     std::unique_ptr<ChemicalSystem>&& chemical_system,
                     std::vector<ReactionRate>&& reaction_rates,
                     std::unique_ptr<UserPunch>&& user_punch,
                     std::unique_ptr<Output>&& output,
                     std::unique_ptr<Dump>&& dump,
                     Knobs&& knobs,
                     bool const use_stream_mode,
                     int const num_chemistry_threads,
                     double const negative_concentration_tolerance)
    : ChemicalSolverInterface(mesh, linear_solver),
      _phreeqc_input_file(specifyFileName(project_file_name, ".inp")),
      _database(std::move(database)),
      _knobs(std::move(knobs)),
      _reaction_rates(std::move(reaction_rates)),
      _chemical_system(std::move(chemical_system)),
      _user_punch(std::move(user_punch)),
      _output(std::move(output)),
      _dump(std::move(dump)),
      _num_chemical_systems(0),
      _negative_concentration_tolerance(negative_concentration_tolerance),
      num_chemistry_threads_(num_chemistry_threads),
      _use_stream_mode(use_stream_mode)
{
    INFO("Chemistry threads per MPI rank: {}.", num_chemistry_threads_);

    if (_use_stream_mode)
    {
        // Stream mode runs exclusively on the instance pool (even for a single
        // thread) so there is one unified execution path. The standalone
        // phreeqc_instance_id is not used here and stays -1.
        if (num_chemistry_threads_ > 1)
        {
            INFO(
                "Parallel chemistry enabled: {} threads will be used for "
                "PHREEQC calculations.",
                num_chemistry_threads_);
        }
        INFO(
            "PhreeqcIO is configured for stream-based data exchange: input "
            "and output will be exchanged via in-memory strings.");
        instance_pool_ = std::make_unique<PhreeqcInstancePool>(
            _database, std::max(1, num_chemistry_threads_));
    }
    else
    {
        // File mode: create and load the standalone PHREEQC instance and
        // enable file-based selected output. PhreeqcInstancePool::
        // createInstance() OGS_FATALs on failure, so any returned id is valid.
        phreeqc_instance_id =
            PhreeqcIOData::PhreeqcInstancePool::createInstance(_database);
        if (SetSelectedOutputFileOn(phreeqc_instance_id, 1) != IPQ_OK)
        {
            OGS_FATAL(
                "Failed to fly the flag for the specified file {:s} where "
                "phreeqc will write output.",
                _output->basic_output_setups.output_file);
        }
        if (_dump)
        {
            SetDumpFileOn(phreeqc_instance_id, 1);
        }
    }
}

PhreeqcIO::~PhreeqcIO()
{
    // The standalone instance is only created in file mode; in stream mode
    // phreeqc_instance_id stays -1 and there is nothing to destroy.
    if (phreeqc_instance_id >= 0)
    {
        DestroyIPhreeqc(phreeqc_instance_id);
    }
}

void PhreeqcIO::initialize()
{
    _num_chemical_systems = chemical_system_index_map.size();

    _chemical_system->initialize(_num_chemical_systems);

    if (_user_punch)
    {
        _user_punch->initialize(_num_chemical_systems);
    }
}

void PhreeqcIO::initializeChemicalSystemConcrete(
    std::vector<double> const& concentrations,
    GlobalIndexType const& chemical_system_id,
    MaterialPropertyLib::Medium const& medium,
    ParameterLib::SpatialPosition const& pos,
    double const t)
{
    setAqueousSolution(concentrations, chemical_system_id,
                       *_chemical_system->aqueous_solution,
                       _negative_concentration_tolerance);

    auto const& solid_phase =
        medium.phase(MaterialPropertyLib::PhaseName::Solid);
    auto const& liquid_phase =
        medium.phase(MaterialPropertyLib::PhaseName::AqueousLiquid);

    for (auto& kinetic_reactant : _chemical_system->kinetic_reactants)
    {
        initializeReactantMolality(kinetic_reactant, chemical_system_id,
                                   solid_phase, liquid_phase, medium, pos, t);
    }

    for (auto& equilibrium_reactant : _chemical_system->equilibrium_reactants)
    {
        initializeReactantMolality(equilibrium_reactant, chemical_system_id,
                                   solid_phase, liquid_phase, medium, pos, t);
    }

    for (auto& exchanger : _chemical_system->exchangers)
    {
        initializeSiteMolality(exchanger, chemical_system_id, solid_phase, pos,
                               t);
    }

    for (auto& surface_site : _chemical_system->surface)
    {
        if (auto const surface_site_ptr =
                std::get_if<MoleBasedSurfaceSite>(&surface_site))
        {
            initializeSiteMolality(*surface_site_ptr, chemical_system_id,
                                   solid_phase, pos, t);
        }
    }
}

void PhreeqcIO::setChemicalSystemConcrete(
    std::vector<double> const& concentrations,
    GlobalIndexType const& chemical_system_id,
    MaterialPropertyLib::Medium const* medium,
    MaterialPropertyLib::VariableArray const& vars,
    ParameterLib::SpatialPosition const& pos, double const t, double const dt)
{
    auto const clamping_stats = setAqueousSolution(
        concentrations, chemical_system_id, *_chemical_system->aqueous_solution,
        _negative_concentration_tolerance);
    if (clamping_stats.n_clamped > 0)
    {
        DBUG(
            "PhreeqcIO: {:d} concentration(s) clamped to zero before "
            "speciation at chemical system {:d} (total clamped: {:g} mol/L).",
            clamping_stats.n_clamped, chemical_system_id,
            clamping_stats.total_clamped_amount);
    }

    auto const& solid_phase =
        medium->phase(MaterialPropertyLib::PhaseName::Solid);
    auto const& liquid_phase =
        medium->phase(MaterialPropertyLib::PhaseName::AqueousLiquid);

    for (auto& kinetic_reactant : _chemical_system->kinetic_reactants)
    {
        setReactantMolality(kinetic_reactant, chemical_system_id, solid_phase,
                            liquid_phase, vars, pos, t, dt);
    }

    for (auto& equilibrium_reactant : _chemical_system->equilibrium_reactants)
    {
        setReactantMolality(equilibrium_reactant, chemical_system_id,
                            solid_phase, liquid_phase, vars, pos, t, dt);
    }
}

void PhreeqcIO::executeSpeciationCalculation(double const dt)
{
    if (_use_stream_mode)
    {
        // Stream mode always uses the pool (serial or parallel), giving one
        // unified execution path. The pool was created with max(1, threads).
        executeSpeciationCalculationParallel(dt);
        return;
    }

    // File-based data exchange.
    DBUG("Executing speciation with file-based data exchange.");
    writeInputsToFile(dt);
    callPhreeqc();
    readOutputsFromFile();
}

double PhreeqcIO::getConcentration(
    int const component_id, GlobalIndexType const chemical_system_id) const
{
    auto const& aqueous_solution = *_chemical_system->aqueous_solution;
    auto const& components = aqueous_solution.components;
    auto const& H_plus_activity = aqueous_solution.H_plus_activity;

    if (component_id < static_cast<int>(components.size()))
    {
        return components[component_id].amount[chemical_system_id];
    }

    if (component_id != static_cast<int>(components.size()))
    {
        OGS_FATAL(
            "Invalid component_id {:d}: must be in [0, {:d}] "
            "(the last index represents H+ activity).",
            component_id, components.size());
    }

    // H+ activity 10^-pH (the transport-side state vector stores this
    // activity for the pH "component" — see setAqueousSolution).
    return H_plus_activity[chemical_system_id];
}

void PhreeqcIO::setAqueousSolutionsPrevFromDumpFile()
{
    if (!_dump)
    {
        return;
    }

    auto const& dump_file = _dump->dump_file;
    std::ifstream in(dump_file);
    if (!in)
    {
        // return if phreeqc dump file doesn't exist. This happens in
        // the first time step when no dump file is provided by the user.
        return;
    }

    _dump->readDumpFile(in, _num_chemical_systems);

    if (!in)
    {
        OGS_FATAL("Error when reading phreeqc dump file '{:s}'", dump_file);
    }

    in.close();
}

void PhreeqcIO::setAqueousSolutionsPrevFromDumpString(
    std::string_view const dump_content)
{
    if (!_dump)
    {
        return;
    }

    if (dump_content.empty())
    {
        // return if dump content is empty. This happens in
        // the first time step when no dump data is available.
        DBUG(
            "Dump content is empty, skipping aqueous solutions initialization "
            "from dump.");
        return;
    }

    _dump->readDumpFromString(dump_content, _num_chemical_systems);
}

void PhreeqcIO::writeInputsToFile(double const dt)
{
    DBUG("Writing phreeqc inputs into file '{:s}'.", _phreeqc_input_file);
    std::ofstream out(_phreeqc_input_file, std::ofstream::out);

    if (!out)
    {
        OGS_FATAL("Could not open file '{:s}' for writing phreeqc inputs.",
                  _phreeqc_input_file);
    }

    out << std::scientific
        << std::setprecision(std::numeric_limits<double>::max_digits10);
    *this << dt;
    out << *this;

    if (!out)
    {
        OGS_FATAL("Failed in generating phreeqc input file '{:s}'.",
                  _phreeqc_input_file);
    }

    out.close();
}

void PhreeqcIO::writeInputHeader(std::ostream& os) const
{
    bool const fixing_pe = _chemical_system->aqueous_solution->fixing_pe;
    if (fixing_pe)
    {
        os << "PHASES\n"
           << "Fix_pe\n"
           << "e- = e-\n"
           << "log_k 0.0\n\n";
    }

    os << _knobs << "\n";
    os << *_output << "\n";

    if (_user_punch)
    {
        os << *_user_punch << "\n";
    }

    if (!_reaction_rates.empty())
    {
        os << "RATES\n";
        os << _reaction_rates << "\n";
    }
}

void PhreeqcIO::writeSystemBlock(std::ostream& os,
                                 std::size_t const chemical_system_id,
                                 std::size_t const solution_id,
                                 std::size_t const prev_solution_id,
                                 double const dt) const
{
    bool const fixing_pe = _chemical_system->aqueous_solution->fixing_pe;

    os << "SOLUTION " << solution_id << "\n";
    _chemical_system->aqueous_solution->print(os, chemical_system_id);

    if (_dump && !_dump->aqueous_solutions_prev.empty())
    {
        os << _dump->aqueous_solutions_prev[chemical_system_id] << "\n\n";
    }

    os << "USE solution none\n";
    os << "END\n\n";

    os << "USE solution " << solution_id << "\n\n";

    auto const& equilibrium_reactants = _chemical_system->equilibrium_reactants;
    if (!equilibrium_reactants.empty() || fixing_pe)
    {
        os << "EQUILIBRIUM_PHASES " << solution_id << "\n";
        for (auto const& r : equilibrium_reactants)
        {
            r.print(os, chemical_system_id);
        }
        fixing_pe ? os << "Fix_pe " << -_chemical_system->aqueous_solution->pe0
                       << " O2(g)\n\n"
                  : os << "\n";
    }

    auto const& kinetic_reactants = _chemical_system->kinetic_reactants;
    if (!kinetic_reactants.empty())
    {
        os << "KINETICS " << solution_id << "\n";
        for (auto const& k : kinetic_reactants)
        {
            k.print(os, chemical_system_id);
        }
        os << "-steps " << dt << "\n\n";
    }

    auto const& surface = _chemical_system->surface;
    if (!surface.empty())
    {
        // To get the amount of surface species from the previous time step,
        // an equilibration calculation with the previous aqueous solution
        // is needed. The previous aqueous solution is saved using the
        // PHREEQC keyword "DUMP" and stored as SOLUTION_RAW within
        // aqueous_solutions_prev. Along with the PHREEQC keyword 'SURFACE',
        // distinguish between the current and previous solutions via
        // prev_solution_id.
        os << "SURFACE " << solution_id << "\n";
        std::size_t const aq_id =
            (_dump && !_dump->aqueous_solutions_prev.empty()) ? prev_solution_id
                                                              : solution_id;
        os << "-equilibrate with solution " << aq_id << "\n";

        if (std::holds_alternative<DensityBasedSurfaceSite>(surface.front()))
        {
            os << "-sites_units density\n";
        }
        else
        {
            os << "-sites_units absolute\n";
        }

        for (auto const& surface_site : surface)
        {
            std::visit(
                overloaded{
                    [&os](DensityBasedSurfaceSite const& s)
                    {
                        os << s.name << " " << s.site_density << " "
                           << s.specific_surface_area << " " << s.mass << "\n";
                    },
                    [&os, chemical_system_id](MoleBasedSurfaceSite const& s)
                    {
                        os << s.name << " " << (*s.molality)[chemical_system_id]
                           << "\n";
                    }},
                surface_site);
        }

        if (std::holds_alternative<MoleBasedSurfaceSite>(surface.front()))
        {
            os << "-no_edl\n";
        }
        os << "SAVE solution " << solution_id << "\n";
    }

    auto const& exchangers = _chemical_system->exchangers;
    if (!exchangers.empty())
    {
        os << "EXCHANGE " << solution_id << "\n";
        std::size_t const aq_id =
            (_dump && !_dump->aqueous_solutions_prev.empty()) ? prev_solution_id
                                                              : solution_id;
        os << "-equilibrate with solution " << aq_id << "\n";
        for (auto const& exchanger : exchangers)
        {
            exchanger.print(os, chemical_system_id);
        }
        os << "SAVE solution " << solution_id << "\n";
    }

    os << "END\n\n";
}

void PhreeqcIO::updateSystemFromOutputLine(std::string_view const line,
                                           std::size_t const chemical_system_id)
{
    auto const& output = *_output;
    auto accepted_items = parseAndFilterChemicalData(
        line, output.dropped_item_ids, chemical_system_id);
    assert(accepted_items.size() == output.accepted_items.size());
    updateChemicalSystemFromOutput(accepted_items, chemical_system_id);
}

std::ostream& operator<<(std::ostream& os, PhreeqcIO const& phreeqc_io)
{
    phreeqc_io.writeInputHeader(os);

    for (std::size_t chemical_system_id = 0;
         chemical_system_id < phreeqc_io._num_chemical_systems;
         ++chemical_system_id)
    {
        std::size_t const solution_id = chemical_system_id + 1;
        std::size_t const prev_solution_id =
            phreeqc_io._num_chemical_systems + chemical_system_id + 1;
        phreeqc_io.writeSystemBlock(os, chemical_system_id, solution_id,
                                    prev_solution_id, phreeqc_io._dt);
    }

    if (phreeqc_io._dump)
    {
        phreeqc_io._dump->print(os, phreeqc_io._num_chemical_systems);
    }

    return os;
}

void PhreeqcIO::callPhreeqc() const
{
    INFO("Phreeqc: Executing chemical calculation.");
    if (RunFile(phreeqc_instance_id, _phreeqc_input_file.c_str()) != IPQ_OK)
    {
        OutputErrorString(phreeqc_instance_id);
        OGS_FATAL(
            "Failed in performing speciation calculation with the generated "
            "phreeqc input file '{:s}'.",
            _phreeqc_input_file);
    }
}

void PhreeqcIO::readOutputsFromFile()
{
    auto const& basic_output_setups = _output->basic_output_setups;
    auto const& phreeqc_result_file = basic_output_setups.output_file;
    DBUG("Reading phreeqc results from file '{:s}'.", phreeqc_result_file);
    std::ifstream in(phreeqc_result_file);

    if (!in)
    {
        OGS_FATAL("Could not open phreeqc result file '{:s}'.",
                  phreeqc_result_file);
    }

    in >> *this;

    if (!in)
    {
        OGS_FATAL("Error when reading phreeqc result file '{:s}'",
                  phreeqc_result_file);
    }

    in.close();
}

void PhreeqcIO::updateChemicalSystemFromOutput(
    std::vector<double> const& accepted_items, std::size_t chemical_system_id)
{
    auto const& output = *_output;
    auto& aqueous_solution = _chemical_system->aqueous_solution;
    auto& components = aqueous_solution->components;
    auto& equilibrium_reactants = _chemical_system->equilibrium_reactants;
    auto& kinetic_reactants = _chemical_system->kinetic_reactants;

    for (int item_id = 0; item_id < static_cast<int>(accepted_items.size());
         ++item_id)
    {
        auto const& accepted_item = output.accepted_items[item_id];
        auto const& item_name = accepted_item.name;

        auto compare_by_name = [&item_name](auto const& item)
        { return item.name == item_name; };

        switch (accepted_item.item_type)
        {
            case ItemType::pH:
            {
                aqueous_solution->H_plus_activity[chemical_system_id] =
                    std::pow(10, -accepted_items[item_id]);
                break;
            }
            case ItemType::pe:
            {
                (*aqueous_solution->pe)[chemical_system_id] =
                    accepted_items[item_id];
                break;
            }
            case ItemType::Component:
            {
                auto& component = BaseLib::findElementOrError(
                    components, compare_by_name,
                    [&]()
                    {
                        OGS_FATAL("Could not find component '{:s}'.",
                                  item_name);
                    });
                component.amount[chemical_system_id] = accepted_items[item_id];
                break;
            }
            case ItemType::EquilibriumReactant:
            {
                auto const& equilibrium_reactant = BaseLib::findElementOrError(
                    equilibrium_reactants, compare_by_name,
                    [&]()
                    {
                        OGS_FATAL("Could not find equilibrium reactant '{:s}'",
                                  item_name);
                    });
                (*equilibrium_reactant.molality)[chemical_system_id] =
                    accepted_items[item_id];
                break;
            }
            case ItemType::KineticReactant:
            {
                auto const& kinetic_reactant = BaseLib::findElementOrError(
                    kinetic_reactants, compare_by_name,
                    [&]()
                    {
                        OGS_FATAL("Could not find kinetic reactant '{:s}'.",
                                  item_name);
                    });
                (*kinetic_reactant.molality)[chemical_system_id] =
                    accepted_items[item_id];
                break;
            }
            case ItemType::SecondaryVariable:
            {
                assert(_user_punch);
                auto const& secondary_variables =
                    _user_punch->secondary_variables;
                auto const& secondary_variable = BaseLib::findElementOrError(
                    secondary_variables, compare_by_name,
                    [&]()
                    {
                        OGS_FATAL("Could not find secondary variable '{:s}'.",
                                  item_name);
                    });
                (*secondary_variable.value)[chemical_system_id] =
                    accepted_items[item_id];
                break;
            }
        }
    }
}

std::istream& operator>>(std::istream& in, PhreeqcIO& phreeqc_io)
{
    // Skip the headline
    in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    std::string line;

    auto const& surface = phreeqc_io._chemical_system->surface;
    auto const& exchangers = phreeqc_io._chemical_system->exchangers;

    int const num_skipped_lines =
        1 + (!surface.empty() ? 1 : 0) + (!exchangers.empty() ? 1 : 0);

    for (std::size_t chemical_system_id = 0;
         chemical_system_id < phreeqc_io._num_chemical_systems;
         ++chemical_system_id)
    {
        for (int i = 0; i < num_skipped_lines; ++i)
        {
            in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        }

        if (!std::getline(in, line))
        {
            OGS_FATAL(
                "Error when reading calculation result of Solution {:d} "
                "after the reaction.",
                chemical_system_id);
        }

        phreeqc_io.updateSystemFromOutputLine(line, chemical_system_id);
    }

    return in;
}

std::vector<std::string> const PhreeqcIO::getComponentList() const
{
    std::vector<std::string> component_names;
    auto const& components = _chemical_system->aqueous_solution->components;
    std::transform(components.begin(), components.end(),
                   std::back_inserter(component_names),
                   [](auto const& c) { return c.name; });

    component_names.push_back("H");

    return component_names;
}

void PhreeqcIO::updateVolumeFractionPostReaction(
    GlobalIndexType const& chemical_system_id,
    MaterialPropertyLib::Medium const& medium,
    ParameterLib::SpatialPosition const& pos, double const porosity,
    double const t, double const dt)
{
    for (auto& kinetic_reactant : _chemical_system->kinetic_reactants)
    {
        updateReactantVolumeFraction(kinetic_reactant, chemical_system_id,
                                     medium, pos, porosity, t, dt);
    }

    for (auto& equilibrium_reactant : _chemical_system->equilibrium_reactants)
    {
        updateReactantVolumeFraction(equilibrium_reactant, chemical_system_id,
                                     medium, pos, porosity, t, dt);
    }
}

void PhreeqcIO::updatePorosityPostReaction(
    GlobalIndexType const& chemical_system_id,
    MaterialPropertyLib::Medium const& medium,
    double& porosity)
{
    for (auto& kinetic_reactant : _chemical_system->kinetic_reactants)
    {
        setPorosityPostReaction(kinetic_reactant, chemical_system_id, medium,
                                porosity);
    }

    for (auto& equilibrium_reactant : _chemical_system->equilibrium_reactants)
    {
        setPorosityPostReaction(equilibrium_reactant, chemical_system_id,
                                medium, porosity);
    }
}

void PhreeqcIO::computeSecondaryVariable(
    std::size_t const ele_id,
    std::vector<GlobalIndexType> const& chemical_system_indices)
{
    for (auto const& kinetic_reactant : _chemical_system->kinetic_reactants)
    {
        (*kinetic_reactant.mesh_prop_molality)[ele_id] =
            averageReactantMolality(kinetic_reactant, chemical_system_indices);
    }

    for (auto const& equilibrium_reactant :
         _chemical_system->equilibrium_reactants)
    {
        (*equilibrium_reactant.mesh_prop_molality)[ele_id] =
            averageReactantMolality(equilibrium_reactant,
                                    chemical_system_indices);
    }
}

std::string PhreeqcIO::generateInputForSystem(
    std::size_t const chemical_system_id, double const dt) const
{
    std::ostringstream os;
    os << std::scientific
       << std::setprecision(std::numeric_limits<double>::max_digits10);

    writeInputHeader(os);

    // Each pool instance processes exactly one system with solution ID 1.
    // The prev_solution_id follows the same remapping convention as in file
    // mode: num_chemical_systems + chemical_system_id + 1.
    std::size_t const prev_solution_id =
        _num_chemical_systems + chemical_system_id + 1;
    writeSystemBlock(os, chemical_system_id, 1, prev_solution_id, dt);

    // Each pool instance dumps its single solution (ID 1) for the next
    // timestep's dump-restore cycle.
    if (_dump)
    {
        os << "DUMP\n";
        os << "-solution 1\n";
        os << "END\n";
    }

    return os.str();
}

void PhreeqcIO::parseOutputForSystem(std::string_view const output_content,
                                     std::size_t const chemical_system_id)
{
    if (output_content.empty())
    {
        OGS_FATAL("Empty output for chemical system {}.", chemical_system_id);
    }

    StringViewLineIterator line_iter(output_content);
    std::string_view line;

    line_iter.getline(line);  // skip headline

    int const num_skipped_lines =
        1 + (!_chemical_system->surface.empty() ? 1 : 0) +
        (!_chemical_system->exchangers.empty() ? 1 : 0);
    line_iter.skip(num_skipped_lines);

    if (!line_iter.getline(line))
    {
        OGS_FATAL(
            "Error when reading calculation result of Solution {} after the "
            "reaction.",
            chemical_system_id);
    }

    updateSystemFromOutputLine(line, chemical_system_id);
}

void PhreeqcIO::executeSpeciationCalculationParallel(double const dt)
{
    INFO("Phreeqc: Executing parallel chemical calculation with {} threads.",
         num_chemistry_threads_);

    // 1. Pre-generate inputs for all chemical systems (sequential)
    // This must be sequential because it reads from shared data structures
    std::vector<std::string> inputs(_num_chemical_systems);
    for (std::size_t i = 0; i < _num_chemical_systems; ++i)
    {
        inputs[i] = generateInputForSystem(i, dt);
    }

    // 2. Storage for outputs
    std::vector<std::string> outputs(_num_chemical_systems);
    std::vector<std::string> dumps(_num_chemical_systems);

    // 3. Parallel execution
    std::vector<std::size_t> failed_systems;
    std::mutex failed_mutex;

#pragma omp parallel num_threads(num_chemistry_threads_)
    {
#ifdef _OPENMP
        int const thread_id = omp_get_thread_num();
#else
        int const thread_id = 0;
#endif
        int const phreeqc_id = instance_pool_->getInstanceForThread(thread_id);

#pragma omp for schedule(dynamic)
        for (std::size_t i = 0; i < _num_chemical_systems; ++i)
        {
            if (RunString(phreeqc_id, inputs[i].c_str()) != IPQ_OK)
            {
                OutputErrorString(phreeqc_id);
                std::lock_guard<std::mutex> guard(failed_mutex);
                failed_systems.push_back(i);
                continue;
            }

            // Retrieve output string
            const char* output_ptr = GetSelectedOutputString(phreeqc_id);
            if (output_ptr)
            {
                outputs[i] = output_ptr;
            }

            // Retrieve dump string if needed
            if (_dump)
            {
                const char* dump_ptr = GetDumpString(phreeqc_id);
                if (dump_ptr)
                {
                    dumps[i] = dump_ptr;
                }
            }
        }
    }

    std::exception_ptr local_error;
    if (!failed_systems.empty())
    {
        std::sort(failed_systems.begin(), failed_systems.end());
        std::string ids;
        for (auto const id : failed_systems)
        {
            if (!ids.empty())
            {
                ids += ", ";
            }
            ids += std::to_string(id);
        }
        local_error = std::make_exception_ptr(std::runtime_error(
            "Failed in performing speciation calculation for "
            "chemical system(s) " +
            ids + "."));
    }
    BaseLib::MPI::allRanksThrowOrNone<std::runtime_error>(local_error);

    // 4. Parse results (sequential, updates shared state)
    for (std::size_t i = 0; i < _num_chemical_systems; ++i)
    {
        parseOutputForSystem(outputs[i], i);
    }

    // 5. Handle dump data if surfaces/exchangers exist.
    // Each pool instance dumped solution 1; remap to the correct system ID.
    if (_dump && !dumps.empty())
    {
        _dump->aqueous_solutions_prev.resize(_num_chemical_systems);
        for (std::size_t i = 0; i < _num_chemical_systems; ++i)
        {
            if (!dumps[i].empty())
            {
                _dump->readDumpFromStringForSystem(dumps[i], i,
                                                   _num_chemical_systems);
            }
        }
    }
}
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
