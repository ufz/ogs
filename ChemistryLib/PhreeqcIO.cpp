/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "PhreeqcIO.h"

#include <IPhreeqc.h>

#include <boost/algorithm/string.hpp>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <numeric>

#include "BaseLib/Algorithm.h"
#include "BaseLib/ConfigTreeUtil.h"
#include "MaterialLib/MPL/Medium.h"
#include "MathLib/LinAlg/Eigen/EigenVector.h"
#include "MathLib/LinAlg/LinAlg.h"
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

template <typename DataBlock>
std::ostream& operator<<(std::ostream& os,
                         std::vector<DataBlock> const& data_blocks)
{
    std::copy(data_blocks.begin(), data_blocks.end(),
              std::ostream_iterator<DataBlock>(os));
    return os;
}

void setAqueousSolution(std::vector<double> const& concentrations,
                        GlobalIndexType const& chemical_system_id,
                        AqueousSolution& aqueous_solution)
{
    GlobalIndexType const offset = aqueous_solution.pH->getRangeBegin();
    GlobalIndexType const global_index = offset + chemical_system_id;

    // components
    auto& components = aqueous_solution.components;
    for (unsigned component_id = 0; component_id < components.size();
         ++component_id)
    {
        MathLib::LinAlg::setLocalAccessibleVector(
            *components[component_id].amount);
        components[component_id].amount->set(global_index,
                                             concentrations[component_id]);
    }

    // pH
    MathLib::LinAlg::setLocalAccessibleVector(*aqueous_solution.pH);
    aqueous_solution.pH->set(global_index, concentrations.back());
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
    auto const& solid_phase = medium.phase("Solid");
    auto const& liquid_phase = medium.phase("AqueousLiquid");

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
    auto const& solid_phase = medium.phase("Solid");

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
                     Knobs&& knobs)
    : ChemicalSolverInterface(mesh, linear_solver),
      _phreeqc_input_file(specifyFileName(project_file_name, ".inp")),
      _database(std::move(database)),
      _chemical_system(std::move(chemical_system)),
      _reaction_rates(std::move(reaction_rates)),
      _user_punch(std::move(user_punch)),
      _output(std::move(output)),
      _dump(std::move(dump)),
      _knobs(std::move(knobs))
{
    // initialize phreeqc instance
    if (CreateIPhreeqc() != phreeqc_instance_id)
    {
        OGS_FATAL(
            "Failed to initialize phreeqc instance, due to lack of memory.");
    }

    // load specified thermodynamic database
    if (LoadDatabase(phreeqc_instance_id, _database.c_str()) != IPQ_OK)
    {
        OGS_FATAL(
            "Failed in loading the specified thermodynamic database file: "
            "{:s}.",
            _database);
    }

    if (SetSelectedOutputFileOn(phreeqc_instance_id, 1) != IPQ_OK)
    {
        OGS_FATAL(
            "Failed to fly the flag for the specified file {:s} where phreeqc "
            "will write output.",
            _output->basic_output_setups.output_file);
    }

    if (_dump)
    {
        // Chemical composition of the aqueous solution of last time step will
        // be written into .dmp file once the second function argument is set to
        // one.
        SetDumpFileOn(phreeqc_instance_id, 1);
    }
}

PhreeqcIO::~PhreeqcIO()
{
    DestroyIPhreeqc(phreeqc_instance_id);
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
                       *_chemical_system->aqueous_solution);

    auto const& solid_phase = medium.phase("Solid");
    auto const& liquid_phase = medium.phase("AqueousLiquid");

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
    setAqueousSolution(concentrations, chemical_system_id,
                       *_chemical_system->aqueous_solution);

    auto const& solid_phase = medium->phase("Solid");
    auto const& liquid_phase = medium->phase("AqueousLiquid");

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
    writeInputsToFile(dt);

    callPhreeqc();

    readOutputsFromFile();
}

double PhreeqcIO::getConcentration(
    int const component_id, GlobalIndexType const chemical_system_id) const
{
    auto const& aqueous_solution = *_chemical_system->aqueous_solution;
    auto& components = aqueous_solution.components;
    auto const& pH = *aqueous_solution.pH;

    if (component_id < static_cast<int>(components.size()))
    {
        MathLib::LinAlg::setLocalAccessibleVector(
            *components[component_id].amount);

        return components[component_id].amount->get(chemical_system_id);
    }

    // pH
    MathLib::LinAlg::setLocalAccessibleVector(*aqueous_solution.pH);

    return pH.get(chemical_system_id);
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
        // return if phreeqc dump file doesn't exist. Normally, this happens in
        // the first time step.
        return;
    }

    _dump->readDumpFile(in, _num_chemical_systems);

    if (!in)
    {
        OGS_FATAL("Error when reading phreeqc dump file '{:s}'", dump_file);
    }

    in.close();
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
        << std::setprecision(std::numeric_limits<double>::digits10);
    out << (*this << dt);

    if (!out)
    {
        OGS_FATAL("Failed in generating phreeqc input file '{:s}'.",
                  _phreeqc_input_file);
    }

    out.close();
}

std::ostream& operator<<(std::ostream& os, PhreeqcIO const& phreeqc_io)
{
    bool const fixing_pe =
        phreeqc_io._chemical_system->aqueous_solution->fixing_pe;
    if (fixing_pe)
    {
        os << "PHASES\n"
           << "Fix_pe\n"
           << "e- = e-\n"
           << "log_k 0.0\n\n";
    }

    os << phreeqc_io._knobs << "\n";

    os << *phreeqc_io._output << "\n";

    auto const& user_punch = phreeqc_io._user_punch;
    if (user_punch)
    {
        os << *user_punch << "\n";
    }

    auto const& reaction_rates = phreeqc_io._reaction_rates;
    if (!reaction_rates.empty())
    {
        os << "RATES\n";
        os << reaction_rates << "\n";
    }

    for (std::size_t chemical_system_id = 0;
         chemical_system_id < phreeqc_io._num_chemical_systems;
         ++chemical_system_id)
    {
        os << "SOLUTION " << chemical_system_id + 1 << "\n";
        phreeqc_io._chemical_system->aqueous_solution->print(
            os, chemical_system_id);

        auto const& dump = phreeqc_io._dump;
        if (dump)
        {
            auto const& aqueous_solutions_prev = dump->aqueous_solutions_prev;
            if (!aqueous_solutions_prev.empty())
            {
                os << aqueous_solutions_prev[chemical_system_id] << "\n\n";
            }
        }

        os << "USE solution none\n";
        os << "END\n\n";

        os << "USE solution " << chemical_system_id + 1 << "\n\n";

        auto const& equilibrium_reactants =
            phreeqc_io._chemical_system->equilibrium_reactants;
        if (!equilibrium_reactants.empty() || fixing_pe)
        {
            os << "EQUILIBRIUM_PHASES " << chemical_system_id + 1 << "\n";
            for (auto const& equilibrium_reactant : equilibrium_reactants)
            {
                equilibrium_reactant.print(os, chemical_system_id);
            }
            fixing_pe
                ? os << "Fix_pe "
                     << -phreeqc_io._chemical_system->aqueous_solution->pe0
                     << " O2(g)\n\n"
                : os << "\n";
        }

        auto const& kinetic_reactants =
            phreeqc_io._chemical_system->kinetic_reactants;
        if (!kinetic_reactants.empty())
        {
            os << "KINETICS " << chemical_system_id + 1 << "\n";
            for (auto const& kinetic_reactant : kinetic_reactants)
            {
                kinetic_reactant.print(os, chemical_system_id);
            }
            os << "-steps " << phreeqc_io._dt << "\n\n";
        }

        auto const& surface = phreeqc_io._chemical_system->surface;
        if (!surface.empty())
        {
            os << "SURFACE " << chemical_system_id + 1 << "\n";
            std::size_t aqueous_solution_id =
                dump->aqueous_solutions_prev.empty()
                    ? chemical_system_id + 1
                    : phreeqc_io._num_chemical_systems + chemical_system_id + 1;
            os << "-equilibrate with solution " << aqueous_solution_id << "\n";

            // print unit
            if (std::holds_alternative<DensityBasedSurfaceSite>(
                    surface.front()))
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
                    overloaded{[&os](DensityBasedSurfaceSite const& s)
                               {
                                   os << s.name << " " << s.site_density << " "
                                      << s.specific_surface_area << " "
                                      << s.mass << "\n";
                               },
                               [&os, chemical_system_id](
                                   MoleBasedSurfaceSite const& s) {
                                   os << s.name << " "
                                      << (*s.molality)[chemical_system_id]
                                      << "\n";
                               }},
                    surface_site);
            }

            // overlook the effect of the buildup of charges onto the surface
            if (std::holds_alternative<MoleBasedSurfaceSite>(surface.front()))
            {
                os << "-no_edl\n";
            }

            os << "SAVE solution " << chemical_system_id + 1 << "\n";
        }

        auto const& exchangers = phreeqc_io._chemical_system->exchangers;
        if (!exchangers.empty())
        {
            os << "EXCHANGE " << chemical_system_id + 1 << "\n";
            std::size_t const aqueous_solution_id =
                dump->aqueous_solutions_prev.empty()
                    ? chemical_system_id + 1
                    : phreeqc_io._num_chemical_systems + chemical_system_id + 1;
            os << "-equilibrate with solution " << aqueous_solution_id << "\n";
            for (auto const& exchanger : exchangers)
            {
                exchanger.print(os, chemical_system_id);
            }
            os << "SAVE solution " << chemical_system_id + 1 << "\n";
        }

        os << "END\n\n";
    }

    auto const& dump = phreeqc_io._dump;
    if (dump)
    {
        dump->print(os, phreeqc_io._num_chemical_systems);
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

std::istream& operator>>(std::istream& in, PhreeqcIO& phreeqc_io)
{
    // Skip the headline
    in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    std::string line;
    auto const& output = *phreeqc_io._output;
    auto const& dropped_item_ids = output.dropped_item_ids;

    auto const& surface = phreeqc_io._chemical_system->surface;
    auto const& exchangers = phreeqc_io._chemical_system->exchangers;

    int const num_skipped_lines =
        1 + (!surface.empty() ? 1 : 0) + (!exchangers.empty() ? 1 : 0);

    auto& equilibrium_reactants =
        phreeqc_io._chemical_system->equilibrium_reactants;
    auto& kinetic_reactants = phreeqc_io._chemical_system->kinetic_reactants;

    for (std::size_t chemical_system_id = 0;
         chemical_system_id < phreeqc_io._num_chemical_systems;
         ++chemical_system_id)
    {
        // Skip equilibrium calculation result of initial solution
        for (int i = 0; i < num_skipped_lines; ++i)
        {
            in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        }

        // Get calculation result of the solution after the reaction
        if (!std::getline(in, line))
        {
            OGS_FATAL(
                "Error when reading calculation result of Solution {:d} "
                "after the reaction.",
                chemical_system_id);
        }

        std::vector<double> accepted_items;
        std::vector<std::string> items;
        boost::trim_if(line, boost::is_any_of("\t "));
        boost::algorithm::split(items, line, boost::is_any_of("\t "),
                                boost::token_compress_on);
        for (int item_id = 0; item_id < static_cast<int>(items.size());
             ++item_id)
        {
            if (std::find(dropped_item_ids.begin(), dropped_item_ids.end(),
                          item_id) == dropped_item_ids.end())
            {
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
                        items[item_id], chemical_system_id + 1, item_id,
                        e.what());
                }
                catch (const std::out_of_range& e)
                {
                    OGS_FATAL(
                        "Out of range error. Could not convert string "
                        "'{:s}' to double for chemical system {:d}, column "
                        "{:d}. Exception '{:s}' was thrown.",
                        items[item_id], chemical_system_id + 1, item_id,
                        e.what());
                }
                accepted_items.push_back(value);
            }
        }
        assert(accepted_items.size() == output.accepted_items.size());

        auto& aqueous_solution = phreeqc_io._chemical_system->aqueous_solution;
        auto& components = aqueous_solution->components;
        auto& user_punch = phreeqc_io._user_punch;

        GlobalIndexType const offset = aqueous_solution->pH->getRangeBegin();
        GlobalIndexType const global_index = offset + chemical_system_id;

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
                    // Update pH value
                    MathLib::LinAlg::setLocalAccessibleVector(
                        *aqueous_solution->pH);
                    aqueous_solution->pH->set(
                        global_index, std::pow(10, -accepted_items[item_id]));
                    break;
                }
                case ItemType::pe:
                {
                    // Update pe value
                    (*aqueous_solution->pe)[chemical_system_id] =
                        accepted_items[item_id];
                    break;
                }
                case ItemType::Component:
                {
                    // Update component concentrations
                    auto const& component = BaseLib::findElementOrError(
                        components, compare_by_name,
                        [&]() {
                            OGS_FATAL("Could not find component '{:s}'.",
                                      item_name);
                        });
                    MathLib::LinAlg::setLocalAccessibleVector(
                        *component.amount);
                    component.amount->set(global_index,
                                          accepted_items[item_id]);
                    break;
                }
                case ItemType::EquilibriumReactant:
                {
                    // Update amounts of equilibrium reactant
                    auto const& equilibrium_reactant =
                        BaseLib::findElementOrError(
                            equilibrium_reactants, compare_by_name,
                            [&]()
                            {
                                OGS_FATAL(
                                    "Could not find equilibrium reactant "
                                    "'{:s}'",
                                    item_name);
                            });
                    (*equilibrium_reactant.molality)[chemical_system_id] =
                        accepted_items[item_id];
                    break;
                }
                case ItemType::KineticReactant:
                {
                    // Update amounts of kinetic reactants
                    auto const& kinetic_reactant = BaseLib::findElementOrError(
                        kinetic_reactants, compare_by_name,
                        [&]() {
                            OGS_FATAL("Could not find kinetic reactant '{:s}'.",
                                      item_name);
                        });
                    (*kinetic_reactant.molality)[chemical_system_id] =
                        accepted_items[item_id];
                    break;
                }
                case ItemType::SecondaryVariable:
                {
                    assert(user_punch);
                    auto const& secondary_variables =
                        user_punch->secondary_variables;
                    // Update values of secondary variables
                    auto const& secondary_variable =
                        BaseLib::findElementOrError(
                            secondary_variables, compare_by_name,
                            [&]() {
                                OGS_FATAL(
                                    "Could not find secondary variable '{:s}'.",
                                    item_name);
                            });
                    (*secondary_variable.value)[chemical_system_id] =
                        accepted_items[item_id];
                    break;
                }
            }
        }
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
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
