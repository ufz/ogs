/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "PhreeqcIO.h"

#include <iphreeqc/src/src/IPhreeqc.h>

#include <boost/algorithm/string.hpp>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <numeric>

#include "BaseLib/Algorithm.h"
#include "BaseLib/ConfigTreeUtil.h"
#include "MathLib/LinAlg/Eigen/EigenVector.h"
#include "MathLib/LinAlg/LinAlg.h"
#include "MeshLib/Mesh.h"
#include "PhreeqcIOData/AqueousSolution.h"
#include "PhreeqcIOData/ChemicalSystem.h"
#include "PhreeqcIOData/Dump.h"
#include "PhreeqcIOData/EquilibriumReactant.h"
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
template <typename DataBlock>
std::ostream& operator<<(std::ostream& os,
                         std::vector<DataBlock> const& data_blocks)
{
    std::copy(data_blocks.begin(), data_blocks.end(),
              std::ostream_iterator<DataBlock>(os));
    return os;
}
}  // namespace

PhreeqcIO::PhreeqcIO(std::string const project_file_name,
                     MeshLib::Mesh const& mesh,
                     std::string&& database,
                     std::unique_ptr<ChemicalSystem>&& chemical_system,
                     std::vector<ReactionRate>&& reaction_rates,
                     std::vector<SurfaceSite>&& surface,
                     std::unique_ptr<UserPunch>&& user_punch,
                     std::unique_ptr<Output>&& output,
                     std::unique_ptr<Dump>&& dump,
                     Knobs&& knobs)
    : _phreeqc_input_file(project_file_name + "_phreeqc.inp"),
      _mesh(mesh),
      _database(std::move(database)),
      _chemical_system(std::move(chemical_system)),
      _reaction_rates(std::move(reaction_rates)),
      _surface(std::move(surface)),
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

void PhreeqcIO::initialize()
{
    _num_chemical_systems = std::accumulate(
        begin(chemical_system_index_map), end(chemical_system_index_map), 0,
        [](auto result, auto const& chemical_system_index) {
            return result + chemical_system_index.size();
        });

    _chemical_system->initialize(_num_chemical_systems);

    if (_user_punch)
    {
        _user_punch->initialize(_num_chemical_systems);
    }
}

void PhreeqcIO::executeInitialCalculation(
    std::vector<GlobalVector> const& interpolated_process_solutions)
{
    setAqueousSolution(interpolated_process_solutions);

    writeInputsToFile();

    execute();

    readOutputsFromFile();
}

std::vector<GlobalVector*> PhreeqcIO::getIntPtProcessSolutions() const
{
    auto const& aqueous_solution = _chemical_system->aqueous_solution;
    std::vector<GlobalVector*> int_pt_process_solutions;
    int_pt_process_solutions.reserve(aqueous_solution->components.size() + 1);

    std::transform(aqueous_solution->components.begin(),
                   aqueous_solution->components.end(),
                   std::back_inserter(int_pt_process_solutions),
                   [](auto const& c) { return c.amount.get(); });

    int_pt_process_solutions.push_back(aqueous_solution->pH.get());

    return int_pt_process_solutions;
}

void PhreeqcIO::doWaterChemistryCalculation(
    std::vector<GlobalVector> const& interpolated_process_solutions,
    double const dt)
{
    setAqueousSolution(interpolated_process_solutions);

    setAqueousSolutionsPrevFromDumpFile();

    writeInputsToFile(dt);

    execute();

    readOutputsFromFile();
}

void PhreeqcIO::setAqueousSolution(
    std::vector<GlobalVector> const& interpolated_process_solutions)
{
    auto& aqueous_solution = _chemical_system->aqueous_solution;

    // components
    auto& components = aqueous_solution->components;
    for (unsigned component_id = 0; component_id < components.size();
         ++component_id)
    {
        MathLib::LinAlg::copy(interpolated_process_solutions[component_id + 1],
                              *components[component_id].amount);
    }

    // pH
    MathLib::LinAlg::copy(interpolated_process_solutions.back(),
                          *aqueous_solution->pH);
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
        OGS_FATAL("Could not open phreeqc dump file '{:s}'.", dump_file);
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
        os << "RATES" << "\n";
        os << reaction_rates << "\n";
    }

    for (std::size_t mesh_item_id = 0;
         mesh_item_id < phreeqc_io.chemical_system_index_map.size();
         ++mesh_item_id)
    {
        auto const& local_chemical_system =
            phreeqc_io.chemical_system_index_map[mesh_item_id];
        for (std::size_t ip = 0; ip < local_chemical_system.size(); ++ip)
        {
            auto const& chemical_system_id = local_chemical_system[ip];

            os << "SOLUTION " << chemical_system_id + 1 << "\n";
            phreeqc_io._chemical_system->aqueous_solution->print(
                os, chemical_system_id);

            auto const& dump = phreeqc_io._dump;
            if (dump)
            {
                auto const& aqueous_solutions_prev =
                    dump->aqueous_solutions_prev;
                if (!aqueous_solutions_prev.empty())
                {
                    os << aqueous_solutions_prev[chemical_system_id] << "\n\n";
                }
            }

            os << "USE solution none" << "\n";
            os << "END" << "\n\n";

            os << "USE solution " << chemical_system_id + 1 << "\n\n";

            auto const& equilibrium_reactants =
                phreeqc_io._chemical_system->equilibrium_reactants;
            if (!equilibrium_reactants.empty())
            {
                os << "EQUILIBRIUM_PHASES " << chemical_system_id + 1 << "\n";
                for (auto const& equilibrium_reactant : equilibrium_reactants)
                {
                    equilibrium_reactant.print(os, chemical_system_id);
                }
                os << "\n";
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
                os << "-steps " << phreeqc_io._dt << "\n" << "\n";
            }

            auto const& surface = phreeqc_io._surface;
            if (!surface.empty())
            {
                os << "SURFACE " << chemical_system_id + 1 << "\n";
                std::size_t aqueous_solution_id =
                    dump->aqueous_solutions_prev.empty()
                        ? chemical_system_id + 1
                        : phreeqc_io._num_chemical_systems +
                              chemical_system_id + 1;
                os << "-equilibrate with solution " << aqueous_solution_id << "\n";
                os << "-sites_units DENSITY" << "\n";
                os << surface << "\n";
                os << "SAVE solution " << chemical_system_id + 1 << "\n";
            }

            os << "END" << "\n\n";
        }
    }

    auto const& dump = phreeqc_io._dump;
    if (dump)
    {
        dump->print(os, phreeqc_io._num_chemical_systems);
    }

    return os;
}

void PhreeqcIO::execute()
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

    auto const& surface = phreeqc_io._surface;
    int const num_skipped_lines = surface.empty() ? 1 : 2;

    auto& equilibrium_reactants =
        phreeqc_io._chemical_system->equilibrium_reactants;
    for (auto& equilibrium_reactant : equilibrium_reactants)
    {
        std::fill(equilibrium_reactant.amount_avg->begin(),
                  equilibrium_reactant.amount_avg->end(), 0.);
    }

    auto& kinetic_reactants = phreeqc_io._chemical_system->kinetic_reactants;
    for (auto& kinetic_reactant : kinetic_reactants)
    {
        std::fill(kinetic_reactant.amount_avg->begin(),
                  kinetic_reactant.amount_avg->end(), 0.);
    }

    for (std::size_t mesh_item_id = 0;
         mesh_item_id < phreeqc_io.chemical_system_index_map.size();
         ++mesh_item_id)
    {
        auto const& local_chemical_system =
            phreeqc_io.chemical_system_index_map[mesh_item_id];
        auto const num_local_chemical_system = local_chemical_system.size();
        for (std::size_t ip = 0; ip < num_local_chemical_system; ++ip)
        {
            auto const& chemical_system_id = local_chemical_system[ip];

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
                    "after "
                    "the reaction.",
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
                            "Invalid argument. Could not convert string '{:s}' "
                            "to double for chemical system {:d}, column {:d}. "
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

            auto& aqueous_solution =
                phreeqc_io._chemical_system->aqueous_solution;
            auto& components = aqueous_solution->components;
            auto& user_punch = phreeqc_io._user_punch;
            for (int item_id = 0;
                 item_id < static_cast<int>(accepted_items.size());
                 ++item_id)
            {
                auto const& accepted_item = output.accepted_items[item_id];
                auto const& item_name = accepted_item.name;

                auto compare_by_name = [&item_name](auto const& item) {
                    return item.name == item_name;
                };

                switch (accepted_item.item_type)
                {
                    case ItemType::pH:
                    {
                        // Update pH value
                        aqueous_solution->pH->set(
                            chemical_system_id,
                            std::pow(10, -accepted_items[item_id]));
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
                        auto& component = BaseLib::findElementOrError(
                            components.begin(), components.end(),
                            compare_by_name,
                            "Could not find component '" + item_name + "'.");
                        component.amount->set(chemical_system_id,
                                              accepted_items[item_id]);
                        break;
                    }
                    case ItemType::EquilibriumReactant:
                    {
                        // Update amounts of equilibrium reactant
                        auto& equilibrium_reactant =
                            BaseLib::findElementOrError(
                                equilibrium_reactants.begin(),
                                equilibrium_reactants.end(), compare_by_name,
                                "Could not find equilibrium reactant '" +
                                    item_name + "'.");
                        (*equilibrium_reactant.amount)[chemical_system_id] =
                            accepted_items[item_id];
                        (*equilibrium_reactant.amount_avg)[mesh_item_id] +=
                            accepted_items[item_id];
                        break;
                    }
                    case ItemType::KineticReactant:
                    {
                        // Update amounts of kinetic reactants
                        auto& kinetic_reactant = BaseLib::findElementOrError(
                            kinetic_reactants.begin(), kinetic_reactants.end(),
                            compare_by_name,
                            "Could not find kinetic reactant '" + item_name +
                                "'.");
                        (*kinetic_reactant.amount)[chemical_system_id] =
                            accepted_items[item_id];
                        (*kinetic_reactant.amount_avg)[mesh_item_id] +=
                            accepted_items[item_id];
                        break;
                    }
                    case ItemType::SecondaryVariable:
                    {
                        assert(user_punch);
                        auto& secondary_variables =
                            user_punch->secondary_variables;
                        // Update values of secondary variables
                        auto& secondary_variable = BaseLib::findElementOrError(
                            secondary_variables.begin(),
                            secondary_variables.end(), compare_by_name,
                            "Could not find secondary variable '" + item_name +
                                "'.");
                        (*secondary_variable.value)[chemical_system_id] =
                            accepted_items[item_id];
                        break;
                    }
                }
            }
        }

        for (auto& equilibrium_reactant : equilibrium_reactants)
        {
            (*equilibrium_reactant.amount_avg)[mesh_item_id] /=
                num_local_chemical_system;
        }

        for (auto& kinetic_reactant : kinetic_reactants)
        {
            (*kinetic_reactant.amount_avg)[mesh_item_id] /=
                num_local_chemical_system;
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
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
