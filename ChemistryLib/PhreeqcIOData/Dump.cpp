// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "Dump.h"

#include <sstream>

#include "BaseLib/Logging.h"

namespace
{
/// Parse a PHREEQC dump stream and return accumulated SOLUTION_RAW blocks.
/// Each block's ID is remapped to num_chemical_systems + (start_id + index) +
/// 1, matching the convention expected by the next timestep's PHREEQC input.
std::vector<std::string> parseDumpContent(
    std::istream& in,
    std::size_t const num_chemical_systems,
    std::size_t const start_id = 0)
{
    std::vector<std::string> results;
    std::string line;
    std::string current;
    std::size_t id = start_id;

    while (std::getline(in, line))
    {
        // Strip trailing \r for robustness when the stream was produced on a
        // platform with CRLF line endings.
        if (!line.empty() && line.back() == '\r')
        {
            line.pop_back();
        }
        if (line.find("USE reaction_pressure none") != std::string::npos)
        {
            break;
        }
        if (line.find("SOLUTION_RAW") != std::string::npos)
        {
            current = "SOLUTION_RAW " +
                      std::to_string(num_chemical_systems + id + 1) + "\n";
            continue;
        }
        current += line;
        current += "\n";
        if (line.find("-gammas") != std::string::npos)
        {
            results.push_back(std::move(current));
            current.clear();
            ++id;
        }
    }
    return results;
}
}  // namespace

namespace ChemistryLib
{
namespace PhreeqcIOData
{
void Dump::print(std::ostream& os, std::size_t const num_chemical_systems) const
{
    os << "DUMP"
       << "\n";
    os << "-file " << dump_file << "\n";
    os << "-append false"
       << "\n";
    os << "-solution 1-" << num_chemical_systems << "\n";
    os << "END"
       << "\n";
}

void Dump::readDumpFile(std::istream& in,
                        std::size_t const num_chemical_systems)
{
    aqueous_solutions_prev = parseDumpContent(in, num_chemical_systems);
}

void Dump::readDumpFromString(std::string_view const dump_content,
                              std::size_t const num_chemical_systems)
{
    std::istringstream in{std::string(dump_content)};
    readDumpFile(in, num_chemical_systems);
}

void Dump::readDumpFromStringForSystem(std::string_view const dump_content,
                                       std::size_t const chemical_system_id,
                                       std::size_t const num_chemical_systems)
{
    if (aqueous_solutions_prev.size() < num_chemical_systems)
    {
        aqueous_solutions_prev.resize(num_chemical_systems);
    }

    std::istringstream in{std::string(dump_content)};
    auto parsed =
        parseDumpContent(in, num_chemical_systems, chemical_system_id);
    if (parsed.empty())
    {
        WARN("No SOLUTION_RAW found in dump for chemical system {}.",
             chemical_system_id);
        return;
    }
    aqueous_solutions_prev[chemical_system_id] = std::move(parsed[0]);
}
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
