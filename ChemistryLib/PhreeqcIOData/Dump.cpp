// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "Dump.h"

#include <iostream>

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
    aqueous_solutions_prev.clear();
    aqueous_solutions_prev.reserve(num_chemical_systems);

    std::string line;
    std::string aqueous_solution_prev;
    std::size_t chemical_system_id = 0;
    while (std::getline(in, line))
    {
        if (line.find("USE reaction_pressure none") != std::string::npos)
        {
            break;
        }

        if (line.find("SOLUTION_RAW") != std::string::npos)
        {
            aqueous_solution_prev =
                "SOLUTION_RAW " +
                std::to_string(num_chemical_systems + chemical_system_id + 1) +
                "\n";
            continue;
        }

        aqueous_solution_prev += line + "\n";

        if (line.find("-gammas") != std::string::npos)
        {
            aqueous_solutions_prev.push_back(aqueous_solution_prev);
            aqueous_solution_prev.clear();
            ++chemical_system_id;
        }
    }
}

void Dump::readDumpFromString(std::string_view dump_content,
                              std::size_t const num_chemical_systems)
{
    aqueous_solutions_prev.clear();
    aqueous_solutions_prev.reserve(num_chemical_systems);

    std::string_view line;
    std::string aqueous_solution_prev;
    std::size_t chemical_system_id = 0;
    std::size_t pos = 0;

    while (pos < dump_content.size())
    {
        // Extract next line
        if (const auto newline_pos = dump_content.find('\n', pos);
            newline_pos == std::string_view::npos)
        {
            line = dump_content.substr(pos);
            pos = dump_content.size();
        }
        else
        {
            line = dump_content.substr(pos, newline_pos - pos);
            pos = newline_pos + 1;
        }

        // Remove trailing \r if present (Windows line endings)
        if (!line.empty() && line.back() == '\r')
        {
            line.remove_suffix(1);
        }

        if (line.find("USE reaction_pressure none") != std::string::npos)
        {
            break;
        }

        if (line.find("SOLUTION_RAW") != std::string::npos)
        {
            aqueous_solution_prev =
                "SOLUTION_RAW " +
                std::to_string(num_chemical_systems + chemical_system_id + 1) +
                "\n";
            continue;
        }

        aqueous_solution_prev += line;
        aqueous_solution_prev += "\n";

        if (line.find("-gammas") != std::string::npos)
        {
            aqueous_solutions_prev.push_back(aqueous_solution_prev);
            aqueous_solution_prev.clear();
            ++chemical_system_id;
        }
    }
}
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
