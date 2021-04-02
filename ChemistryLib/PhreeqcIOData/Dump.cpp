/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

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
    std::string aqueous_solution_prev_;
    std::size_t chemical_system_id = 0;
    while (std::getline(in, line))
    {
        if (line.find("USE reaction_pressure none") != std::string::npos)
        {
            break;
        }

        if (line.find("SOLUTION_RAW") != std::string::npos)
        {
            aqueous_solution_prev_ =
                "SOLUTION_RAW " +
                std::to_string(num_chemical_systems + chemical_system_id + 1) +
                "\n";
            continue;
        }

        aqueous_solution_prev_ += line + "\n";

        if (line.find("-gammas") != std::string::npos)
        {
            aqueous_solutions_prev.push_back(aqueous_solution_prev_);
            aqueous_solution_prev_.clear();
            ++chemical_system_id;
        }
    }
}
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
