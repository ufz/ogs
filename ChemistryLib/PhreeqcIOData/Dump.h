/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <filesystem>
#include <iosfwd>
#include <string>
#include <vector>

#include "BaseLib/Error.h"

namespace ChemistryLib
{
namespace PhreeqcIOData
{
struct Dump
{
    explicit Dump(std::string const& project_file_name)
        : dump_file(project_file_name + ".dmp")
    {
        try
        {
            if (std::filesystem::remove(dump_file))
            {
                INFO("Deleted the redundant phreeqc dump file {:s}\n",
                     dump_file);
            }
        }
        catch (std::filesystem::filesystem_error const& e)
        {
            ERR("filesystem error: {:s}\n", e.what());
        }
    }

    void print(std::ostream& os, std::size_t const num_chemical_systems) const;

    void readDumpFile(std::istream& in, std::size_t const num_chemical_systems);

    std::string const dump_file;
    std::vector<std::string> aqueous_solutions_prev;
};
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
