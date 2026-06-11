// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
extern std::string specifyFileName(std::string const& project_file_name,
                                   std::string const& file_extension);

struct Dump
{
    Dump() = default;

    explicit Dump(std::string const& project_file_name)
        : dump_file(specifyFileName(project_file_name, ".dmp"))
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

    void readDumpFromString(std::string_view dump_content,
                            std::size_t const num_chemical_systems);

    /// Parse dump data for a single chemical system from parallel execution.
    /// In parallel mode, each PHREEQC instance produces dump for solution 1,
    /// which needs to be mapped to the correct chemical_system_id.
    void readDumpFromStringForSystem(std::string_view dump_content,
                                     std::size_t const chemical_system_id,
                                     std::size_t const num_chemical_systems);

    std::string dump_file;
    std::vector<std::string> aqueous_solutions_prev;
};
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
