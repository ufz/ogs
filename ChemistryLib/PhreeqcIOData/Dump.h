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

    std::string const dump_file;
    std::vector<std::string> aqueous_solutions_prev;
    std::string dump_content_stream;  // In-memory dump data from PhreeqC
};
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
