// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <string>
#include <vector>

#pragma once

struct CommandLineArguments final
{
    std::string reference_path;
    std::string project;
    std::vector<std::string> xml_patch_file_names;
    std::string outdir;
    std::string mesh_dir;
    std::string script_dir;
    std::string log_level;
    bool log_parallel;
    bool write_prj;
    bool nonfatal;
    bool reference_path_is_set;
#ifndef _WIN32
    bool enable_fpe_is_set;
#endif  // _WIN32
};

CommandLineArguments parseCommandLineArguments(
    int argc, char* argv[], bool const exit_on_exception = true);
