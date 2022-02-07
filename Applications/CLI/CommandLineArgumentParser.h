/**
 * \brief  Declaration of CommandLineArgumentParser.
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <string>
#include <vector>

#pragma once

struct CommandLineArgumentParser final
{
    CommandLineArgumentParser(int argc, char* argv[]);

    std::string reference_path;
    std::string project;
    std::vector<std::string> xml_patch_file_names;
    std::string outdir;
    std::string mesh_dir;
    std::string log_level;
    bool write_prj;
    bool nonfatal;
    bool reference_path_is_set;
#ifndef _WIN32
    bool enable_fpe_is_set;
#endif  // _WIN32
};
