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

#include <sstream>
#include <vector>

namespace BaseLib
{
/**
 * @brief Applies includes and patch files to project file.
 *
 * @param prj_stream The processed prj as a stringstream.
 * @param filepath    The prj file.
 * @param patch_files Optional patch files.
 */
void prepareProjectFile(std::stringstream& prj_stream,
                        const std::string& filepath,
                        const std::vector<std::string>& patch_files = {});
}  // namespace BaseLib
