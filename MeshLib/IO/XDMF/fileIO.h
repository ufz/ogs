/**
 * \file
 * \author Tobias Meisel
 * \date   2020-12-08
 * \brief  Dispatches HDF5 functions specific to execution platform (w/o MPI).
 * There are multiple implementation to this interface!
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <filesystem>

namespace MeshLib::IO
{
int64_t createFile(std::filesystem::path const& filepath, unsigned int n_files);
int64_t openHDF5File(std::filesystem::path const& filepath,
                     unsigned int n_files);
int64_t createHDF5TransferPolicy();
}  // namespace MeshLib::IO
