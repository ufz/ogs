// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <filesystem>

namespace MeshLib::IO
{
int64_t createFile(std::filesystem::path const& filepath, unsigned int n_files);
int64_t openHDF5File(std::filesystem::path const& filepath,
                     unsigned int n_files);
int64_t createHDF5TransferPolicy();
}  // namespace MeshLib::IO
