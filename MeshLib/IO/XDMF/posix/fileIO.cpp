// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "MeshLib/IO/XDMF/fileIO.h"

#include <hdf5.h>
namespace MeshLib::IO
{
int64_t createFile(std::filesystem::path const& filepath, unsigned int)
{
    return H5Fcreate(filepath.string().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
                     H5P_DEFAULT);
}

int64_t openHDF5File(std::filesystem::path const& filepath)
{
    return H5Fopen(filepath.string().c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
}

int64_t createHDF5TransferPolicy()
{
    return H5P_DEFAULT;
}
}  // namespace MeshLib::IO
