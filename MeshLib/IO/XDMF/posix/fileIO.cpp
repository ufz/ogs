/**
 * \file
 * \author Tobias Meisel
 * \date   2020-12-08
 * \brief  Function specific to execution without MPI
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MeshLib/IO/XDMF/fileIO.h"

#include <hdf5.h>
namespace MeshLib::IO
{
int64_t createFile(std::filesystem::path const& filepath)
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
