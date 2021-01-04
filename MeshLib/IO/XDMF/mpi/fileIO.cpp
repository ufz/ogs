/**
 * \file
 * \author Tobias Meisel
 * \date   2020-12-08
 * \brief  Function specific to execution with MPI, never include directly!!
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "../fileIO.h"

#include <hdf5.h>
#include <mpi.h>

#include "BaseLib/Logging.h"

namespace MeshLib::IO
{
hid_t createFile(std::filesystem::path const& filepath)
{
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Info info = MPI_INFO_NULL;
    hid_t const plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, comm, info);
    hid_t file = H5Fcreate(filepath.string().c_str(), H5F_ACC_TRUNC,
                           H5P_DEFAULT, plist_id);
    H5Pclose(plist_id);
    return file;
}

hid_t openHDF5File(std::filesystem::path const& filepath)
{
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Info info = MPI_INFO_NULL;
    hid_t const plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, comm, info);
    hid_t file = H5Fopen(filepath.string().c_str(), H5F_ACC_RDWR, plist_id);
    H5Pclose(plist_id);
    return file;
}

hid_t createHDF5TransferPolicy()
{
    // property list for collective dataset write
    hid_t io_transfer_property = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(io_transfer_property, H5FD_MPIO_COLLECTIVE);
    return io_transfer_property;
}

}  // namespace MeshLib::IO