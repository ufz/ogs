/**
 * \file
 * \author Tobias Meisel
 * \date   2020-12-08
 * \brief  Function specific to execution with MPI, never include directly!!
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MeshLib/IO/XDMF/fileIO.h"

#include <hdf5.h>
#include <mpi.h>

#include "BaseLib/Logging.h"
#include "getCommunicator.h"

using namespace std::string_literals;
namespace MeshLib::IO
{
std::filesystem::path partitionFilename(
    std::filesystem::path const& basic_filepath, int const file_group)
{
    std::string const filename = (file_group > 0)
                                     ? basic_filepath.stem().string() + "_"s +
                                           std::to_string(file_group) +
                                           basic_filepath.extension().string()
                                     : basic_filepath.filename().string();
    std::filesystem::path const filepathwithextension =
        basic_filepath.parent_path() / filename;
    DBUG("HDF Filepath: {:s}.", filepathwithextension.string());
    return filepathwithextension;
};

hid_t createFile(std::filesystem::path const& filepath,
                 unsigned int const n_files)
{
    auto const communicator = getCommunicator(n_files);
    MPI_Comm const comm = communicator.mpi_communicator;
    MPI_Info const info = MPI_INFO_NULL;
    hid_t const plist_id = H5Pcreate(H5P_FILE_ACCESS);

    H5Pset_fapl_mpio(plist_id, comm, info);
    std::filesystem::path const partition_filename =
        partitionFilename(filepath, communicator.color);
    hid_t file = H5Fcreate(partition_filename.string().c_str(), H5F_ACC_TRUNC,
                           H5P_DEFAULT, plist_id);
    H5Pclose(plist_id);

    return file;
}

hid_t openHDF5File(std::filesystem::path const& filepath,
                   unsigned int const n_files)
{
    MPI_Comm const comm = getCommunicator(n_files).mpi_communicator;
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
