/**
 * \file
 * \author Tobias Meisel
 * \date   2020-12-15
 * \brief  Writes vectorized data to HDF File
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once
#include <filesystem.h>
#include <hdf5.h>

#include <map>
#include <memory>
#include <vector>

#include "HdfData.h"

namespace MeshLib::IO
{
using HDFAttributes = std::vector<HdfData>;
struct MeshHdfData
{
    HDFAttributes constant_attributes;
    HDFAttributes variable_attributes;
    std::string name;
};

class HdfWriter final
{
public:
    /**
     * \brief Write file with geometry and topology data. The data
     * itself is held by a structure outside of this class. The writer assumes
     * the data holder to not change during writing
     * @param meshes meta data of meshes to be written
     * @param initial_step number of the step (temporal collection), usually 0,
     * greater 0 with continuation of simulation
     * @param filepath absolute or relative filepath to the hdf5 file
     * @param use_compression if true gzip compression is enabled
     * @param is_file_manager True if process (in parallel execution) is
     * @param num_of_files Number of outputfiles
     * File_Manager
     */
    HdfWriter(std::vector<MeshHdfData> meshes,
              unsigned long long initial_step,
              std::filesystem::path const& filepath,
              bool use_compression,
              bool is_file_manager,
              unsigned int num_of_files);
    /**
     * \brief Writes attributes. The data
     * itself is hold by a structure outside of this class. The writer assumes
     * the data holder to not change during writing and HdfData given to
     * constructor to be still valid
     * @param time time_value of step to be written to temporal collection
     */
    void writeStep(double time);
    ~HdfWriter();

private:
    // internal data holder
    struct HdfMesh;

    std::filesystem::path const _hdf5_filepath;
    hid_t const _file;
    hid_t const _meshes_group;
    std::vector<std::unique_ptr<HdfMesh>> _hdf_meshes;
    std::vector<double> _step_times;
    bool const _use_compression;
    bool const _is_file_manager;
};
}  // namespace MeshLib::IO