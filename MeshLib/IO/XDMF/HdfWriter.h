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
#include <vector>

#include "HdfData.h"

namespace MeshLib::IO
{
class HdfWriter final
{
public:
    /**
     * \brief Write file with geometry and topology data. The data
     * itself is held by a structure outside of this class. The writer assumes
     * the data holder to not change during writing
     * @param constant_attributes vector of constant attributes (each attribute
     * is a OGS mesh property), geometry and topology are considered as constant
     * attributes
     * @param variable_attributes vector of variable attributes (each attribute
     * is a OGS mesh property
     * @param step number of the step (temporal collection)
     * @param filepath absolute or relative filepath to the hdf5 file
     * @param use_compression if true gzip compression is enabled
     */
    HdfWriter(std::vector<HdfData> constant_attributes,
              std::vector<HdfData>
                  variable_attributes,
              int const step,
              std::filesystem::path const& filepath,
              bool const use_compression);

    /**
     * \brief Writes attributes. The data
     * itself is hold by a structure outside of this class. The writer assumes
     * the data holder to not change during writing and HdfData given to
     * constructor to be still valid
     * @param step number of the step (temporal collection)
     * @return true = success, false = error
     */
    bool writeStep(int step) const;
    ~HdfWriter();

private:
    std::vector<HdfData> const _variable_attributes;
    std::filesystem::path const _hdf5_filepath;
    bool const _use_compression;
    hid_t const _file;
    hid_t _group;
    std::map<std::string, hid_t> _datasets;
};
}  // namespace MeshLib::IO