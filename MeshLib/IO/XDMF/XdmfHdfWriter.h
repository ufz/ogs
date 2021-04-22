/**
 * \file
 * \author Tobias Meisel
 * \date   2020-11-13
 * \brief  XdmfWriter which create contiguous data for geometry and topology
 * and writes this and all attributes to 1 xdmf + 1 hdf file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <filesystem.h>

#include <set>

#include "HdfWriter.h"
#include "MeshLib/Mesh.h"
#include "Xdmf3Writer.h"

namespace MeshLib::IO
{
class XdmfHdfWriter final
{
public:
    /**
     * \brief Write xdmf and h5 file with geometry and topology data.
     * @param mesh Mesh or NodePartitionedMesh to be written to file(s)
     * @param filepath absolute or relative filepath to the hdf5 file
     * @param time_step number of the step (temporal collection)
     */
    XdmfHdfWriter(MeshLib::Mesh const& mesh,
                  std::filesystem::path const& filepath, int time_step,
                  std::set<std::string> variable_output_names, bool use_compression);
    /**
     * \brief Write attribute data that has modified to previous time step or
     * initial
     * @param time_step number of the step (temporal collection)
     * @param time time value of the current time_step
     */
    void writeStep(int time_step, double time) const;

private:
    // hdf_writer must be destructed before xdmf_writer
    std::unique_ptr<Xdmf3Writer> _xdmf_writer;
    std::unique_ptr<HdfWriter> _hdf_writer;
};
}  // namespace MeshLib::IO