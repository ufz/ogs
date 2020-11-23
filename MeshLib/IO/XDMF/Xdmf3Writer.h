/**
 * \file
 * \author Tobias Meisel
 * \date   2020-11-13
 * \brief  XdmfWriter which takes contiguous data and writes 1 xdmf + 1 hdf file
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <filesystem.h>

#include <boost/shared_ptr.hpp>

#include "XdmfData.h"

class XdmfGridCollection;
class XdmfTopology;
class XdmfGeometry;
class XdmfWriter;
class XdmfDomain;

namespace MeshLib::IO
{
class Geometry;
class Topology;
}  // namespace MeshLib::IO

namespace MeshLib::IO
{
class Xdmf3Writer
{
public:
    /**
     * \brief Write xdmf and h5 file with geometry and topology data.
     * @param filepath absolute or relativ filepath to the hdf5 file
     * @param time_step number of the step (temporal collection)
     * @param geometry contains point coordinates
     * @param topology contains cells
     * @param attributes vector of attributes (each attribute is a OGS property)
     */
    Xdmf3Writer(std::filesystem::path const& filepath,
                Geometry const& geometry,
                Topology const& topology,
                std::vector<AttributeMeta> const& attributes,
                int time_step);
    /**
     * \brief Write attribute data that has modified to previous time step or
     * initial
     * @param filepath absolute or relativ filepath to the hdf5 file
     * @param time_step number of the step (temporal collection)
     * @param time time value of the current time_step
     */

    void WriteStep(int time_step, double time);

private:
    boost::shared_ptr<XdmfGridCollection> _gridCollection;
    boost::shared_ptr<XdmfTopology> _initial_topology;
    boost::shared_ptr<XdmfGeometry> _initial_geometry;
    std::vector<AttributeMeta> _attributes;
    boost::shared_ptr<XdmfWriter> _writer;
    boost::shared_ptr<XdmfDomain> _root;
    std::filesystem::path const _hdf5filepath;
    bool _use_compression;
};
}  // namespace MeshLib::IO