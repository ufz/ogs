/**
 * \file
 * \author Tobias Meisel
 * \date   2020-11-13
 * \brief  XdmfWriter which takes contiguous data and writes 1 xdmf + 1 hdf file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <filesystem.h>

#include <boost/shared_ptr.hpp>

#include "MeshLib/Mesh.h"
#include "XdmfData.h"

class XdmfGridCollection;
class XdmfTopology;
class XdmfGeometry;
class XdmfWriter;
class XdmfDomain;

namespace MeshLib::IO
{
struct Geometry;
struct Topology;
}  // namespace MeshLib::IO

namespace MeshLib::IO
{
class Xdmf3Writer final
{
public:
    /**
     * \brief Write xdmf and h5 file with geometry and topology data.
     * The data itself is held by a structure outside of this class.
     * The writer assumes a constant data holder (data itself can change, memory
     * address is the same)
     * @param geometry contains meta data to coordinates
     * @param topology contains meta data cells
     * @param attributes vector of attributes (each attribute is a OGS property)
     * @param filepath absolute or relative filepath to the hdf5 file
     * @param time_step number of the step (temporal collection)
     */
    Xdmf3Writer(XdmfData const& geometry, XdmfData const& topology,
                std::vector<XdmfData> attributes,
                std::filesystem::path const& filepath, int time_step);
    /**
     * \brief Write attribute data that has modified to previous time step or
     * initial
     * @param time_step number of the step (temporal collection)
     * @param time time value of the current time_step
     */
    void writeStep(int time_step, double time);

private:
    boost::shared_ptr<XdmfGridCollection> _gridCollection;
    boost::shared_ptr<XdmfTopology> _initial_topology;
    boost::shared_ptr<XdmfGeometry> _initial_geometry;
    std::vector<XdmfData> const _attributes;
    boost::shared_ptr<XdmfWriter> _writer;
    boost::shared_ptr<XdmfDomain> _root;
    std::string const _hdf5filename;
};
}  // namespace MeshLib::IO
