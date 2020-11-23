/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "Xdmf3Writer.h"

#include <Xdmf.hpp>
#include <XdmfAttribute.hpp>
#include <XdmfDomain.hpp>
#include <XdmfGeometryType.hpp>
#include <XdmfGridCollection.hpp>
#include <XdmfGridCollectionType.hpp>
#include <XdmfHDF5Controller.hpp>
#include <XdmfHeavyDataDescription.hpp>
#include <XdmfInformation.hpp>
#include <XdmfTopologyType.hpp>
#include <XdmfUnstructuredGrid.hpp>
#include <XdmfWriter.hpp>
#include <string>

#include "InfoLib/GitInfo.h"
#include "XdmfData.h"
#include "writeHDF5.h"

using namespace std::string_literals;

static std::string getTimeSection(int const step, std::string const& name)
{
    return "t_"s + std::to_string(step) + "/"s + name;
}

using namespace MeshLib::IO;

static boost::shared_ptr<XdmfGeometry> getLightGeometry(
    std::string const& hdf5filename, int const step, Geometry const& geometry)
{
    auto xdmf_geometry = XdmfGeometry::New();
    xdmf_geometry->setType(XdmfGeometryType::XYZ());
    boost::shared_ptr<XdmfHDF5Controller> geometry_controller =
        XdmfHDF5Controller::New(hdf5filename,
                                getTimeSection(step, "geometry"),
                                XdmfArrayType::Float64(),
                                geometry.starts,
                                geometry.strides,
                                geometry.vdims,
                                geometry.vdims);
    xdmf_geometry->setHeavyDataController(geometry_controller);
    return xdmf_geometry;
}

static boost::shared_ptr<XdmfTopology> getLightTopology(
    std::string const& hdf5filename, int const step, Topology const& topology)
{
    auto xdmf_topology = XdmfTopology::New();
    xdmf_topology->setType(XdmfTopologyType::Mixed());
    auto topology_controller =
        XdmfHDF5Controller::New(hdf5filename,
                                getTimeSection(step, "topology"),
                                XdmfArrayType::Int32(),
                                topology.starts,
                                topology.strides,
                                topology.vdims,
                                topology.vdims);
    xdmf_topology->setHeavyDataController(topology_controller);
    return xdmf_topology;
}

static boost::shared_ptr<XdmfAttribute> getLightAttribute(
    std::string const& hdf5filename,
    int const step,
    AttributeMeta const& attribute)
{
    auto attribute_controller =
        XdmfHDF5Controller::New(hdf5filename,
                                getTimeSection(step, attribute.name),
                                attribute.data_type,
                                attribute.starts,
                                attribute.strides,
                                attribute.vdims,
                                attribute.vdims);

    auto xdmf_attribute = XdmfAttribute::New();
    xdmf_attribute->setCenter(attribute.attribute_center);
    xdmf_attribute->setName(attribute.name);
    xdmf_attribute->setHeavyDataController(attribute_controller);
    return xdmf_attribute;
}

namespace MeshLib::IO
{
Xdmf3Writer::Xdmf3Writer(std::filesystem::path const& filepath,
                         Geometry const& geometry,
                         Topology const& topology,
                         std::vector<AttributeMeta> const& attributes,
                         int const timestep)
    : _attributes(std::move(attributes)),
      _hdf5filepath(filepath.parent_path() /
                    (std::string(filepath.stem()) + ".h5"))
{
    std::filesystem::path xdmf_filepath =
        filepath.parent_path() / (std::string(filepath.stem()) + ".xdmf");
    auto ret_hdf5 = writeHDF5Initial(geometry.flattend_values,
                                     geometry.vldims,
                                     topology.flattend_values,
                                     timestep,
                                     _hdf5filepath);
    // If we find a library for compression we use it
    _use_compression = ret_hdf5.second;

    _initial_geometry =
        getLightGeometry(_hdf5filepath.filename(), timestep, geometry);
    _initial_topology =
        getLightTopology(_hdf5filepath.filename(), timestep, topology);

    _writer = XdmfWriter::New(xdmf_filepath);
    _writer->setMode(XdmfWriter::DistributedHeavyData);

    auto version = XdmfInformation::New();
    version->setKey(GitInfoLib::GitInfo::OGS_VERSION);
    version->setValue(GitInfoLib::GitInfo::ogs_version);

    auto grid_collection = XdmfGridCollection::New();
    grid_collection->setType(XdmfGridCollectionType::Temporal());

    _root = XdmfDomain::New();
    _root->insert(version);
    _root->insert(grid_collection);
}

void Xdmf3Writer::WriteStep(int const time_step, double const time)
{
    auto grid = XdmfUnstructuredGrid::New();
    grid->setGeometry(_initial_geometry);
    grid->setTopology(_initial_topology);
    grid->setTime(XdmfTime::New(time));

    for (auto const& attribute_data : _attributes)
    {
        writeHDF5Step(_hdf5filepath,
                      time_step,
                      attribute_data.name,
                      attribute_data.data_start,
                      attribute_data.vldims,
                      attribute_data.data_type,
                      _use_compression);
        grid->insert(getLightAttribute(
            _hdf5filepath.filename(), time_step, attribute_data));
    }

    auto gridcollection = _root->getGridCollection(0);
    gridcollection->insert(grid);
    _root->accept(_writer);
}
}  // namespace MeshLib::IO