/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
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

using namespace MeshLib::IO;
using namespace std::string_literals;

boost::shared_ptr<const XdmfAttributeCenter> elemTypeOGS2XDMF(
    MeshLib::MeshItemType const elem_type)
{
    std::map<MeshLib::MeshItemType,
             boost::shared_ptr<const XdmfAttributeCenter>>
        mesh_item_type_ogs2xdmf = {
            {MeshLib::MeshItemType::Cell, XdmfAttributeCenter::Cell()},
            {MeshLib::MeshItemType::Edge, XdmfAttributeCenter::Edge()},
            {MeshLib::MeshItemType::Face, XdmfAttributeCenter::Face()},
            {MeshLib::MeshItemType::Node, XdmfAttributeCenter::Node()},
            {MeshLib::MeshItemType::IntegrationPoint,
             XdmfAttributeCenter::Other()}};

    return mesh_item_type_ogs2xdmf.at(elem_type);
}

static std::string getTimeSection(int const step, std::string const& name)
{
    return "t_"s + std::to_string(step) + "/"s + name;
}

static boost::shared_ptr<XdmfGeometry> getLightGeometry(
    std::string const& hdf5filename, int const step, XdmfData const& geometry)
{
    auto xdmf_geometry = XdmfGeometry::New();
    xdmf_geometry->setType(XdmfGeometryType::XYZ());
    boost::shared_ptr<XdmfHDF5Controller> geometry_controller =
        XdmfHDF5Controller::New(hdf5filename,
                                getTimeSection(step, "geometry"),
                                XdmfArrayType::Float64(),
                                geometry.starts,
                                geometry.strides,
                                geometry.block_dims,
                                geometry.block_dims);
    xdmf_geometry->setHeavyDataController(geometry_controller);
    return xdmf_geometry;
}

static boost::shared_ptr<XdmfTopology> getLightTopology(
    std::string const& hdf5filename, int const step, XdmfData const& topology)
{
    auto xdmf_topology = XdmfTopology::New();
    xdmf_topology->setType(XdmfTopologyType::Mixed());
    auto topology_controller =
        XdmfHDF5Controller::New(hdf5filename,
                                getTimeSection(step, "topology"),
                                XdmfArrayType::Int32(),
                                topology.starts,
                                topology.strides,
                                topology.block_dims,
                                topology.block_dims);
    xdmf_topology->setHeavyDataController(topology_controller);
    return xdmf_topology;
}

static boost::shared_ptr<XdmfAttribute> getLightAttribute(
    std::string const& hdf5filename, int const step, XdmfData const& attribute)
{
    auto const attribute_controller =
        XdmfHDF5Controller::New(hdf5filename,
                                getTimeSection(step, attribute.name),
                                attribute.data_type,
                                attribute.starts,
                                attribute.strides,
                                attribute.block_dims,
                                attribute.block_dims);

    auto const xdmf_attribute = XdmfAttribute::New();
    auto const center = elemTypeOGS2XDMF(*(attribute.attribute_center));
    xdmf_attribute->setCenter(center);
    xdmf_attribute->setName(attribute.name);
    xdmf_attribute->setHeavyDataController(attribute_controller);
    return xdmf_attribute;
}

namespace MeshLib::IO
{
Xdmf3Writer::Xdmf3Writer(XdmfData const& geometry, XdmfData const& topology,
                         std::vector<XdmfData> constant_attributes,
                         std::vector<XdmfData> variable_attributes,
                         std::filesystem::path const& filepath,
                         int const time_step)
    : _variable_attributes(std::move(variable_attributes)),
      _hdf5filename(filepath.stem().string() + ".h5")
{
    _initial_geometry = getLightGeometry(_hdf5filename, time_step, geometry);
    _initial_topology = getLightTopology(_hdf5filename, time_step, topology);

    for (auto const& attribute : constant_attributes)
    {
        _constant_attributes.push_back(
            getLightAttribute(_hdf5filename, time_step, attribute));
    }

    _writer = XdmfWriter::New(filepath.string());
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

Xdmf3Writer::~Xdmf3Writer()
{
    _root->accept(_writer);
}

void Xdmf3Writer::writeStep(int const time_step, double const time)
{
    auto grid = XdmfUnstructuredGrid::New();
    grid->setGeometry(_initial_geometry);
    grid->setTopology(_initial_topology);

    for (auto const& constant_attribute : _constant_attributes)
    {
        grid->insert(constant_attribute);
    }

    grid->setTime(XdmfTime::New(time));

    for (auto const& attribute : _variable_attributes)
    {
        grid->insert(getLightAttribute(_hdf5filename, time_step, attribute));
    }

    auto grid_collection = _root->getGridCollection(0);
    grid_collection->insert(grid);

}
}  // namespace MeshLib::IO
