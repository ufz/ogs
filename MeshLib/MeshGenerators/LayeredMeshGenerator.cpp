/**
 * \file   LayeredMeshGenerator.cpp
 * \author Karsten Rink
 * \date   2014-09-18
 * \brief  Implementation of the SubsurfaceMapper class.
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "LayeredMeshGenerator.h"

#include <vector>
#include <fstream>

#include <logog/include/logog.hpp>

#include "FileIO/AsciiRasterInterface.h"

#include "GeoLib/Raster.h"

#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/PropertyVector.h"
#include "MeshLib/Properties.h"
#include "MeshLib/MeshQuality/MeshValidation.h"
#include "MeshLib/MeshSearch/NodeSearch.h"
#include "MeshLib/MeshEditing/RemoveMeshComponents.h"

LayeredMeshGenerator::LayeredMeshGenerator()
: _elevation_epsilon(0.0001), _minimum_thickness(std::numeric_limits<float>::epsilon())
{
}

bool LayeredMeshGenerator::createLayers(MeshLib::Mesh const& mesh,
                                        std::vector<std::string> const& raster_paths,
                                        double minimum_thickness,
                                        double noDataReplacementValue)
{
    if (mesh.getDimension() != 2 || !allRastersExist(raster_paths))
        return false;

    std::vector<GeoLib::Raster const*> rasters;
    rasters.reserve(raster_paths.size());
    for (auto path = raster_paths.begin(); path != raster_paths.end(); ++path)
        rasters.push_back(FileIO::AsciiRasterInterface::getRasterFromASCFile(*path));

    bool result = createRasterLayers(mesh, rasters, minimum_thickness, noDataReplacementValue);
    std::for_each(rasters.begin(), rasters.end(), [](GeoLib::Raster const*const raster){ delete raster; });
    return result;
}

std::unique_ptr<MeshLib::Mesh>
LayeredMeshGenerator::getMesh(std::string const& mesh_name) const
{
	if (_nodes.empty() || _elements.empty())
		return nullptr;

	MeshLib::Properties properties;
	if (_materials.size() == _elements.size())
	{
		boost::optional<MeshLib::PropertyVector<int> &> materials =
			properties.createNewPropertyVector<int>("MaterialIDs", MeshLib::MeshItemType::Cell);
		materials->resize(_materials.size());
		std::copy(_materials.cbegin(), _materials.cend(), materials->begin());
	}
	else
		WARN ("Skipping MaterialID information, number of entries does not match element number");

	std::unique_ptr<MeshLib::Mesh> result(new MeshLib::Mesh(mesh_name, _nodes, _elements, properties));
	MeshLib::NodeSearch ns(*result.get());
	if (ns.searchUnused() > 0) {
		std::unique_ptr<MeshLib::Mesh> new_mesh(MeshLib::removeNodes(
			*result.get(), ns.getSearchedNodeIDs(), mesh_name));
		return new_mesh;
	}
	return result;
}

double LayeredMeshGenerator::calcEpsilon(GeoLib::Raster const& low, GeoLib::Raster const& high)
{
    const double max (*std::max_element(high.begin(), high.end()));
    const double min (*std::min_element( low.begin(),  low.end()));
    return ((max-min)*1e-07);
}

MeshLib::Node* LayeredMeshGenerator::getNewLayerNode(MeshLib::Node const& dem_node,
                                                     MeshLib::Node const& last_layer_node,
                                                     GeoLib::Raster const& raster,
                                                     std::size_t new_node_id) const
{
    double const elevation = std::min(raster.interpolateValueAtPoint(dem_node), dem_node[2]);

    if ((std::abs(elevation - raster.getNoDataValue()) < std::numeric_limits<double>::epsilon()) ||
        (elevation - last_layer_node[2] < _minimum_thickness))
        return new MeshLib::Node(last_layer_node);
    else
        return new MeshLib::Node(dem_node[0], dem_node[1], elevation, new_node_id);
}

bool LayeredMeshGenerator::allRastersExist(std::vector<std::string> const& raster_paths) const
{
    for (auto raster = raster_paths.begin(); raster != raster_paths.end(); ++raster)
    {
        std::ifstream file_stream (*raster, std::ifstream::in);
        if (!file_stream.good())
            return false;
        file_stream.close();
    }
    return true;
}

void LayeredMeshGenerator::cleanUpOnError()
{
    std::for_each(_nodes.begin(), _nodes.end(), [](MeshLib::Node *node) { delete node; });
    std::for_each(_elements.begin(), _elements.end(), [](MeshLib::Element *elem) { delete elem; });
}
