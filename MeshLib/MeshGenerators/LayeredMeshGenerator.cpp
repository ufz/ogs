/**
 * \file
 * \author Karsten Rink
 * \date   2014-09-18
 * \brief  Implementation of the SubsurfaceMapper class.
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "LayeredMeshGenerator.h"

#include <vector>
#include <fstream>

#include "BaseLib/Logging.h"

#include "GeoLib/Raster.h"

#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/PropertyVector.h"
#include "MeshLib/Properties.h"
#include "MeshLib/MeshSearch/NodeSearch.h"
#include "MeshLib/MeshEditing/RemoveMeshComponents.h"

bool LayeredMeshGenerator::createLayers(
    MeshLib::Mesh const& mesh,
    std::vector<GeoLib::Raster const*> const& rasters,
    double minimum_thickness,
    double noDataReplacementValue)
{
    if (mesh.getDimension() != 2)
    {
        return false;
    }

    bool result = createRasterLayers(mesh, rasters, minimum_thickness, noDataReplacementValue);
    std::for_each(rasters.begin(), rasters.end(), [](GeoLib::Raster const*const raster){ delete raster; });
    return result;
}

std::unique_ptr<MeshLib::Mesh>
LayeredMeshGenerator::getMesh(std::string const& mesh_name) const
{
    if (nodes_.empty() || elements_.empty())
    {
        return nullptr;
    }

    MeshLib::Properties properties;
    if (materials_.size() == elements_.size())
    {
        auto* const materials = properties.createNewPropertyVector<int>(
            "MaterialIDs", MeshLib::MeshItemType::Cell);
        assert(materials != nullptr);
        materials->reserve(materials_.size());
        std::copy(materials_.cbegin(),
                  materials_.cend(),
                  std::back_inserter(*materials));
    }
    else
    {
        WARN(
            "Skipping MaterialID information, number of entries does not match "
            "element number");
    }

    std::unique_ptr<MeshLib::Mesh> result(new MeshLib::Mesh(mesh_name, nodes_, elements_, properties));
    MeshLib::NodeSearch ns(*result);
    if (ns.searchUnused() > 0) {
        std::unique_ptr<MeshLib::Mesh> new_mesh(
            MeshLib::removeNodes(*result, ns.getSearchedNodeIDs(), mesh_name));
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

    if ((std::abs(elevation - raster.getHeader().no_data) <
         std::numeric_limits<double>::epsilon()) ||
        (elevation - last_layer_node[2] < minimum_thickness_))
    {
        return new MeshLib::Node(last_layer_node);
    }

    return new MeshLib::Node(dem_node[0], dem_node[1], elevation, new_node_id);
}
