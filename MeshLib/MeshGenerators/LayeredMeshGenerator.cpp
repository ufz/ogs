/**
 * \file   LayeredMeshGenerator.cpp
 * \author Karsten Rink
 * \date   2014-09-18
 * \brief  Implementation of the SubsurfaceMapper class.
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "LayeredMeshGenerator.h"

#include <vector>
#include <fstream>

#include "FileIO/AsciiRasterInterface.h"

#include "GeoLib/Raster.h"

#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/MeshQuality/MeshValidation.h"

LayeredMeshGenerator::LayeredMeshGenerator()
: _elevation_epsilon(0.0001)
{
}

bool LayeredMeshGenerator::createLayers(MeshLib::Mesh const& mesh, std::vector<std::string> const& raster_paths, double noDataReplacementValue)
{
    if (mesh.getDimension() != 2 || !allRastersExist(raster_paths))
        return false;

    std::vector<GeoLib::Raster const*> rasters;
    rasters.reserve(raster_paths.size());
    for (auto path = raster_paths.begin(); path != raster_paths.end(); ++path)
        rasters.push_back(FileIO::AsciiRasterInterface::getRasterFromASCFile(*path));

    bool result = createRasterLayers(mesh, rasters, noDataReplacementValue);
    std::for_each(rasters.begin(), rasters.end(), [](GeoLib::Raster const*const raster){ delete raster; });
    return result;
}

MeshLib::Mesh* LayeredMeshGenerator::getMesh(std::string const& mesh_name) const
{
    if (_nodes.empty() || _elements.empty())
        return nullptr;

    MeshLib::Mesh* result (new MeshLib::Mesh(mesh_name, _nodes, _elements));
    MeshLib::MeshValidation::removeUnusedMeshNodes(*result);
    return result;
}

double LayeredMeshGenerator::calcEpsilon(GeoLib::Raster const& low, GeoLib::Raster const& high)
{
    const double max (*std::max_element(high.begin(), high.end()));
    const double min (*std::min_element( low.begin(),  low.end()));
    return ((max-min)*1e-07);
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
