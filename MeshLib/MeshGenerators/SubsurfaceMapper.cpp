/**
 * \file   SubsurfaceMapper.cpp
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

#include "SubsurfaceMapper.h"

#include <vector>
#include <fstream>

#include "FileIO/AsciiRasterInterface.h"

#include "Raster.h"

#include "Mesh.h"
#include "Node.h"
#include "Elements/Element.h"
#include "MeshQuality/MeshValidation.h"

SubsurfaceMapper::SubsurfaceMapper()
: _elevation_epsilon(0.0001)
{
}

bool SubsurfaceMapper::createLayers(MeshLib::Mesh const& mesh, std::vector<std::string> const& raster_paths, double noDataReplacementValue)
{
    if (mesh.getDimension() != 2 || !allRastersExist(raster_paths))
        return false;

    std::vector<GeoLib::Raster const*> rasters;
    rasters.reserve(raster_paths.size());
    for (auto path = raster_paths.begin(); path != raster_paths.end(); ++path)
        rasters.push_back(FileIO::AsciiRasterInterface::getRasterFromASCFile(*path));

    bool result = createRasterLayers(mesh, rasters, noDataReplacementValue);
    std::for_each(rasters.begin(), rasters.end(), [](GeoLib::Raster const*const raster){ delete raster; });
    return true;
}

MeshLib::Mesh* SubsurfaceMapper::getMesh(std::string const& mesh_name) const 
{
    if (_nodes.empty() || _elements.empty())
        return nullptr;

    MeshLib::Mesh* result (new MeshLib::Mesh(mesh_name, _nodes, _elements)); 
    MeshLib::MeshValidation::removeUnusedMeshNodes(*result);
    return result;
}

double SubsurfaceMapper::calcEpsilon(GeoLib::Raster const& high, GeoLib::Raster const& low)
{
    const double max (*std::max_element(high.begin(), high.end()));
    const double min (*std::min_element( low.begin(),  low.end()));
    return ((max-min)*1e-07);
}

bool SubsurfaceMapper::allRastersExist(std::vector<std::string> const& raster_paths) const
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

void SubsurfaceMapper::cleanUpOnError()
{
    std::for_each(_nodes.begin(), _nodes.end(), [](MeshLib::Node *node) { delete node; });
    std::for_each(_elements.begin(), _elements.end(), [](MeshLib::Element *elem) { delete elem; });
}
