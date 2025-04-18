/**
 * \file
 * \author Karsten Rink
 * \date   2012-09-25
 * \brief  Definition of the GeoMapper class.
 *
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <cstddef>
#include <vector>

#include "GeoLib/GEOObjects.h"
#include "GeoLib/Grid.h"
#include "GeoLib/Point.h"
#include "MathLib/Point3d.h"

namespace MeshLib {
    class Mesh;
    class Node;
}

namespace GeoLib {
    class Raster;
}

namespace MeshGeoToolsLib {

/**
 * \brief A set of tools for mapping the elevation of geometric objects
 */
class GeoMapper
{
public:
    GeoMapper(GeoLib::GEOObjects &geo_objects, const std::string &geo_name);
    ~GeoMapper();

    /// Maps geometry based on a raster file
    void mapOnDEM(std::unique_ptr<GeoLib::Raster const> raster);

    /**
     * Maps the geometry based on the given mesh file. The elevation value of all geometric
     * points are modified such that they are located on the mesh surface.
     */
    void mapOnMesh(MeshLib::Mesh const*const mesh);

    /// Maps geometry to a constant elevation value
    void mapToConstantValue(double value);

    /**
     * Maps the geometry based on the given mesh file. I.e., all geometric
     * points are assigned an elevation value on the mesh surface. Additional
     * points are inserted whenever a polyline from the original geometry
     * intersects a mesh node or the edge of a mesh element.
     * A new geometry with the given name is inserted into _geo_objects.
     * \param mesh          Mesh the geometry is mapped on
     */
    void advancedMapOnMesh(MeshLib::Mesh const& mesh);

private:
    /// Mapping stations, boreholes on a raster or mesh.
    void mapStationData(std::vector<GeoLib::Point*> const& points);

    /// Mapping points on a raster.
    void mapPointDataToDEM(std::vector<GeoLib::Point*> const& points) const;

    /// Mapping points on mesh.
    void mapPointDataToMeshSurface(std::vector<GeoLib::Point*> const& pnts);

    /// Returns the elevation at Point (x,y) based on a mesh. This uses collision detection for triangles and nearest neighbor for quads.
    /// NOTE: This method only returns correct values if the node numbering of the elements is correct!
    double getMeshElevation(double x, double y, double min_val, double max_val) const;

    /// Returns the elevation at Point (x,y) based on a raster
    float getDemElevation(GeoLib::Point const& pnt) const;

    GeoLib::GEOObjects& _geo_objects;
    std::string& _geo_name;

    /// only necessary for mapping on mesh
    MeshLib::Mesh* _surface_mesh;
    GeoLib::Grid<MeshLib::Node>* _grid;

    /// only necessary for mapping on DEM
    std::unique_ptr<GeoLib::Raster const> _raster;
};

} // end namespace MeshGeoToolsLib
