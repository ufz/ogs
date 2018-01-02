/**
 * \file   MeshLayerMapper.h
 * \author Karsten Rink
 * \date   2010-11-01
 * \brief  Definition of the MeshLayerMapper class.
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "LayeredMeshGenerator.h"

namespace MeshLib
{

/**
 * \brief Manipulating and adding prism element layers to an existing 2D mesh
 */
class MeshLayerMapper : public LayeredMeshGenerator
{
public:
    ~MeshLayerMapper() override = default;

    /**
    * Based on a 2D triangle-or quad mesh this method creates a 3D mesh with a given number of prism- or hex-layers
    * \param mesh The triangle/quad mesh that is the basis for the new prism/hex mesh
    * \param layer_thickness_vector The size of the vector equals the number of
        layers of prism/hex elements that will be extruded from the
        triangle/quad elements of the original mesh. The thickness of the
        \f$i\f$-th layer corresponds to the \f$i\f$ entry of the vector.
    * \param mesh_name The name of the newly created mesh
    * \return A mesh with the requested number of layers of prism/hex elements
    */
    static MeshLib::Mesh* createStaticLayers(MeshLib::Mesh const& mesh, std::vector<float> const& layer_thickness_vector, std::string const& mesh_name = "SubsurfaceMesh");

    /**
    * Based on a 2D triangle mesh this method creates a 3D mesh with a given number of prism-layers.
    * Note: While this method would technically also work with quad meshes, this is discouraged as quad elements will most likely not
    * be coplanar after the mapping process which result in invaled mesh elements.
    * \param mesh                    The 2D triangle mesh that is the basis for the new 3D prism mesh
    * \param rasters                 Containing all the raster-data for the subsurface layers from bottom to top (starting with the bottom of the oldest layer and ending with the DEM)
    * \param minimum_thickness       Minimum thickness of each of the newly created layers (i.e. nodes with a vertical distance smaller than this will be collapsed)
    * \param noDataReplacementValue  Default z-coordinate if there are mesh nodes not located on the DEM raster (i.e. raster_paths[0])
    * \return A mesh with the requested number of layers of prism elements (also including Tet- & Pyramid-elements in case of degenerated prisms)
    */
    bool createRasterLayers(MeshLib::Mesh const& mesh,
                            std::vector<GeoLib::Raster const*> const& rasters,
                            double minimum_thickness,
                            double noDataReplacementValue = 0.0) override;

    /**
    * Maps the elevation of nodes of a given 2D mesh according to the raster. At
    * locations where no information is given, node elevation is set to
    * noDataReplacementValue.
    */
    static bool layerMapping(MeshLib::Mesh &mesh, const GeoLib::Raster &raster, double noDataReplacementValue);

    /// Maps the elevation of all mesh nodes to the specified static value.
    static bool mapToStaticValue(MeshLib::Mesh &mesh, double value);

private:
    /// Adds another layer to a subsurface mesh
    void addLayerToMesh(const MeshLib::Mesh& mesh_layer,
                        unsigned layer_id,
                        GeoLib::Raster const& raster) override;
};

} // end namespace MeshLib
