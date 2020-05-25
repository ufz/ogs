/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */
#pragma once

#include <vector>

namespace GeoLib
{
class Surface;
}

namespace MeshLib
{
class Mesh;
class Element;
}

namespace MeshGeoToolsLib
{
class MeshNodeSearcher;

/**
 * This class collects element faces located along a surface.
 * Note that internal faces are not collected in this class.
 */
class BoundaryElementsOnSurface
{
public:
    /**
     * Constructor
     * @param mesh             a mesh object
     * @param mshNodeSearcher  a MeshNodeSearcher object which is internally used to search mesh nodes
     * @param sfc              a surface object where face elements are searched for
     */
    BoundaryElementsOnSurface(MeshLib::Mesh const& mesh,
                              MeshNodeSearcher const& mshNodeSearcher,
                              GeoLib::Surface const& sfc);

    /// destructor
    virtual ~BoundaryElementsOnSurface();

    /// return the mesh object
    MeshLib::Mesh const& getMesh() const {return mesh_;}

    /**
     * Deploying this method the user can get access to the underlying
     * GeoLib::Surface.
     * @return the underlying GeoLib::Surface
     */
    GeoLib::Surface const& getSurface() const {return sfc_;}

    /**
     * Return the vector of boundary elements (i.e. faces). The elements are unsorted.
     */
    std::vector<MeshLib::Element*> const& getBoundaryElements() const {return boundary_elements_;}

private:
    MeshLib::Mesh const& mesh_;
    GeoLib::Surface const& sfc_;
    std::vector<MeshLib::Element*> boundary_elements_;
};

} // end namespace MeshGeoToolsLib
