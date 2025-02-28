/**
 * \file
 * \author Thomas Fischer
 * \date   Oct 12, 2012
 * \brief  Implementation of the Mesh2MeshPropertyInterpolation class.
 *
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "MeshLib/PropertyVector.h"

namespace MeshLib
{

class Mesh;
}

namespace MeshToolsLib
{

/**
 * Class Mesh2MeshPropertyInterpolation transfers properties of
 * mesh elements of a (source) mesh to mesh elements of another
 * (destination) mesh deploying weighted interpolation. The two
 * meshes must have the same dimension.
 */
class Mesh2MeshPropertyInterpolation final
{
public:
    /**
     * Constructor taking the source or input mesh and properties.
     * @param src_mesh the mesh the given property information is assigned to.
     * @param property_name is the name of a PropertyVector in the \c
     * source_mesh
     */
    Mesh2MeshPropertyInterpolation(MeshLib::Mesh const& src_mesh,
                                   std::string const& property_name);

    /**
     * Calculates entries for the property vector and sets appropriate indices
     * in the mesh elements.
     * @param dest_mesh the mesh the property information will be calculated and
     * set via weighted interpolation
     * @return true if the operation was successful, false on error
     */
    bool setPropertiesForMesh(MeshLib::Mesh& dest_mesh) const;

private:
    /**
     * @param dest_mesh
     * @param dest_properties
     */
    void interpolatePropertiesForMesh(
        MeshLib::Mesh const& dest_mesh,
        MeshLib::PropertyVector<double>& dest_properties) const;

    /**
     * Method interpolates the element wise given properties to the nodes of the
     * element
     * @param interpolated_properties the vector must have the same number of
     * entries as the source mesh has number of nodes, the content of the
     * particular entries will be overwritten
     */
    void interpolateElementPropertiesToNodeProperties(
        std::vector<double>& interpolated_properties) const;

    MeshLib::Mesh const& _src_mesh;
    std::string const& _property_name;
};

}  // namespace MeshToolsLib
