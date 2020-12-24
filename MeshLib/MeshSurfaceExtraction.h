/**
 * \file
 * \author Karsten Rink
 * \date   2013-04-04
 * \brief  Definition of the MeshSurfaceExtraction class
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <cstddef>
#include <vector>

#include "MathLib/Vector3.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Properties.h"

namespace MeshLib
{
// forward declarations
class Mesh;
class Element;
class Node;

/**
 * \brief A set of tools concerned with extracting nodes and elements from a
 * mesh surface
 */
class MeshSurfaceExtraction
{
public:
    /// Returns a vector of the areas assigned to each node on a surface mesh.
    static std::vector<double> getSurfaceAreaForNodes(
        const MeshLib::Mesh& mesh);

    /// Returns the surface nodes of a mesh.
    static std::vector<MeshLib::Node*> getSurfaceNodes(
        const MeshLib::Mesh& mesh, const MathLib::Vector3& dir, double angle);

    /**
     * Returns the 2d-element mesh representing the surface of the given mesh.
     * \param subsfc_mesh The original mesh
     * \param dir The direction in which face normals have to
     * point to be considered surface elements
     * \param angle The angle of the allowed deviation from the
     * given direction (0 <= angle <= 90 degrees)
     * \param subsfc_node_id_prop_name The name of the \c PropertyVector in
     * the surface mesh the subsurface mesh node ids are stored to. If the
     * string is empty, there won't be such a \c PropertyVector.
     * \param subsfc_element_id_prop_name The name of the PropertyVector in
     * the surface mesh that stores the subsurface element ids. If the string
     * is empty, there won't be such a \c PropertyVector.
     * \param face_id_prop_name The name of the \c PropertyVector in the surface
     * mesh that stores the face number of the subsurface element that belongs
     * to the boundary. If the string is empty, there won't be such a \c
     * PropertyVector.
     * \return A 2D mesh representing the surface in direction dir
     */
    static MeshLib::Mesh* getMeshSurface(
        const MeshLib::Mesh& subsfc_mesh,
        Eigen::Vector3d const& dir,
        double angle,
        std::string const& subsfc_node_id_prop_name = "",
        std::string const& subsfc_element_id_prop_name = "",
        std::string const& face_id_prop_name = "");

private:
    /// Functionality needed for getSurfaceNodes() and getMeshSurface()
    static void get2DSurfaceElements(
        const std::vector<MeshLib::Element*>& all_elements,
        std::vector<MeshLib::Element*>& sfc_elements,
        std::vector<std::size_t>& element_to_bulk_element_id_map,
        std::vector<std::size_t>& element_to_bulk_face_id_map,
        Eigen::Vector3d const& dir,
        double angle,
        unsigned mesh_dimension);
};

namespace BoundaryExtraction
{
std::unique_ptr<MeshLib::Mesh> getBoundaryElementsAsMesh(
    MeshLib::Mesh const& bulk_mesh,
    std::string const& subsfc_node_id_prop_name,
    std::string const& subsfc_element_id_prop_name,
    std::string const& face_id_prop_name);
}

void addBulkIDPropertiesToMesh(
    MeshLib::Mesh& surface_mesh,
    std::string const& node_to_bulk_node_id_map_name,
    std::vector<std::size_t> const& node_to_bulk_node_id_map,
    std::string const& element_to_bulk_element_id_map_name,
    std::vector<std::size_t> const& element_to_bulk_element_id_map,
    std::string const& element_to_bulk_face_id_map_name,
    std::vector<std::size_t> const& element_to_bulk_face_id_map);

}  // namespace MeshLib
