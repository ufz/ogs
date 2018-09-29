/**
 * \file   DuplicateMeshComponents.h
 * \author Karsten Rink
 * \date   2014-03-25
 * \brief  Definition of Duplicate functions
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <vector>

namespace MeshLib
{

    class Mesh;
    class Node;
    class Element;

    /// Creates a deep copy of a Node vector
    std::vector<MeshLib::Node*> copyNodeVector(const std::vector<MeshLib::Node*> &nodes);

    /** Creates a deep copy of an element vector using the given Node vector.
     * @param elements The element vector that should be duplicated.
     * @param nodes    The new node vector used for the duplicated element vector. This should be consistent with the original node vector.
     * @return A deep copy of the elements vector using the new nodes vector.
     */
    std::vector<MeshLib::Element*> copyElementVector(const std::vector<MeshLib::Element*> &elements, const std::vector<MeshLib::Node*> &nodes);

    /// Copies an element without change, using the nodes vector from the result mesh.
    MeshLib::Element* copyElement(MeshLib::Element const*const element,
                                  const std::vector<MeshLib::Node*> &nodes);

    /// Copies an element without change, using the nodes vector from the result mesh.
    template <typename E>
    MeshLib::Element* copyElement(MeshLib::Element const*const element, const std::vector<MeshLib::Node*> &nodes);

    /// Clones a vector of elements using the Element::clone() function.
    std::vector<MeshLib::Element*> cloneElements(
        std::vector<MeshLib::Element*> const& elements);
} // end namespace MeshLib
