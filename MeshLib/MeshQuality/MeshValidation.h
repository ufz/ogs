/**
 * \file
 * \author Karsten Rink
 * \date   2013-04-04
 * \brief  Definition of the MeshValidation class
 *
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <array>
#include <vector>
#include <limits>

#include "MeshLib/Elements/ElementErrorCode.h"

namespace MeshLib {
    class Mesh;
    class Element;

/**
 * \brief A collection of methods for testing mesh quality and correctness
 */
struct MeshValidation final
{
    /// Constructor
    /// \warning This might change the mesh when removing unused mesh nodes.
    explicit MeshValidation(MeshLib::Mesh& mesh);

    /**
     * Tests if elements are geometrically correct.
     * @param mesh The mesh that is tested
     * @param min_volume The minimum required volume for a mesh element, so it is NOT considered faulty
     * @return Vector of error codes for each mesh element
     */
    static std::vector<ElementErrorCode> testElementGeometry(
        const MeshLib::Mesh &mesh,
        double min_volume = std::numeric_limits<double>::epsilon());

    /**
     * Detailed output which ElementID is associated with which error(s)
     * @return String containing the report
     */
    static std::array<std::string, static_cast<std::size_t>(ElementErrorFlag::MaxValue)>
        ElementErrorCodeOutput(const std::vector<ElementErrorCode> &error_codes);

    /**
     * Tests if holes are located within the mesh.
     * In this context, a hole is a boundary of an element with no neighbor that cannot be reached from
     * the actual boundary of the mesh. Examples include a missing triangle in a 2D mesh or a missing
     * Tetrahedron in a 3D mesh. The method does not work for 1d-meshes.
     * Note, that this method does not work when complex 3D structures are build from 2D mesh elements,
     * e.g. using the LayeredVolume-class, where more than two 2D elements may share an edge.
     * @param mesh The mesh that is tested
     * @return The number of holes that have been found.
     */
    static unsigned detectHoles(MeshLib::Mesh const& mesh);
};

} // end namespace MeshLib
