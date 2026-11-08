// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <cstddef>
#include <vector>

namespace MeshLib
{
class Element;
}

namespace MeshToolsLib
{
/**
 * \brief Reverses order of nodes. In particular, this fixes issues between
 * (Gmsh or OGS5) and OGS6 meshes.
 *
 * The elements are modified in place: only those whose node order is wrong
 * according to the OGS6 element definitions are touched.
 *
 * \param elements        Mesh elements whose nodes should be reordered
 * \param mesh_space_dimension Dimension of the embedding space
 * \param forced          If true, additionally verify the Jacobian determinant
 *                        of every reordered element and abort if it is still
 *                        negative. Whether an element is reordered at all is
 *                        decided by its node order in either case.
 * \param volume_check    If true, verify that reordering preserved the element
 *                        volume within 100 machine epsilons, absolute or
 *                        relative to the original volume.
 * \return The number of corrected elements.
 * \note Terminates via OGS_FATAL (which raises std::runtime_error) for
 *       unsupported element types, when the volume check fails, and, if
 *       \p forced is set, when the Jacobian determinant of a reordered element
 *       is negative.
 */
std::size_t reverseNodeOrder(std::vector<MeshLib::Element*> const& elements,
                             int mesh_space_dimension, bool forced,
                             bool volume_check);
}  // namespace MeshToolsLib
