/**
 * \file
 * \date 2014-09-19
 * \brief Interface for heuristic search length strategy.
 *
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include "MeshGeoToolsLib/SearchLength.h"

namespace MeshLib
{
class Mesh;
}

namespace MeshGeoToolsLib
{

/// HeuristicSearchLength implements a mesh dependent criterion for searching
/// mesh nodes near a geometry. For this purpose it computes the average
/// \f$\mu\f$ and the standard deviation \f$\sigma\f$ of edge length or node distance of mesh
/// elements. The search length is set to \f$\mu-2\sigma\f$. This strategy
/// is useful for meshes with different sizes of elements.
class HeuristicSearchLength : public SearchLength
{
public:
    /// Type of length to be sampled
    enum class LengthType
    {
        Edge, /// edge length of elements, which is recommended for meshes without nonlinear nodes
        Node  /// distance between nodes
    };

    /**
     * Constructor
     * @param mesh  mesh object
     * @param length_type  length type to be sampled
     */
    explicit HeuristicSearchLength(MeshLib::Mesh const& mesh,
                                   LengthType length_type = LengthType::Edge);

private:
    MeshLib::Mesh const& _mesh;
};

} // end namespace MeshGeoToolsLib
