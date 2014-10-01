/**
 * @date 2014-09-19
 * @brief Interface for heuristic search length strategy.
 *
 * @copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#ifndef HEURISTICSEARCHLENGTH_H_
#define HEURISTICSEARCHLENGTH_H_

#include "MeshGeoToolsLib/SearchLength.h"

namespace MeshGeoToolsLib
{

/// HeuristicSearchLength implements a mesh dependent criterion for searching
/// mesh nodes near a geometry. For this purpose it computes the average
/// \f$\mu\f$ and the standard deviation \f$\sigma\f$ of edge length of mesh
/// elements. The search lenght is set to \f$\mu-2\sigma\f$. This strategy
/// is usefull for meshes with different sizes of elements.
class HeuristicSearchLength : public SearchLength
{
public:
	explicit HeuristicSearchLength(MeshLib::Mesh const& mesh);
private:
	MeshLib::Mesh const& _mesh;
};

} // end namespace MeshGeoToolsLib

#endif

