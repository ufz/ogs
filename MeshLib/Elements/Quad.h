/**
 * \file
 * \author Thomas Fischer
 * \date   Sep 27, 2012
 * \brief  Definition of the Quad class.
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef QUAD_H_
#define QUAD_H_

#include "TemplateQuad.h"

namespace MeshLib
{

namespace detail
{
class QuadEdgeLinearNodes
{
protected:
	static constexpr unsigned _edge_nodes[4][2] = {
		{0, 1}, // Edge 0
		{1, 2}, // Edge 1
		{2, 3}, // Edge 2
		{0, 3}  // Edge 3
	};
};

class QuadEdgeQuadraticNodes
{
protected:
	static constexpr unsigned _edge_nodes[4][3] = {
		{0, 1, 4}, // Edge 0
		{1, 2, 5}, // Edge 1
		{2, 3, 6}, // Edge 2
		{0, 3, 7}  // Edge 3
	};
};
} // end detail

typedef TemplateQuad<4, CellType::QUAD4, detail::QuadEdgeLinearNodes> Quad;
typedef TemplateQuad<8, CellType::QUAD8, detail::QuadEdgeQuadraticNodes> Quad8;
typedef TemplateQuad<9, CellType::QUAD9, detail::QuadEdgeQuadraticNodes> Quad9;

}

#endif /* QUAD_H_ */
