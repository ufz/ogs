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
#include "TemplateQuadQuadratic.h"

namespace MeshLib
{

typedef TemplateQuad<4, CellType::QUAD4> Quad;
typedef TemplateQuadQuadratic<8, CellType::QUAD8> Quad8;
typedef TemplateQuadQuadratic<9, CellType::QUAD9> Quad9;

}

#endif /* QUAD_H_ */
