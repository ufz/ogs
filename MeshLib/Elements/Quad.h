/**
 * \file
 * \author Thomas Fischer
 * \date   Sep 27, 2012
 * \brief  Definition of the Quad class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef QUAD_H_
#define QUAD_H_

#include "TemplateQuad.h"

namespace MeshLib {

typedef TemplateQuad<4,CellType::QUAD4> Quad;
typedef TemplateQuad<8, CellType::QUAD8> Quad8;
typedef TemplateQuad<9, CellType::QUAD9> Quad9;

}

#endif /* QUAD_H_ */
