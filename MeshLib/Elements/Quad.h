/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 *
 * \file Quad.h
 *
 *  Created on  Sep 27, 2012 by Thomas Fischer
 */

#ifndef QUAD_H_
#define QUAD_H_

#include "TemplateQuad.h"

namespace MeshLib {

typedef TemplateQuad<4,FEMElemType::QUAD4> Quad;

}

#endif /* QUAD_H_ */
