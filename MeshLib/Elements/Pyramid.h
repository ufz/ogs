/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 *
 * \file Pyramid.h
 *
 *  Created on  Sep 27, 2012 by Thomas Fischer
 */

#ifndef PYRAMID_H_
#define PYRAMID_H_

#include "TemplatePyramid.h"

namespace MeshLib {

typedef TemplatePyramid<5,CellType::PYRAMID5> Pyramid;

}

#endif /* PYRAMID_H_ */
