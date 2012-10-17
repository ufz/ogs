/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 *
 * \file Prism.h
 *
 *  Created on  Sep 27, 2012 by Thomas Fischer
 */

#ifndef PRISM_H_
#define PRISM_H_

#include "TemplatePrism.h"

namespace MeshLib {

typedef TemplatePrism<6, CellType::PRISM6> Prism;

}

#endif /* PRISM_H_ */
