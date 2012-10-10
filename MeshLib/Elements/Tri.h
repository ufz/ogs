/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 *
 * \file Tri.h
 *
 *  Created on  Sep 27, 2012 by Thomas Fischer
 */

#ifndef TRI_H_
#define TRI_H_

#include "TemplateTri.h"

namespace MeshLib {

typedef TemplateTri<3,FEMElemType::TRI3> Tri;

}

#endif /* TRI_H_ */
