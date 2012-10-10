/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 *
 * \file Hex.h
 *
 *  Created on  Sep 27, 2012 by Thomas Fischer
 */

#ifndef HEX_H_
#define HEX_H_

#include "TemplateHex.h"

namespace MeshLib {
typedef TemplateHex<8, FEMElemType::HEX8> Hex;
}

#endif /* HEX_H_ */
