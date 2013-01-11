/**
 * \file
 * \author Thomas Fischer
 * \date   Sep 27, 2012
 * \brief  Definition of the Hex class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef HEX_H_
#define HEX_H_

#include "TemplateHex.h"

namespace MeshLib {
typedef TemplateHex<8, CellType::HEX8> Hex;
}

#endif /* HEX_H_ */
