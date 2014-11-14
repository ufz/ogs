/**
 * \file
 * \author Thomas Fischer
 * \date   Sep 27, 2012
 * \brief  Definition of the Hex class.
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef HEX_H_
#define HEX_H_

//#include "TemplateHex.h"
#include "TemplateElement.h"
#include "HexRule8.h"
#include "HexRule20.h"
#include "Cell.h"

namespace MeshLib
{

typedef TemplateElement<Cell, HexRule8> Hex;
typedef TemplateElement<Cell, HexRule20> Hex20;
//typedef TemplateHex<8, CellType::HEX8> Hex;

}

#endif /* HEX_H_ */
