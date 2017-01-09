/**
 * \file
 * \author Thomas Fischer
 * \date   Sep 27, 2012
 * \brief  Definition of the Hex class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef HEX_H_
#define HEX_H_

#include "TemplateElement.h"
#include "HexRule8.h"
#include "HexRule20.h"

extern template class MeshLib::TemplateElement<MeshLib::HexRule20>;
extern template class MeshLib::TemplateElement<MeshLib::HexRule8>;

namespace MeshLib {
typedef TemplateElement<HexRule8> Hex;
typedef TemplateElement<HexRule20> Hex20;
}

#endif /* HEX_H_ */
