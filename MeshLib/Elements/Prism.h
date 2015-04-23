/**
 * \file
 * \author Thomas Fischer
 * \date   Sep 27, 2012
 * \brief  Definition of the Prism class.
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PRISM_H_
#define PRISM_H_

#include "TemplateElement.h"
#include "PrismRule6.h"
#include "PrismRule15.h"

namespace MeshLib {

typedef TemplateElement<PrismRule6> Prism;
typedef TemplateElement<PrismRule15> Prism15;

}

#endif /* PRISM_H_ */
