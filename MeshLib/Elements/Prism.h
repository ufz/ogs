/**
 * \file
 * \author Thomas Fischer
 * \date   Sep 27, 2012
 * \brief  Definition of the Prism class.
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "TemplateElement.h"
#include "PrismRule6.h"
#include "PrismRule15.h"

extern template class MeshLib::TemplateElement<MeshLib::PrismRule15>;
extern template class MeshLib::TemplateElement<MeshLib::PrismRule6>;

namespace MeshLib {
using Prism = TemplateElement<MeshLib::PrismRule6>;
using Prism15 = TemplateElement<MeshLib::PrismRule15>;
}
