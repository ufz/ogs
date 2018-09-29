/**
 * \file
 * \author Thomas Fischer
 * \date   Sep 27, 2012
 * \brief  Definition of the Hex class.
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
#include "HexRule8.h"
#include "HexRule20.h"

extern template class MeshLib::TemplateElement<MeshLib::HexRule20>;
extern template class MeshLib::TemplateElement<MeshLib::HexRule8>;

namespace MeshLib {
using Hex = TemplateElement<MeshLib::HexRule8>;
using Hex20 = TemplateElement<MeshLib::HexRule20>;
}
