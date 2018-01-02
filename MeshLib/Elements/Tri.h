/**
 * \file
 * \author Thomas Fischer
 * \date   Sep 27, 2012
 * \brief  Definition of the Tri class.
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
#include "TriRule3.h"
#include "TriRule6.h"


extern template class MeshLib::TemplateElement<MeshLib::TriRule3>;
extern template class MeshLib::TemplateElement<MeshLib::TriRule6>;

namespace MeshLib {
using Tri = TemplateElement<MeshLib::TriRule3>;
using Tri6 = TemplateElement<MeshLib::TriRule6>;
}
