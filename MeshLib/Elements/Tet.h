/**
 * \file
 * \author Thomas Fischer
 * \date   Sep 27, 2012
 * \brief  Definition of the Tet class.
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
#include "TetRule4.h"
#include "TetRule10.h"

extern template class MeshLib::TemplateElement<MeshLib::TetRule10>;
extern template class MeshLib::TemplateElement<MeshLib::TetRule4>;

namespace MeshLib {
using Tet = TemplateElement<MeshLib::TetRule4>;
using Tet10 = TemplateElement<MeshLib::TetRule10>;
}
