/**
 * \file
 * \author Thomas Fischer
 * \date   Sep 27, 2012
 * \brief  Definition of the Quad class.
 *
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "TemplateElement.h"
#include "QuadRule4.h"
#include "QuadRule8.h"
#include "QuadRule9.h"

extern template class MeshLib::TemplateElement<MeshLib::QuadRule4>;
extern template class MeshLib::TemplateElement<MeshLib::QuadRule8>;
extern template class MeshLib::TemplateElement<MeshLib::QuadRule9>;

namespace MeshLib
{
using Quad = TemplateElement<MeshLib::QuadRule4>;
using Quad8 = TemplateElement<MeshLib::QuadRule8>;
using Quad9 = TemplateElement<MeshLib::QuadRule9>;
}
