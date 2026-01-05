// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
