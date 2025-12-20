// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
