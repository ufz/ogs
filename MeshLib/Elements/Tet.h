// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
