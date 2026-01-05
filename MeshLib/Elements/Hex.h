// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
