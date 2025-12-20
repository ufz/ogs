// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
