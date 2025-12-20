// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "TemplateElement.h"
#include "PointRule1.h"

extern template class MeshLib::TemplateElement<MeshLib::PointRule1>;

namespace MeshLib
{
    using Point = TemplateElement<PointRule1>;
}
