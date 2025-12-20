// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "TemplateElement.h"
#include "LineRule2.h"
#include "LineRule3.h"

extern template class MeshLib::TemplateElement<MeshLib::LineRule2>;
extern template class MeshLib::TemplateElement<MeshLib::LineRule3>;

namespace MeshLib {
using Line = TemplateElement<MeshLib::LineRule2>;
using Line3 = TemplateElement<MeshLib::LineRule3>;
}
