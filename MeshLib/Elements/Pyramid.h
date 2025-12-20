// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "TemplateElement.h"
#include "PyramidRule5.h"
#include "PyramidRule13.h"

extern template class MeshLib::TemplateElement<MeshLib::PyramidRule13>;
extern template class MeshLib::TemplateElement<MeshLib::PyramidRule5>;

namespace MeshLib {
using Pyramid = TemplateElement<MeshLib::PyramidRule5>;
using Pyramid13 = TemplateElement<MeshLib::PyramidRule13>;
}
