// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <range/v3/algorithm/max.hpp>

#include "MeshLib/PropertyVector.h"

namespace MeshLib
{
/// Returns one past the largest material ID present in \c material_ids, or zero
/// if there are no material IDs. Gaps in the used IDs are not filled.
inline int nextUnusedMaterialId(PropertyVector<int> const& material_ids)
{
    return material_ids.empty() ? 0 : ranges::max(material_ids) + 1;
}
}  // namespace MeshLib
