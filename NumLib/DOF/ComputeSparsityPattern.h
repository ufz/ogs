// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <vector>

#include "NumLib/NumericsConfig.h"

namespace MeshLib
{
class Mesh;
}

namespace NumLib
{
class LocalToGlobalIndexMap;

/**
 * \brief Computes a sparsity pattern for the given inputs.
 *
 * \param dof_table            maps mesh nodes to global indices
 * \param mesh                 mesh for which the two parameters above are
 * defined
 *
 * \return The computed sparsity pattern.
 */
GlobalSparsityPattern computeSparsityPattern(
    LocalToGlobalIndexMap const& dof_table, MeshLib::Mesh const& mesh);
}  // namespace NumLib
