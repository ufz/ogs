/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef NUMLIB_COMPUTESPARSITYPATTERN_H
#define NUMLIB_COMPUTESPARSITYPATTERN_H

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
 * @brief Computes a sparsity pattern for the given inputs.
 *
 * @param dof_table            maps mesh nodes to global indices
 * @param mesh                 mesh for which the two parameters above are defined
 *
 * @return The computed sparsity pattern.
 */
GlobalSparsityPattern computeSparsityPattern(
    LocalToGlobalIndexMap const& dof_table, MeshLib::Mesh const& mesh);
}

#endif // NUMLIB_COMPUTESPARSITYPATTERN_H

