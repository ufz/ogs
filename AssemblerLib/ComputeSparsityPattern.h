/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef ASSEMBLERLIB_COMPUTESPARSITYPATTERN_H
#define ASSEMBLERLIB_COMPUTESPARSITYPATTERN_H

#include <vector>

#include "MathLib/LinAlg/SparsityPattern.h"
#include "ProcessLib/NumericsConfig.h"

namespace MeshLib
{
class Mesh;
}

namespace AssemblerLib
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
MathLib::SparsityPattern computeSparsityPattern(
    LocalToGlobalIndexMap const& dof_table, MeshLib::Mesh const& mesh);
}

#endif // ASSEMBLERLIB_COMPUTESPARSITYPATTERN_H

