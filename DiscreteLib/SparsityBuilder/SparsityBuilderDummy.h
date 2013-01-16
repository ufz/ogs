/**
 * \file
 * \author Norihiro Watanabe
 * \date   2012-08-03
 * \brief  Helper macros.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef SPARSITYBUILDERDUMMY_H_
#define SPARSITYBUILDERDUMMY_H_

#include "MathLib/LinAlg/Sparse/Sparsity.h"

namespace MeshLib
{
class Mesh;
}

namespace DiscreteLib
{
class DofEquationIdTable;

/**
 *
 */
class SparsityBuilderDummy
{
public:
    SparsityBuilderDummy(const MeshLib::Mesh&, const DofEquationIdTable&, MathLib::RowMajorSparsity&) {};
};

}

#endif //SPARSITYBUILDERDUMMY_H_
