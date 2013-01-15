/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file SparsityBuilder.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "MathLib/LinAlg/Sparse/Sparsity.h"

namespace MeshLib
{
class Mesh;
}

namespace DiscreteLib
{
class DofEquationIdTable;

class SparsityBuilderDummy
{
public:
    SparsityBuilderDummy(const MeshLib::Mesh&, const DofEquationIdTable&, MathLib::RowMajorSparsity&) {};
};

}
