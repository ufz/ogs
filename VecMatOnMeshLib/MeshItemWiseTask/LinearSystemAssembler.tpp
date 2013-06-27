/**
 * \file
 * \author Norihiro Watanabe
 * \date   2013-04-16
 * \brief
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef LINEARSYSTEMASSEMBLER_TPP_
#define LINEARSYSTEMASSEMBLER_TPP_

#include "MathLib/LinAlg/Dense/DenseMatrix.h"
#include "MathLib/LinAlg/Dense/DenseVector.h"

namespace VecMatOnMeshLib
{

template<class T_MAT, class T_VEC, class T_MESH_ITEM, class T_LOCAL_ASSEMBLY>
void LinearSystemAssembler<T_MAT, T_VEC, T_MESH_ITEM, T_LOCAL_ASSEMBLY>::operator()(const T_MESH_ITEM* item, std::size_t id)
{
    assert(_data_pos.size() > id);

    auto pos = _data_pos[id];
    MathLib::DenseMatrix<double> local_A(pos.size(), pos.size());
    MathLib::DenseVector<double> local_rhs(pos.size());
    _local_assembler(*item, local_A, local_rhs);
    _A.addSubMatrix(pos, pos, local_A);
    _rhs.addSubVector(pos, local_rhs);
}

} //end VecMatOnMeshLib

#endif /* LINEARSYSTEMASSEMBLER_TPP_ */
