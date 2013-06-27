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

#ifndef MATRIXASSEMBLER_TPP_
#define MATRIXASSEMBLER_TPP_

#include "MathLib/LinAlg/Dense/DenseMatrix.h"

namespace VecMatOnMeshLib
{

template<class T_MAT, class T_MESH_ITEM, class T_LOCAL_ASSEMBLY>
void MatrixAssembler<T_MAT,T_MESH_ITEM,T_LOCAL_ASSEMBLY>::operator()(const T_MESH_ITEM* item, std::size_t id)
{
    assert(_data_pos.size() > id);

    auto pos = _data_pos[id];
    MathLib::DenseMatrix<double> local_mat(pos.size(), pos.size());
    _local_assembler(*item, local_mat);
    _mat.addSubMatrix(pos, pos, local_mat);
}

} // VecMatOnMeshLib

#endif /* MATRIXASSEMBLER_TPP_ */
