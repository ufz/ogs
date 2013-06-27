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

#ifndef VECTORASSEMBLER_TPP_
#define VECTORASSEMBLER_TPP_

#include <cassert>

#include "MathLib/LinAlg/Dense/DenseVector.h"

namespace VecMatOnMeshLib
{

template<class T_VEC, class T_MESH_ITEM, class T_LOCAL_ASSEMBLY>
VectorAssembler<T_VEC, T_MESH_ITEM, T_LOCAL_ASSEMBLY>::VectorAssembler(T_VEC &vec, T_LOCAL_ASSEMBLY &local_assembler, const std::vector<std::vector<std::size_t> > &data_pos )
: _global_vec(vec), _local_assembler(local_assembler), _data_pos(data_pos) {}

template<class T_VEC, class T_MESH_ITEM, class T_LOCAL_ASSEMBLY>
void VectorAssembler<T_VEC, T_MESH_ITEM, T_LOCAL_ASSEMBLY>::operator()(const T_MESH_ITEM* item, std::size_t id)
{
    assert(_data_pos.size() > id);

    auto vec_pos = _data_pos[id];

    MathLib::DenseVector<double> local_vec(vec_pos.size());
    _local_assembler(*item, local_vec);
    _global_vec.addSubVector(vec_pos, local_vec);
}

} // VecMatOnMeshLib

#endif /* VECTORASSEMBLER_TPP_ */
