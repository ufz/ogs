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

#ifndef MATRIXASSEMBLER_H_
#define MATRIXASSEMBLER_H_

#include <vector>

#include "ITask.h"

namespace VecMatOnMeshLib
{


template<class T_MAT, class T_MESH_ITEM, class T_LOCAL_ASSEMBLY>
class MatrixAssembler : public ITask<T_MESH_ITEM>
{
public:
    MatrixAssembler(T_MAT &mat, T_LOCAL_ASSEMBLY &local_assembler, const std::vector<std::vector<std::size_t> > &data_pos)
    : _mat(mat), _local_assembler(local_assembler), _data_pos(data_pos) {}

    virtual ~MatrixAssembler() {}

    virtual void operator()(const T_MESH_ITEM* item, std::size_t id);

protected:
    T_MAT &_mat;
    T_LOCAL_ASSEMBLY &_local_assembler;
    const std::vector<std::vector<std::size_t> > &_data_pos;
};


}

#include "MatrixAssembler.tpp"

#endif /* MATRIXASSEMBLER_H_ */
