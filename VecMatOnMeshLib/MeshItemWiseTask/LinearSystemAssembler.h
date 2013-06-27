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

#ifndef LINEARSYSTEMASSEMBLER_H_
#define LINEARSYSTEMASSEMBLER_H_

#include <vector>

#include "ITask.h"

namespace VecMatOnMeshLib
{

template<class T_MAT, class T_VEC, class T_MESH_ITEM, class T_LOCAL_ASSEMBLY>
class LinearSystemAssembler : public ITask<T_MESH_ITEM>
{
public:
    LinearSystemAssembler(T_MAT &A, T_VEC &rhs, T_LOCAL_ASSEMBLY &local_assembler, const std::vector<std::vector<std::size_t> > &data_pos)
    : _A(A), _rhs(rhs), _local_assembler(local_assembler), _data_pos(data_pos) {}

    virtual ~LinearSystemAssembler() {}

    virtual void operator()(const T_MESH_ITEM* item, std::size_t id);

protected:
    T_MAT &_A;
    T_VEC &_rhs;
    T_LOCAL_ASSEMBLY &_local_assembler;
    const std::vector<std::vector<std::size_t> > &_data_pos;
};

} // VecMatOnMeshLib

#include "LinearSystemAssembler.tpp"

#endif /* LINEARSYSTEMASSEMBLER_H_ */
