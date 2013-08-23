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

#include "LocalToGlobalIndexMap.h"

namespace VecMatOnMeshLib
{

template<class T_MAT, class T_VEC, class T_MESH_ITEM, class T_LOCAL_ASSEMBLY>
class LinearSystemAssembler
{
public:
    LinearSystemAssembler(T_MAT &A, T_VEC &rhs,
        T_LOCAL_ASSEMBLY &local_assembler,
        LocalToGlobalIndexMap const& data_pos)
    : _A(A), _rhs(rhs), _local_assembler(local_assembler), _data_pos(data_pos) {}

    virtual ~LinearSystemAssembler() {}

    virtual void operator()(const T_MESH_ITEM* item, std::size_t id) const
    {
        assert(_data_pos.size() > id);

        LocalToGlobalIndexMap::RowColumnIndices const& indices = _data_pos[id];
        MathLib::DenseMatrix<double> local_A(indices.rows.size(), indices.columns.size());
        MathLib::DenseVector<double> local_rhs(indices.rows.size());
        _local_assembler(*item, local_A, local_rhs);
        _A.add(indices, local_A);
        _rhs.add(indices.rows, local_rhs);
    }

protected:
    T_MAT &_A;
    T_VEC &_rhs;
    T_LOCAL_ASSEMBLY &_local_assembler;
    LocalToGlobalIndexMap const& _data_pos;
};

} // VecMatOnMeshLib

#endif /* LINEARSYSTEMASSEMBLER_H_ */
