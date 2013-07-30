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

#ifndef VECTORASSEMBLER_H_
#define VECTORASSEMBLER_H_

#include <vector>

#include "LocalToGlobalIndexMap.h"

namespace VecMatOnMeshLib
{

/**
 * Get a local linear algebra object and assemble it into the global one.
 *
 * \tparam LINALG_OBJ_T     Linear algebra object type
 * \tparam T_MESH_ITEM      Mesh item type in MeshItemType.h
 * \tparam T_LOCAL_ASSEMBLY Local assembler class
 */
template<class LINALG_OBJ_T, class T_MESH_ITEM, class T_LOCAL_ASSEMBLY>
class VectorAssembler
{
public:
	/**
	 * @param vec               Global linear algebra object
	 * @param local_assembler   Local assembler object
	 * @param data_pos          Mapping table from DoFs to position in the global linear algebra object.
	 * The number of rows should equal to the number of mesh items. The number
	 * of columns should equal to the number of components on that mesh item.
	 */
    VectorAssembler(LINALG_OBJ_T &linalg_obj,
		T_LOCAL_ASSEMBLY &local_assembler,
		LocalToGlobalIndexMap const& data_pos)
    : _linalg_obj(linalg_obj), _local_assembler(local_assembler), _data_pos(data_pos) {}

    virtual ~VectorAssembler() {}

    /**
     * do some task for the given mesh item and update the global linear algebra
     * object
     *
     * @param item  Pointer to a mesh item
     * @param id    Index of the mesh item. The index is used to search a mapping in data_pos
     */
    virtual void operator()(const T_MESH_ITEM* item, std::size_t id)
	{
		assert(_data_pos.size() > id);

		LocalToGlobalIndexMap::RowColumnIndices const& indices = _data_pos[id];
		MathLib::DenseVector<double> local_linalg_obj(indices.rows.size());
		_local_assembler(*item, local_linalg_obj);
		_linalg_obj.add(indices.rows, local_linalg_obj);
	}

protected:
    LINALG_OBJ_T &_linalg_obj;
    T_LOCAL_ASSEMBLY &_local_assembler;
    LocalToGlobalIndexMap const& _data_pos;
};

} // VecMatOnMeshLib

#endif /* VECTORASSEMBLER_H_ */
