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
#include "ITask.h"

namespace VecMatOnMeshLib
{

/**
 * Get a local vector and assemble it into a global vector
 *
 * \tparam T_VEC				Vector class
 * \tparam T_MESH_ITEM			Mesh item type in MeshItemType.h
 * \tparam T_LOCAL_ASSEMBLY		Local assembler class
 */
template<class T_VEC, class T_MESH_ITEM, class T_LOCAL_ASSEMBLY>
class VectorAssembler : public ITask<T_MESH_ITEM>
{
public:
	/**
	 * constructor
	 *
	 * @param vec               Global vector object
	 * @param local_assembler   Local assembler object
	 * @param data_pos          Mapping table from DoFs to position in the global vector.
	 * The number of rows should equal to the number of mesh items. The number
	 * of columns should equal to the number of components on that mesh item.
	 */
    VectorAssembler(T_VEC &vec, T_LOCAL_ASSEMBLY &local_assembler, const std::vector<std::vector<std::size_t> > &data_pos );

    virtual ~VectorAssembler() {}

    /**
     * do some task for the given mesh item and update the global vector
     *
     * @param item  Pointer to a mesh item
     * @param id    Index of the mesh item. The index is used to search a mapping in data_pos
     */
    virtual void operator()(const T_MESH_ITEM* item, std::size_t id);

protected:
    T_VEC &_global_vec;
    T_LOCAL_ASSEMBLY &_local_assembler;
    const std::vector<std::vector<std::size_t> > &_data_pos;
};

} // VecMatOnMeshLib

#include "VectorAssembler.tpp"

#endif /* VECTORASSEMBLER_H_ */
