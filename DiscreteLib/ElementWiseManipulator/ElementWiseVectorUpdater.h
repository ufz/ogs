/**
 * \file   ElementWiseVectorUpdater.h
 * \author Norihiro Watanabe
 * \date   2012-08-03
 * \brief
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef ELEMENTWISEVECTORUPDATER_H_
#define ELEMENTWISEVECTORUPDATER_H_

#include <vector>
#include "MeshLib/Elements/Element.h"
#include "DiscreteLib/DoF/DofEquationIdTable.h"
#include "DiscreteLib/Core/IDiscreteVectorAssembler.h"

namespace DiscreteLib
{

/**
 * \brief Element-wise vector updater
 *
 * \tparam T_VALUE  vector data type
 * \tparam T_LOCAL  local assembler
 */
template <class T_VALUE, class T_LOCAL>
class ElementWiseVectorUpdater
{
public:
    typedef T_LOCAL LocalAssemblerType;
    typedef typename IDiscreteVectorAssembler<T_VALUE>::VectorType GlobalVectorType;

    /**
     * constructor
     * @param msh
     * @param a
     */
    ElementWiseVectorUpdater(std::size_t mesh_id, const LocalAssemblerType &a)
    : _msh_id(mesh_id), _e_assembler(a)
    {
    }

    /**
     * update a discrete vector for the given element
     * @param e             Mesh element
     * @param dofManager    DoF mapping table
     * @param globalVec     Discrete vector to be updated
     */
    void update(const MeshLib::Element &e, const DofEquationIdTable &dofManager, GlobalVectorType &globalVec)
    {
        const std::size_t e_nnodes = e.getNNodes();
        LocalVector localVec;
        std::vector<std::size_t> ele_node_ids(e_nnodes);
        std::vector<std::size_t> local_dofmap_row;

        std::vector<T_VALUE> local_u_n;
        // get dof map
        for (std::size_t i=0; i<e_nnodes; ++i)
            ele_node_ids[i] = e.getNodeIndex(i);
        dofManager.mapEqsID(_msh_id, ele_node_ids, local_dofmap_row);
        // local assembly
        localVec.resize(local_dofmap_row.size(), .0);
        _e_assembler.assembly(e, localVec);
        // update global
        globalVec.addSubvector(local_dofmap_row, &localVec[0]);
    }

private:
    const std::size_t _msh_id;
    LocalAssemblerType _e_assembler;
};

} //end

#endif //ELEMENTWISEVECTORUPDATER_H_
