/**
 * \file   ElementWiseLinearSystemUpdater.h
 * \author Norihiro Watanabe
 * \date   2012-08-03
 * \brief  Linear equation updater with element-wise operations
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef ELEMENTWISELINEARSYSTEMUPDATER_H_
#define ELEMENTWISELINEARSYSTEMUPDATER_H_

#include <vector>

#include "MeshLib/Mesh.h"
#include "MeshLib/Elements/Element.h"
#include "DiscreteLib/DoF/DofEquationIdTable.h"

namespace DiscreteLib
{

/**
 * \brief Element-wise global equation updater
 *
 * \tparam T_LOCAL Local assembler
 */
template <class T_LOCAL, class T_SOLVER>
class ElementWiseLinearSystemUpdater
{
public:
    typedef T_LOCAL LocalAssemblerType;
    typedef T_SOLVER SolverType;

    /**
     *
     * @param msh
     * @param local_assembler
     */
    ElementWiseLinearSystemUpdater(MeshLib::Mesh* msh, LocalAssemblerType* local_assembler)
    : _msh(msh), _e_assembler(local_assembler)
    {};

    /**
     *
     * @param e     mesh element
     * @param eqs   global equation
     */
    void update(const MeshLib::Element &e, const DofEquationIdTable &dofManager, SolverType &eqs);

private:
    MeshLib::Mesh* _msh;
    LocalAssemblerType* _e_assembler;
};


template <class T_LOCAL, class T_SOLVER>
void ElementWiseLinearSystemUpdater<T_LOCAL, T_SOLVER>::update(const MeshLib::Element &e, const DofEquationIdTable &dofManager, SolverType &eqs)
{
    std::vector<size_t> ele_node_ids, ele_node_size_order;
    std::vector<size_t> local_dofmap_row;
    std::vector<size_t> local_dofmap_column;

    // get dof map
    for (std::size_t i=0; i<e.getNNodes(); ++i)
        ele_node_ids.push_back(e.getNodeIndex(i));
    dofManager.mapEqsID(_msh->getID(), ele_node_ids, local_dofmap_row, local_dofmap_column);

    // local assembly
    LocalLinearSystem localEQS(local_dofmap_row.size());
    _e_assembler->assembly(e, localEQS);

    // update global
    eqs.addAsub(local_dofmap_row, local_dofmap_column, localEQS.getMat());
    eqs.addRHSsub(local_dofmap_row, localEQS.getRHSVec());
}
} //end

#endif //ELEMENTWISELINEAREQUATIONUPDATER_H_

