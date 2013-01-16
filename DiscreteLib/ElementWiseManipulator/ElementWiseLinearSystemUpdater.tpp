/**
 * \file   ElementWiseLinearSystemUpdater.tpp
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

#ifndef ELEMENTWISELINEARSYSTEMUPDATER_TPP_
#define ELEMENTWISELINEARSYSTEMUPDATER_TPP_

namespace DiscreteLib
{

template <class T_LOCAL, class T_SOLVER>
void ElementWiseLinearSystemUpdater<T_LOCAL, T_SOLVER>::update(const MeshLib::Element &e, const DofEquationIdTable &dofManager, SolverType &eqs)
{
    const std::size_t e_nnodes = e.getNNodes();
    std::vector<std::size_t> ele_node_ids(e_nnodes);
    std::vector<std::size_t> local_dofmap_row;
    std::vector<std::size_t> local_dofmap_column;

    // get dof map
    for (std::size_t i=0; i<e_nnodes; ++i)
        ele_node_ids[i] = e.getNodeIndex(i);
    dofManager.mapEqsID(_msh_id, ele_node_ids, local_dofmap_row, local_dofmap_column);

    // local assembly
    LocalLinearSystem localEQS(local_dofmap_row.size());
    _e_assembler.assembly(e, localEQS);

    // update global
    eqs.addSubMat(local_dofmap_row, local_dofmap_column, localEQS.getMat());
    eqs.addSubRHS(local_dofmap_row, localEQS.getRHSVec());
}

} //DiscreteLib

#endif //ELEMENTWISELINEARSYSTEMUPDATER_TPP_

