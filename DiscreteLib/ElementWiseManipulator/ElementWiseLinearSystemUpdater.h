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

#include "MeshLib/Elements/Element.h"
#include "DiscreteLib/DoF/DofEquationIdTable.h"

namespace DiscreteLib
{

/**
 * \brief Element-wise linear system updater
 *
 * \tparam T_LOCAL_ASSEMBLER    Local assembler
 * \tparam T_SOLVER             Linear solver
 */
template <class T_LOCAL_ASSEMBLER, class T_SOLVER>
class ElementWiseLinearSystemUpdater
{
public:
    typedef T_LOCAL_ASSEMBLER LocalAssemblerType;
    typedef T_SOLVER SolverType;

    /**
     * Constructor
     * @param msh_id            Mesh ID
     * @param local_assembler   Local assembler object
     */
    ElementWiseLinearSystemUpdater(std::size_t msh_id, const LocalAssemblerType &local_assembler)
    : _msh_id(msh_id), _e_assembler(local_assembler)
    {};

    /**
     * Update a linear system for the given element
     * @param e     mesh element
     * @param eqs   linear system
     */
    void update(const MeshLib::Element &e, const DofEquationIdTable &dofManager, SolverType &eqs);

private:
    const std::size_t _msh_id;
    LocalAssemblerType _e_assembler;
};

} //end

#include "ElementWiseLinearSystemUpdater.tpp"

#endif //ELEMENTWISELINEAREQUATIONUPDATER_H_

