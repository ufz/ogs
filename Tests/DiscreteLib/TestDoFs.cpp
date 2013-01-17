/**
 * \file
 * \author Norihiro Watanabe
 * \date   2012-08-03
 * \brief  Implementation tests.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <vector>

#include <gtest/gtest.h>

#include "BaseLib/CodingTools.h"
#include "DiscreteLib/DoF/DofEquationIdTable.h"

using namespace DiscreteLib;

//# DoF ###################################################################################################
TEST(DiscreteLib, DoF_var1)
{
    const std::size_t mesh_id = 0;
    const std::size_t pt_id0 = 0;
    const std::size_t n_pt = 10;

    DofEquationIdTable dofManagerA;
    const std::size_t varId = dofManagerA.addVariableDoFs(mesh_id, pt_id0, n_pt);
    dofManagerA.construct();

    //dofManagerA.printout();

    ASSERT_EQ(1u, dofManagerA.getNumberOfVariables());
    ASSERT_EQ(1u, dofManagerA.getNumberOfMeshes());
    ASSERT_EQ(0u, dofManagerA.mapEqsID(DoF(varId, mesh_id, 0)));
    ASSERT_EQ(9u, dofManagerA.mapEqsID(DoF(varId, mesh_id, 9)));
};

TEST(DiscreteLib, DoF_var2_byDof)
{
    const std::size_t mesh_id = 0;
    const std::size_t pt_id0 = 0;
    const std::size_t n_pt = 10;

    DofEquationIdTable dofManagerB;
    const std::size_t varIdB1 = dofManagerB.addVariableDoFs(mesh_id, pt_id0, n_pt);
    const std::size_t varIdB2 = dofManagerB.addVariableDoFs(mesh_id, pt_id0, n_pt);
    dofManagerB.setNumberingType(DofNumberingType::BY_VARIABLE);
    dofManagerB.construct();
    //dofManagerB.printout();

    ASSERT_EQ(2u, dofManagerB.getNumberOfVariables());
    //ASSERT_EQ(dofManagerB.getTotalNumberOfDiscretePoints(), 20);
    ASSERT_EQ(0u, dofManagerB.mapEqsID(DoF(varIdB1, mesh_id, 0)));
    ASSERT_EQ(9u, dofManagerB.mapEqsID(DoF(varIdB1, mesh_id, 9)));
    ASSERT_EQ(10u, dofManagerB.mapEqsID(DoF(varIdB2, mesh_id, 0)));
    ASSERT_EQ(19u, dofManagerB.mapEqsID(DoF(varIdB2, mesh_id, 9)));
};

TEST(DiscreteLib, DoF_var2_byPoint)
{
    DofEquationIdTable dofManagerB;
    std::size_t varIdB1 = dofManagerB.addVariableDoFs(0, 0, 10);
    std::size_t varIdB2 = dofManagerB.addVariableDoFs(0, 0, 10);
    dofManagerB.setNumberingType(DofNumberingType::BY_POINT);
    dofManagerB.construct();
    //dofManagerB.printout();
    ASSERT_EQ(2u, dofManagerB.getNumberOfVariables());
    ASSERT_EQ(1u, dofManagerB.getNumberOfMeshes());
    ASSERT_EQ(0u, dofManagerB.mapEqsID(DoF(varIdB1, 0, 0)));
    ASSERT_EQ(18u, dofManagerB.mapEqsID(DoF(varIdB1, 0, 9)));
    ASSERT_EQ(1u, dofManagerB.mapEqsID(DoF(varIdB2, 0, 0)));
    ASSERT_EQ(19u, dofManagerB.mapEqsID(DoF(varIdB2, 0, 9)));
}

TEST(DiscreteLib, DoF_inactive)
{
    DofEquationIdTable *dofManager = new DofEquationIdTable();
    std::size_t varId = dofManager->addVariableDoFs(0, 0, 9);
    dofManager->activateDoF(DoF(varId, 0, 8), false);
    //dofManager->activateDoF(varId, 0, 8, false);
    dofManager->setNumberingType(DofNumberingType::BY_POINT);
    dofManager->construct();

    ASSERT_EQ(1u, dofManager->getNumberOfVariables());
    ASSERT_EQ(8u, dofManager->getTotalNumberOfActiveDoFs());
    ASSERT_EQ(8u, dofManager->getTotalNumberOfActiveDoFsWithoutGhost());
    ASSERT_EQ(0u, dofManager->mapEqsID(DoF(varId, 0, 0)));
    ASSERT_EQ((std::size_t)-1, dofManager->mapEqsID(DoF(varId, 0, 8)));

    delete dofManager;
}

TEST(DiscreteLib, DoF_ghost_nodes)
{
    const std::size_t mesh_id = 0;
    {
        DofEquationIdTable *dofManager = new DofEquationIdTable();
        std::size_t ghost_nodes[] = {5, 6, 7, 8};
        const std::vector<std::size_t> vec_ghost_nodes(ghost_nodes, ghost_nodes+4);
        const std::size_t varId = dofManager->addVariableDoFs(mesh_id, 0, 9);
        dofManager->setGhostPoints(mesh_id, vec_ghost_nodes);
        dofManager->setNumberingType(DofNumberingType::BY_POINT);
        dofManager->construct();

        ASSERT_EQ(1u, dofManager->getNumberOfVariables());
        ASSERT_EQ(9u, dofManager->getTotalNumberOfActiveDoFs());
        ASSERT_EQ(5u, dofManager->getTotalNumberOfActiveDoFsWithoutGhost());
        ASSERT_EQ(0u, dofManager->mapEqsID(DoF(varId, mesh_id, 0)));
        ASSERT_TRUE(dofManager->isGhostPoint(mesh_id, 8));
        ASSERT_EQ(8u, dofManager->mapEqsID(DoF(varId, mesh_id, 8)));

        delete dofManager;
    }
    {
        DofEquationIdTable *dofManager = new DofEquationIdTable();
        std::size_t ghost_nodes[] = {0, 1, 2, 3, 4};
        const std::vector<std::size_t> vec_ghost_nodes(ghost_nodes, ghost_nodes+5);
        const std::size_t varId = dofManager->addVariableDoFs(mesh_id, 0, 9);
        dofManager->setGhostPoints(mesh_id, vec_ghost_nodes);
        dofManager->setNumberingType(DofNumberingType::BY_POINT);
        dofManager->construct();
        //dofManager->printout();

        ASSERT_EQ(1u, dofManager->getNumberOfVariables());
        ASSERT_EQ(9u, dofManager->getTotalNumberOfActiveDoFs());
        ASSERT_EQ(4u, dofManager->getTotalNumberOfActiveDoFsWithoutGhost());
        ASSERT_TRUE(dofManager->isGhostPoint(mesh_id, 0));
        ASSERT_EQ(8u, dofManager->mapEqsID(DoF(varId, mesh_id, 4)));
        ASSERT_EQ(0u, dofManager->mapEqsID(DoF(varId, mesh_id, 5)));
        delete dofManager;
    }
}
