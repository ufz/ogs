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
    DofEquationIdTable dofManagerA;
    size_t varId = dofManagerA.addVariableDoFs(0, 0, 10);
    dofManagerA.construct();

    ASSERT_EQ(1u, dofManagerA.getNumberOfVariables());
    ASSERT_EQ(1u, dofManagerA.getNumberOfMeshes());
    //ASSERT_EQ(dofMap1->getNumberOfDiscretePoints(), 10);
    ASSERT_EQ(0u, dofManagerA.mapEqsID(varId, 0, 0));
    ASSERT_EQ(9u, dofManagerA.mapEqsID(varId, 0, 9));
};

TEST(DiscreteLib, DoF_var2_byDof)
{
    DofEquationIdTable dofManagerB;
    size_t dofIdB1 = dofManagerB.addVariableDoFs(0, 0, 10);
    size_t dofIdB2 = dofManagerB.addVariableDoFs(0, 0, 10);
    dofManagerB.setNumberingType(DofNumberingType::BY_VARIABLE);
    dofManagerB.construct();
    ASSERT_EQ(2u, dofManagerB.getNumberOfVariables());
    //ASSERT_EQ(dofManagerB.getTotalNumberOfDiscretePoints(), 20);
    ASSERT_EQ(0u, dofManagerB.mapEqsID(dofIdB1, 0, 0));
    ASSERT_EQ(9u, dofManagerB.mapEqsID(dofIdB1, 0, 9));
    ASSERT_EQ(10u, dofManagerB.mapEqsID(dofIdB2, 0, 0));
    ASSERT_EQ(19u, dofManagerB.mapEqsID(dofIdB2, 0, 9));
};

TEST(DiscreteLib, DoF_var2_byPoint)
{
    DofEquationIdTable dofManagerB;
    size_t dofIdB1 = dofManagerB.addVariableDoFs(0, 0, 10);
    size_t dofIdB2 = dofManagerB.addVariableDoFs(0, 0, 10);
    dofManagerB.setNumberingType(DofNumberingType::BY_POINT);
    dofManagerB.construct();
    ASSERT_EQ(2u, dofManagerB.getNumberOfVariables());
    ASSERT_EQ(1u, dofManagerB.getNumberOfMeshes());
    ASSERT_EQ(0u, dofManagerB.mapEqsID(dofIdB1, 0, 0));
    ASSERT_EQ(18u, dofManagerB.mapEqsID(dofIdB1, 0, 9));
    ASSERT_EQ(1u, dofManagerB.mapEqsID(dofIdB2, 0, 0));
    ASSERT_EQ(19u, dofManagerB.mapEqsID(dofIdB2, 0, 9));
}

TEST(DiscreteLib, DoF_inactive)
{
    DofEquationIdTable *dofManager = new DofEquationIdTable();
    size_t varId = dofManager->addVariableDoFs(0, 0, 9);
    dofManager->deactivateDoFs(varId, 0, 8);
    dofManager->setNumberingType(DofNumberingType::BY_POINT);
    dofManager->construct();

    ASSERT_EQ(1u, dofManager->getNumberOfVariables());
    ASSERT_EQ(8u, dofManager->getTotalNumberOfActiveDoFs());
    ASSERT_EQ(8u, dofManager->getTotalNumberOfActiveDoFsWithoutGhost());
    ASSERT_EQ(0u, dofManager->mapEqsID(varId, 0, 0));
    ASSERT_EQ((size_t)-1, dofManager->mapEqsID(varId, 0, 8));

    delete dofManager;
}

TEST(DiscreteLib, DoF_ghost_nodes)
{
    {
        DofEquationIdTable *dofManager = new DofEquationIdTable();
        int ghost_nodes[] = {5, 6, 7, 8};
        std::vector<size_t> vec_ghost_nodes(ghost_nodes, ghost_nodes+4);
        size_t varId = dofManager->addVariableDoFs(0, 0, 9);
        dofManager->setGhostPoints(0, vec_ghost_nodes);
        dofManager->setNumberingType(DofNumberingType::BY_POINT);
        dofManager->construct();

        ASSERT_EQ(1u, dofManager->getNumberOfVariables());
        ASSERT_EQ(9u, dofManager->getTotalNumberOfActiveDoFs());
        ASSERT_EQ(5u, dofManager->getTotalNumberOfActiveDoFsWithoutGhost());
        ASSERT_EQ(0u, dofManager->mapEqsID(varId, 0, 0));
        ASSERT_TRUE(dofManager->isGhostPoint(0, 8));
        ASSERT_EQ(8u, dofManager->mapEqsID(varId, 0, 8));

        delete dofManager;
    }
    {
        DofEquationIdTable *dofManager = new DofEquationIdTable();
        int ghost_nodes[] = {0, 1, 2, 3, 4};
        std::vector<size_t> vec_ghost_nodes(ghost_nodes, ghost_nodes+5);
        size_t varId = dofManager->addVariableDoFs(0, 0, 9);
        dofManager->setGhostPoints(0, vec_ghost_nodes);
        dofManager->setNumberingType(DofNumberingType::BY_POINT);
        dofManager->construct();

        ASSERT_EQ(1u, dofManager->getNumberOfVariables());
        ASSERT_EQ(9u, dofManager->getTotalNumberOfActiveDoFs());
        ASSERT_EQ(4u, dofManager->getTotalNumberOfActiveDoFsWithoutGhost());
        ASSERT_TRUE(dofManager->isGhostPoint(0, 0));
        ASSERT_EQ(8u, dofManager->mapEqsID(varId, 0, 4));
        ASSERT_EQ(0u, dofManager->mapEqsID(varId, 0, 5));
        delete dofManager;
    }
}
