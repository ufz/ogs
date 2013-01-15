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

#include <gtest/gtest.h>

#include <vector>

#include "BaseLib/CodingTools.h"
#include "MathLib/LinAlg/SystemOfLinearEquations/LisLinearSystem.h"
#include "MeshLib/MeshGenerator.h"
#include "DiscreteLib/ElementWiseManipulator/IElemenetWiseLinearEquationLocalAssembler.h"
#include "DiscreteLib/ElementWiseManipulator/ElementWiseLinearSystemUpdater.h"
#include "DiscreteLib/Serial/DiscreteSystem.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace MathLib;
using namespace MeshLib;
using namespace DiscreteLib;

namespace
{

inline void ASSERT_DOUBLE_ARRAY_EQ(const double* Expected, const double* Actual, size_t N, double epsilon=1.0e-8) {
    for (size_t i=0; i<N; i++) \
        ASSERT_NEAR(Expected[i], Actual[i], epsilon);
}
}

struct DiscreteExample1
{
    std::vector<std::size_t> list_dirichlet_bc_id;
    std::vector<double> list_dirichlet_bc_value;
    static const std::size_t dim_eqs = 9;
    std::vector<double> exH;

    DiscreteExample1()
    {
        std::size_t int_dirichlet_bc_id[] = {2,5,8,0,3,6};
        list_dirichlet_bc_id.assign(int_dirichlet_bc_id, int_dirichlet_bc_id+6);
        list_dirichlet_bc_value.resize(6);
        fill(list_dirichlet_bc_value.begin(), list_dirichlet_bc_value.begin()+3, .0);
        fill(list_dirichlet_bc_value.begin()+3, list_dirichlet_bc_value.end(), 1.0);
        exH.resize(9);
        for (std::size_t i=0; i<9; i++) {
            if (i%3==0) exH[i] = 1.0;
            if (i%3==1) exH[i] = 0.5;
            if (i%3==2) exH[i] = 0.;
        }
    }

    class TestElementAssembler : public DiscreteLib::IElemenetWiseLinearSystemLocalAssembler
    {
        LocalMatrix _m;
    public:
        TestElementAssembler()
        : _m(4,4)
        {
            //_m.resize(4,4);
            _m(0,0) = 4.0; _m(0,1) = -1.0; _m(0,2) = -2.0; _m(0,3) = -1.0; 
            _m(1,1) = 4.0; _m(1,2) = -1.0; _m(1,3) = -2.0;
            _m(2,2) = 4.0; _m(2,3) = -1.0;
            _m(3,3) = 4.0;
            for (std::size_t i=0; i<4; i++)
                for (std::size_t j=0; j<i; j++) _m(i,j) = _m(j,i);
            //_m *= 1.e-11/6.0;
        }
        void assembly(const MeshLib::Element &/*e*/, LocalLinearSystem &eqs)
        {
            //(*eqs.getA()) = _m;
        }
    };
};

TEST(Discrete, VecSingle1)
{
    MeshLib::Mesh *msh = MeshGenerator::generateRegularQuadMesh(2.0, 2, .0, .0, .0);
    DiscreteSystem dis(msh);
    DiscreteVector<double> *v = dis.createVector<double>(msh->getNNodes());

    double expected[] = {0, 1, 2, 3, 4, 5, 6, 7, 8};
    std::copy(expected, expected+9, v->begin());
    
    ASSERT_EQ(9u, v->size());
    ASSERT_DOUBLE_ARRAY_EQ(expected, &(*v)[0], 9);
}

#ifdef USE_LIS
TEST(Discrete, Lis1)
{
    DiscreteExample1 ex1;
    DiscreteExample1::TestElementAssembler ele_assembler;
    LisLinearSystem lis(ex1.dim_eqs);
    lis.getOption().solver_type = LisOption::SolverType::CG;
    lis.getOption().precon_type = LisOption::PreconType::NONE;
    MeshLib::Mesh *msh = MeshGenerator::generateRegularQuadMesh(2.0, 2, .0, .0, .0);

    // define discrete system
    DiscreteSystem dis(msh);
    {
        // DoF?
        DofEquationIdTable dofManager;
        dofManager.addVariableDoFs(0, 0, msh->getNNodes());
        dofManager.construct(DofNumberingType::BY_VARIABLE);
        // create a linear problem
        IDiscreteLinearSystem *linear_eq = dis.createLinearEquation<LisLinearSystem, SparsityBuilderDummy>(&lis, &dofManager);
        // solve the equation
        linear_eq->initialize();
        linear_eq->setPrescribedDoF(0, ex1.list_dirichlet_bc_id, ex1.list_dirichlet_bc_value);
        typedef DiscreteLib::ElementWiseLinearSystemUpdater<DiscreteExample1::TestElementAssembler,LisLinearSystem> MyUpdater;
        typedef DiscreteSystem::MyLinearEquationAssembler<MyUpdater,LisLinearSystem>::type MyGlobalAssembler;
        MyUpdater updater(msh, &ele_assembler);
        MyGlobalAssembler assembler(&updater);
        linear_eq->construct(assembler);
        //linear_eq->getLinearEquation()->printout();
        linear_eq->solve();

        ASSERT_DOUBLE_ARRAY_EQ(&ex1.exH[0], linear_eq->getLocalX(), 9, 1.e-5);
    }
}
#endif

