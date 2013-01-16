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
#include <memory>

#include "BaseLib/CodingTools.h"
#include "MathLib/LinAlg/SystemOfLinearEquations/LisLinearSystem.h"
#include "MeshLib/MeshGenerator.h"
#include "DiscreteLib/ElementWiseManipulator/IElemenetWiseLinearSystemLocalAssembler.h"
#include "DiscreteLib/ElementWiseManipulator/ElementWiseLinearSystemUpdater.h"
#include "DiscreteLib/Serial/SerialDiscreteSystem.h"

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

struct SteadyDiffusion2DExample1
{
    class LocalAssembler : public DiscreteLib::IElemenetWiseLinearSystemLocalAssembler
    {
        LocalMatrix _m;
    public:
        LocalAssembler()
        : _m(4,4)
        {
            _m(0,0) = 4.0; _m(0,1) = -1.0; _m(0,2) = -2.0; _m(0,3) = -1.0; 
            _m(1,1) = 4.0; _m(1,2) = -1.0; _m(1,3) = -2.0;
            _m(2,2) = 4.0; _m(2,3) = -1.0;
            _m(3,3) = 4.0;
            // copy upper triangle to lower to make symmetric
            for (std::size_t i=0; i<4; i++)
                for (std::size_t j=0; j<i; j++)
                    _m(i,j) = _m(j,i);
            //_m *= 1.e-11/6.0;
            for (std::size_t i=0; i<4; i++)
                for (std::size_t j=0; j<4; j++)
                    _m(i,j) *= 1.e-11/6.0;
        }

        void assembly(const MeshLib::Element &/*e*/, LocalLinearSystem &eqs)
        {
            //eqs.getMat() = _m;
            for (std::size_t i=0; i<4; i++)
                for (std::size_t j=0; j<4; j++)
                    eqs.getMat()(i,j) = _m(i,j);
        }
    };

    static const std::size_t dim_eqs = 9;
    MeshLib::Mesh* msh;
    std::vector<std::size_t> vec_DirichletBC_id;
    std::vector<double> vec_DirichletBC_value;
    std::vector<double> exact_solutions;

    SteadyDiffusion2DExample1()
    {
        msh = MeshGenerator::generateRegularQuadMesh(2.0, 2, .0, .0, .0);
        std::size_t int_dirichlet_bc_id[] = {2,5,8,0,3,6};
        vec_DirichletBC_id.assign(int_dirichlet_bc_id, int_dirichlet_bc_id+6);
        vec_DirichletBC_value.resize(6);
        fill(vec_DirichletBC_value.begin(), vec_DirichletBC_value.begin()+3, .0);
        fill(vec_DirichletBC_value.begin()+3, vec_DirichletBC_value.end(), 1.0);
        exact_solutions.resize(9);
        for (std::size_t i=0; i<9; i++) {
            if (i%3==0) exact_solutions[i] = 1.0;
            if (i%3==1) exact_solutions[i] = 0.5;
            if (i%3==2) exact_solutions[i] = 0.;
        }
    }

    ~SteadyDiffusion2DExample1()
    {
        delete msh;
    }
}; //SteadyDiffusion2DExample1

} //namespace


TEST(DiscreteLib, SerialVec1)
{
    std::unique_ptr<MeshLib::Mesh> msh(MeshGenerator::generateRegularQuadMesh(2.0, 2, .0, .0, .0));

    SerialDiscreteSystem dis(msh.get());
    SerialDiscreteVector<double> *v = dis.createVector<double>(msh->getNNodes());

    double expected[] = {0, 1, 2, 3, 4, 5, 6, 7, 8};
    for (std::size_t i=0; i<9; i++)
        (*v)[i] = expected[i];
    
    ASSERT_EQ(9u, v->size());
    ASSERT_DOUBLE_ARRAY_EQ(expected, &(*v)[0], 9);
}

#ifdef USE_LIS
TEST(DiscreteLib, Lis1)
{
    // example
    SteadyDiffusion2DExample1 ex1;
    typedef SteadyDiffusion2DExample1::LocalAssembler MyLocalAssembler;
    MyLocalAssembler ele_assembler;
    const MeshLib::Mesh* msh = ex1.msh;

    //----------------------------------------------------------------------
    // select DiscreteSystem and LinearSolver
    //----------------------------------------------------------------------
    typedef SerialDiscreteSystem MyDiscreteSystem;
    typedef LisLinearSystem MyLinearSolver;

    {
        // define a discrete system
        MyDiscreteSystem dis(msh);

        //----------------------------------------------------------------------
        // the following codes should be independent of a parallelization scheme
        //----------------------------------------------------------------------
        typedef MyDiscreteSystem::MyLinearSystem<MyLinearSolver,SparsityBuilderDummy>::type MyDiscreteLinearSystem;
        typedef DiscreteLib::ElementWiseLinearSystemUpdater<MyLocalAssembler,MyLinearSolver> MyElementWiseLinearSystemUpdater;
        typedef MyDiscreteSystem::MyLinearSystemAssembler<MyElementWiseLinearSystemUpdater,MyLinearSolver>::type MyGlobalAssembler;

        // set up a DoF table
        DofEquationIdTable dofMappingTable;
        const std::size_t varId = dofMappingTable.addVariableDoFs(msh->getID(), 0, msh->getNNodes());
        dofMappingTable.construct(DofNumberingType::BY_VARIABLE);

        // create a linear system
        MyDiscreteLinearSystem *linear_eq = dis.createLinearSystem<MyLinearSolver, SparsityBuilderDummy>(&dofMappingTable);
        MyLinearSolver* lis = linear_eq->getLinearSolver();
        lis->getOption().solver_type = LisOption::SolverType::CG;
        lis->getOption().precon_type = LisOption::PreconType::NONE;

        // construct the system
        MyElementWiseLinearSystemUpdater updater(msh->getID(), &ele_assembler);
        MyGlobalAssembler assembler(&updater);
        linear_eq->construct(assembler);

        // apply Dirichlet BC
        linear_eq->setKnownSolution(varId, ex1.vec_DirichletBC_id, ex1.vec_DirichletBC_value);
        //linear_eq->getLinearEquation()->printout();

        // solve the linear
        linear_eq->solve();

        // check a solution
        std::vector<double> x;
        linear_eq->getSolVec(x);
        ASSERT_DOUBLE_ARRAY_EQ(&ex1.exact_solutions[0], &x[0], 9, 1.e-5);
    }
}
#endif

