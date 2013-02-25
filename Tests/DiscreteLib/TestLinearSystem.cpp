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
#include <boost/property_tree/ptree.hpp>

#include <vector>
#include <memory>

#include "BaseLib/CodingTools.h"
#ifdef USE_LIS
    #include "MathLib/LinAlg/SystemOfLinearEquations/LisLinearSystem.h"
#endif  // USE_LIS
#include "MathLib/LinAlg/SystemOfLinearEquations/DenseLinearSystem.h"
#include "MeshLib/MeshGenerator.h"
#include "DiscreteLib/ElementWiseManipulator/IElemenetWiseLinearSystemLocalAssembler.h"
#include "DiscreteLib/ElementWiseManipulator/IElemenetWiseVectorLocalAssembler.h"
#include "DiscreteLib/ElementWiseManipulator/ElementWiseLinearSystemUpdater.h"
#include "DiscreteLib/ElementWiseManipulator/ElementWiseVectorUpdater.h"
#include "DiscreteLib/Serial/SerialDiscreteSystem.h"
#include "DiscreteLib/SparsityBuilder/SparsityBuilderDummy.h"

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


template <class MyDiscreteSystem, class MyLinearSolver, class MySparsityBuilderType, class MyLocalAssembler>
void solveLinear(   boost::property_tree::ptree &t_root,
                    const Mesh* msh,
                    const MyLocalAssembler &ele_assembler,
                    const std::vector<std::size_t> &vec_DirichletBC_id,
                    const std::vector<double> &vec_DirichletBC_value,
                    std::vector<double> &x)
{
    // define a discrete system
    MyDiscreteSystem dis(msh);

    //----------------------------------------------------------------------
    // the following codes should be independent of a parallelization scheme
    //----------------------------------------------------------------------
    typedef typename MyDiscreteSystem::template MyLinearSystem<MyLinearSolver,MySparsityBuilderType>::type MyDiscreteLinearSystem;
    typedef DiscreteLib::ElementWiseLinearSystemUpdater<MyLocalAssembler,MyLinearSolver> MyElementWiseLinearSystemUpdater;
    typedef typename MyDiscreteSystem::template MyLinearSystemAssembler<MyElementWiseLinearSystemUpdater,MyLinearSolver>::type MyGlobalAssembler;

    // set up a DoF table
    DofEquationIdTable dofMappingTable;
    const std::size_t varId = dofMappingTable.addVariableDoFs(msh->getID(), 0, msh->getNNodes());
    dofMappingTable.construct(DofNumberingType::BY_VARIABLE);

    // create a linear system
    MyDiscreteLinearSystem *linear_eq = dis.template createLinearSystem<MyLinearSolver, MySparsityBuilderType>(&dofMappingTable);
    linear_eq->getLinearSolver()->setOption(t_root);

    // set Dirichlet BC (this should be done before construct)
    linear_eq->setKnownSolution(varId, vec_DirichletBC_id, vec_DirichletBC_value);

    // construct the system
    MyElementWiseLinearSystemUpdater updater(msh->getID(), ele_assembler);
    MyGlobalAssembler assembler(updater);
    linear_eq->construct(assembler);

    // solve the linear
    //linear_eq->getLinearSolver()->printout();
    linear_eq->solve();
    //linear_eq->getLinearSolver()->printout();

    // check a solution
    linear_eq->getSolVec(x);
}

} //namespace

TEST(DiscreteLib, LinearSerialDense)
{
    // example
    SteadyDiffusion2DExample1 ex1;
    typedef SteadyDiffusion2DExample1::LocalAssembler MyLocalAssembler;
    MyLocalAssembler ele_assembler;

    //----------------------------------------------------------------------
    // select DiscreteSystem and LinearSolver
    //----------------------------------------------------------------------
    typedef SerialDiscreteSystem MyDiscreteSystem;
    typedef DenseLinearSystem MyLinearSolver;
    typedef SparsityBuilderDummy MySparsityBuilderType;
    boost::property_tree::ptree t_root;

    std::vector<double> x;
    solveLinear<MyDiscreteSystem, MyLinearSolver, MySparsityBuilderType, MyLocalAssembler>
        (t_root, ex1.msh, ele_assembler, ex1.vec_DirichletBC_id, ex1.vec_DirichletBC_value, x);

    ASSERT_DOUBLE_ARRAY_EQ(&ex1.exact_solutions[0], &x[0], 9, 1.e-5);
}

#ifdef USE_LIS
TEST(DiscreteLib, LinearSerialLis1)
{
    // example
    SteadyDiffusion2DExample1 ex1;
    typedef SteadyDiffusion2DExample1::LocalAssembler MyLocalAssembler;
    MyLocalAssembler ele_assembler;

    //----------------------------------------------------------------------
    // select DiscreteSystem and LinearSolver
    //----------------------------------------------------------------------
    typedef SerialDiscreteSystem MyDiscreteSystem;
    typedef LisLinearSystem MyLinearSolver;
    typedef SparsityBuilderDummy MySparsityBuilderType;

    // set solver options using Boost property tree
    boost::property_tree::ptree t_root;
    boost::property_tree::ptree t_solver;
    t_solver.put("solver_type", "CG");
    t_solver.put("precon_type", "NONE");
    t_solver.put("matrix_type", "CCS");
    t_solver.put("error_tolerance", 1e-15);
    t_solver.put("max_iteration_step", 1000);
    t_root.put_child("LinearSolver", t_solver);

    std::vector<double> x;
    solveLinear<MyDiscreteSystem, MyLinearSolver, MySparsityBuilderType, MyLocalAssembler>
        (t_root, ex1.msh, ele_assembler, ex1.vec_DirichletBC_id, ex1.vec_DirichletBC_value, x);

    ASSERT_DOUBLE_ARRAY_EQ(&ex1.exact_solutions[0], &x[0], 9, 1.e-5);
}
#endif

