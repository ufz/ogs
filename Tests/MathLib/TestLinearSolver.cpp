/**
 * \file
 * \author Norihiro Watanabe
 * \date   2013-04-16
 * \brief  Implementation tests.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include<type_traits> // for is_constructible


#include <gtest/gtest.h>
#include <boost/property_tree/ptree.hpp>

#include "MathLib/LinAlg/Dense/DenseVector.h"
#include "MathLib/LinAlg/Dense/DenseMatrix.h"
#include "MathLib/LinAlg/Dense/GlobalDenseMatrix.h"
#include "MathLib/LinAlg/Dense/DenseTools.h"
#include "MathLib/LinAlg/FinalizeMatrixAssembly.h"
#include "MathLib/LinAlg/Solvers/GaussAlgorithm.h"

#ifdef USE_LIS
#include "MathLib/LinAlg/Lis/LisLinearSolver.h"
#include "MathLib/LinAlg/Lis/LisTools.h"
#endif

#if defined(USE_PETSC)
// Test class MathLib::PETScLinearEquation with  the overload interface
#include "MathLib/LinAlg/PETSc/PETScLinearEquation.h"

// Test class MathLib::PETScLinearEquation, PETScMatrix and PETScVector with  the overload interface
#include "MathLib/LinAlg/PETSc/PETScMatrix.h"
#include "MathLib/LinAlg/PETSc/PETScVector.h"
#include "MathLib/LinAlg/PETSc/PETScLinearSolver.h"
#endif

#include "../TestTools.h"

namespace
{

template<class T_Mat>
void setMatrix9x9(T_Mat &mat)
{
    double d_mat[] = {
        6.66667e-012, -1.66667e-012, 0, -1.66667e-012, -3.33333e-012, 0, 0, 0, 0,
        -1.66667e-012, 1.33333e-011, -1.66667e-012, -3.33333e-012, -3.33333e-012, -3.33333e-012, 0, 0, 0,
        0, -1.66667e-012, 6.66667e-012, 0, -3.33333e-012, -1.66667e-012, 0, 0, 0,
        -1.66667e-012, -3.33333e-012, 0, 1.33333e-011, -3.33333e-012, 0, -1.66667e-012, -3.33333e-012, 0,
        -3.33333e-012, -3.33333e-012, -3.33333e-012, -3.33333e-012, 2.66667e-011, -3.33333e-012, -3.33333e-012, -3.33333e-012, -3.33333e-012,
        0, -3.33333e-012, -1.66667e-012, 0, -3.33333e-012, 1.33333e-011, 0, -3.33333e-012, -1.66667e-012,
        0, 0, 0, -1.66667e-012, -3.33333e-012, 0, 6.66667e-012, -1.66667e-012, 0,
        0, 0, 0, -3.33333e-012, -3.33333e-012, -3.33333e-012, -1.66667e-012, 1.33333e-011, -1.66667e-012,
        0, 0, 0, 0, -3.33333e-012, -1.66667e-012, 0, -1.66667e-012, 6.66667e-012
    };
	for (unsigned i = 0; i < 9; i++)
		for (unsigned j = 0; j < 9; j++)
			mat.setValue(i, j, d_mat[i*9+j]);
}

struct Example1
{
    MathLib::GlobalDenseMatrix<double> mat;
    std::vector<size_t> vec_dirichlet_bc_id;
    std::vector<double> vec_dirichlet_bc_value;
    static const std::size_t dim_eqs = 9;
    double* exH;

    Example1()
    : mat(dim_eqs, dim_eqs), exH(new double[dim_eqs])
    {
        setMatrix9x9(mat);
        std::size_t int_dirichlet_bc_id[] = {2,5,8,0,3,6};
        vec_dirichlet_bc_id.assign(int_dirichlet_bc_id, int_dirichlet_bc_id+6);
        vec_dirichlet_bc_value.resize(6);
        std::fill(vec_dirichlet_bc_value.begin(), vec_dirichlet_bc_value.begin()+3, .0);
        std::fill(vec_dirichlet_bc_value.begin()+3, vec_dirichlet_bc_value.end(), 1.0);
        for (std::size_t i=0; i<9; i++) {
            if (i%3==0) exH[i] = 1.0;
            if (i%3==1) exH[i] = 0.5;
            if (i%3==2) exH[i] = 0.;
        }
    }

    ~Example1()
    {
        delete [] exH;
    }
};


 
  template <class T_MATRIX, class T_VECTOR,   class T_LINEAR_SOVLER 
	    //            class = typename std::enable_if<not std::is_void<T_VECTOR>::value>::type
           >
void checkLinearSolverInterface(T_MATRIX  &A,  boost::property_tree::ptree &ls_option)
{
    Example1 ex1;

    //    if (!std::is_constructible<T_MATRIX_OR_EXTSOLVER, MathLib::PETScLinearEquation>::value)
    // if (!std::is_constructible<T_VECTOR, void>::value)

    {
       // set a coefficient matrix
       A.setZero();
       for (size_t i=0; i<ex1.dim_eqs; i++) {
           for (size_t j=0; j<ex1.dim_eqs; j++) {
               double v = ex1.mat(i, j);
               if (v!=.0)
                   A.add(i, j, v);
           }
       }


      // set RHS and solution vectors
      T_VECTOR rhs(ex1.dim_eqs);
      T_VECTOR x(ex1.dim_eqs);

      // apply BC
      MathLib::applyKnownSolution(A, rhs, ex1.vec_dirichlet_bc_id, ex1.vec_dirichlet_bc_value);

      //Call inside solver.
      MathLib::finalizeMatrixAssembly(A);

      // solve
      T_LINEAR_SOVLER ls(A, &ls_option);
      ls.solve(rhs, x);
      ASSERT_ARRAY_NEAR(ex1.exH, x, ex1.dim_eqs, 1e-5);
    }
}

  //template argument T_VECTOR will be removed if it is not used anymore
 template <class T_LINEAR_EQUATION, typename T_VECTOR >
void checkLinearSolverInterface(T_LINEAR_EQUATION &l_eqs,  boost::property_tree::ptree &ls_option)
{

    // ---------------------------------------------
    // Test case
    Example1 ex1;



    const int msize = l_eqs.getMPI_Size(); 
    const int mrank = l_eqs.getMPI_Rank(); 

    //-------------------------------------------------------------------
    // Assembly test
    double *local_matrix = nullptr;  // nullptr not support by IBM C++11
    int *idx_r = nullptr; // 
    int *idx_c = nullptr;
    local_matrix = new double[msize * ex1.dim_eqs];
    idx_c = new int[ex1.dim_eqs];
    idx_r = new int[msize];

    
    for(int j=0;j<ex1.dim_eqs; j++)
    {
       idx_c[j] = j; 
    }

    for(int i=0; i<msize; i++)
    {
       idx_r[i] = mrank *  msize + i;
       for(int j=0;j<ex1.dim_eqs; j++)
       {
           local_matrix[i*ex1.dim_eqs +j] =  ex1.mat(idx_r[i], j);
       }
    } 

    //-------------------------------------------------------------------
    //
    // Solver configuration   
    l_eqs.Config(ls_option);
    l_eqs.initializeMatVec();


    // local assembly
    l_eqs.addMatrixEntries(msize, idx_r, ex1.dim_eqs, idx_c, local_matrix);
    // No need to change RHS for this example.    
    l_eqs.finalAssemble();

  
    //-------------------------------------------------------------------
    // Apply Dirichlet BC test
    int *bc_eqs_id = nullptr;
    double *bc_eqs_value = nullptr;
    // the following caculation will be removed when a real function about D-BC is ready
    int bc_size_rank = 0;  
    for(size_t i=0; i<ex1.vec_dirichlet_bc_id.size(); i++)
    {
      const int bc_id =  ex1.vec_dirichlet_bc_id[i]; 
      if(bc_id > msize*mrank && bc_id < msize*(mrank+1))
        bc_size_rank++;  
    }


    if(bc_size_rank > 0)
    {

       bc_eqs_id = new int[bc_size_rank];
       bc_eqs_value = new double[bc_size_rank];
       bc_size_rank = 0;
       for(size_t i=0; i<ex1.vec_dirichlet_bc_id.size(); i++)
       {
          const int bc_id =  ex1.vec_dirichlet_bc_id[i]; 



       if(bc_id >= msize*mrank && bc_id < msize*(mrank+1))
       {
	  bc_eqs_id[bc_size_rank] = bc_id;  
          bc_eqs_value[bc_size_rank] =  ex1.vec_dirichlet_bc_value[i];

          bc_size_rank++;

       }  
    }
    }
    //-------------------------------------------------------------------


    // Apply Dirichlet BC
    l_eqs.applyKnownSolutions(bc_size_rank, bc_eqs_id,  bc_eqs_value);


    // Solve the linear equation
    l_eqs.Solver();


     l_eqs.mappingSolution(); 
     double *x = l_eqs.getGlobalSolution();  //T_VECTOR x, also works, template argument T_VECTOR will be removed 


    // Convergence test
     ASSERT_ARRAY_NEAR(ex1.exH, x, ex1.dim_eqs, 1e-5);


    // Test
    if(bc_eqs_id != nullptr)
      {
        delete [] bc_eqs_id;
      } 
    if(bc_eqs_value != nullptr)
      {
        delete [] bc_eqs_value;
      } 
  
    delete [] local_matrix;
    delete [] idx_c;
    delete [] idx_r;

}


} // end namespace


TEST(MathLib, CheckInterface_GaussAlgorithm)
{
    boost::property_tree::ptree t_root;
    boost::property_tree::ptree t_solver;
    //t_solver.put("solver_package", "Dense");
    t_root.put_child("LinearSolver", t_solver);

    typedef MathLib::GaussAlgorithm<MathLib::GlobalDenseMatrix<double>, MathLib::DenseVector<double> > LinearSolverType;
    MathLib::GlobalDenseMatrix<double> A(Example1::dim_eqs, Example1::dim_eqs);

    checkLinearSolverInterface<MathLib::GlobalDenseMatrix<double>, MathLib::DenseVector<double>, LinearSolverType>(A, t_root);
}


#ifdef USE_LIS
TEST(Math, CheckInterface_Lis)
{
    // set solver options using Boost property tree
    boost::property_tree::ptree t_root;
    boost::property_tree::ptree t_solver;
    //t_solver.put("solver_package", "LIS");
    t_solver.put("solver_type", "CG");
    t_solver.put("precon_type", "NONE");
    t_solver.put("error_tolerance", 1e-15);
    t_solver.put("max_iteration_step", 1000);
    t_root.put_child("LinearSolver", t_solver);

    MathLib::LisMatrix A(Example1::dim_eqs);
    checkLinearSolverInterface<MathLib::LisMatrix, MathLib::LisVector, MathLib::LisLinearSolver>(A, t_root);
}
#endif

#if defined(USE_PETSC)

#define test_p1
#ifdef test_p1
// Test class MathLib::PETScLinearEquation with  the overload interface
#include "MathLib/LinAlg/PETSc/PETScLinearEquation.h"
TEST(Math, CheckInterface_PETSc_1)
{
   int mrank, msize;

   MPI_Comm_rank(PETSC_COMM_WORLD, &mrank);
   MPI_Comm_size(PETSC_COMM_WORLD, &msize);

   if(msize != 3)
   {
      PetscSynchronizedPrintf(PETSC_COMM_WORLD, "===\nThis is test of PETSc solver. The numnber of cores must be 3 exactly");

     PetscFinalize();
     exit(EXIT_FAILURE);
   }

   PetscSynchronizedPrintf(PETSC_COMM_WORLD, "===\nUse PETSc solver");
   PetscSynchronizedPrintf(PETSC_COMM_WORLD, "Number of CPUs: %d, rank: %d\n", msize, mrank);
   PetscSynchronizedFlush(PETSC_COMM_WORLD);


    // set solver options using Boost property tree
    boost::property_tree::ptree t_root;
    boost::property_tree::ptree t_solver;
    //t_solver.put("solver_package", "LIS");
    t_solver.put("solver_type", "cg");
    t_solver.put("precon_type", "none");
    t_solver.put("error_tolerance", 1e-15);
    t_solver.put("max_iteration_step", 1000);
    t_root.put_child("LinearSolver", t_solver);



    int sparse_info[4] = {Example1::dim_eqs, Example1::dim_eqs, Example1::dim_eqs, Example1::dim_eqs};

    MathLib::PETScLinearEquation petsc_leq;
    petsc_leq.set_rank_size(mrank, msize);
    petsc_leq.Init(Example1::dim_eqs, sparse_info);
  
    checkLinearSolverInterface<MathLib::PETScLinearEquation, std::vector<double>>(petsc_leq, t_root);

}

#endif

#define test_p2
#ifdef test_p2
// Test class MathLib::PETScLinearEquation, PETScMatrix and PETScVector with  the overload interface
TEST(Math, CheckInterface_PETSc_2)
{
   int mrank, msize;

   MPI_Comm_rank(PETSC_COMM_WORLD, &mrank);
   MPI_Comm_size(PETSC_COMM_WORLD, &msize);

   if(msize != 3)
   {
      PetscSynchronizedPrintf(PETSC_COMM_WORLD, "===\nThis is test of PETSc solver. The numnber of cores must be 3 exactly");

     PetscFinalize();
     exit(EXIT_FAILURE);
   }

   PetscSynchronizedPrintf(PETSC_COMM_WORLD, "===\nUse PETSc solver");
   PetscSynchronizedPrintf(PETSC_COMM_WORLD, "Number of CPUs: %d, rank: %d\n", msize, mrank);
   PetscSynchronizedFlush(PETSC_COMM_WORLD);


    // set solver options using Boost property tree
    boost::property_tree::ptree t_root;
    boost::property_tree::ptree t_solver;
    //t_solver.put("solver_package", "LIS");
    t_solver.put("solver_type", "cg");
    t_solver.put("precon_type", "none");
    t_solver.put("error_tolerance", 1e-15);
    t_solver.put("max_iteration_step", 1000);
    t_root.put_child("LinearSolver", t_solver);



    int sparse_info[4] = {Example1::dim_eqs, Example1::dim_eqs, Example1::dim_eqs, Example1::dim_eqs};


    MathLib::PETScMatrix A;
    A.set_rank_size(mrank, msize);
    A.Init(Example1::dim_eqs, sparse_info);

    MathLib::PETScVector b(Example1::dim_eqs);
    MathLib::PETScVector x;
    b.CloneMe(x);

 
    MathLib::PETScLinearSolver petsc_leq(A, b, x);
    //petsc_leq.set_rank_size(mrank, msize);
    //x.set_rank_size(mrank, msize);
    //b.set_rank_size(mrank, msize);
    petsc_leq.Init(Example1::dim_eqs);

    checkLinearSolverInterface<MathLib::PETScLinearSolver, std::vector<double>>(petsc_leq, t_root);

}
#endif

#define test_p3
#ifdef test_p3
// Test class MathLib::PETScLinearEquation, PETScMatrix and PETScVector with  the overload interface
TEST(Math, CheckInterface_PETSc_3)
{
   int mrank, msize;

   MPI_Comm_rank(PETSC_COMM_WORLD, &mrank);
   MPI_Comm_size(PETSC_COMM_WORLD, &msize);

   if(msize != 3)
   {
      PetscSynchronizedPrintf(PETSC_COMM_WORLD, "===\nThis is test of PETSc solver. The numnber of cores must be 3 exactly");

     PetscFinalize();
     exit(EXIT_FAILURE);
   }

   PetscSynchronizedPrintf(PETSC_COMM_WORLD, "===\nUse PETSc solver");
   PetscSynchronizedPrintf(PETSC_COMM_WORLD, "Number of CPUs: %d, rank: %d\n", msize, mrank);
   PetscSynchronizedFlush(PETSC_COMM_WORLD);


    // set solver options using Boost property tree
    boost::property_tree::ptree t_root;
    boost::property_tree::ptree t_solver;
    //t_solver.put("solver_package", "LIS");
    t_solver.put("solver_type", "cg");
    t_solver.put("precon_type", "none");
    t_solver.put("error_tolerance", 1e-15);
    t_solver.put("max_iteration_step", 1000);
    t_root.put_child("LinearSolver", t_solver);



    int sparse_info[4] = {Example1::dim_eqs, Example1::dim_eqs, Example1::dim_eqs, Example1::dim_eqs};


    MathLib::PETScMatrix A;
    A.set_rank_size(mrank, msize);
    A.Init(Example1::dim_eqs, sparse_info);

    MathLib::PETScVector b(Example1::dim_eqs);
    MathLib::PETScVector x;
    b.CloneMe(x);

 
    MathLib::PETScLinearSolver petsc_leq(A, b, x);
    //petsc_leq.set_rank_size(mrank, msize);
    //x.set_rank_size(mrank, msize);
    //b.set_rank_size(mrank, msize);
    petsc_leq.Init(Example1::dim_eqs);

    checkLinearSolverInterface<MathLib::PETScLinearSolver, std::vector<double>>(petsc_leq, t_root);

}
#endif


#endif

