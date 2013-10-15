/**
 * \file
 * \author Norihiro Watanabe
 * \date   2013-05-15
 * \brief  Interface tests of global matrix classes
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <gtest/gtest.h>

#include "MathLib/LinAlg/Dense/DenseMatrix.h"
#include "MathLib/LinAlg/Dense/GlobalDenseMatrix.h"
#include "MathLib/LinAlg/FinalizeMatrixAssembly.h"
#ifdef USE_LIS
#include "MathLib/LinAlg/Lis/LisMatrix.h"
#endif

namespace
{

template <class T_MATRIX>
void checkGlobalMatrixInterface(T_MATRIX &m)
{
  
    ASSERT_EQ(10u, m.getNRows());
    ASSERT_EQ(10u, m.getNCols());
    //Cannot get a fixed number with PETSc for most of memory allocation types.
    //   ASSERT_EQ(0u,  m.getRangeBegin());
    //   ASSERT_EQ(10u, m.getRangeEnd());

    
    m.setValue(0, 0, 1.0);
    m.add(0, 0, 1.0);
    
    // PETSc need the following line.
    finalizeMatrixAssembly(m);

    m.setZero();

 
    MathLib::DenseMatrix<double> local_m(2,2, 1.0);
    std::vector<std::size_t> vec_pos(2);
    vec_pos[0] = 1;
    vec_pos[1] = 3;
    m.add(vec_pos, vec_pos, local_m);

    

    ASSERT_TRUE(finalizeMatrixAssembly(m));

 
}

} // end namespace

TEST(Math, CheckInterface_GlobalDenseMatrix)
{
    MathLib::GlobalDenseMatrix<double> m(10, 10);
    checkGlobalMatrixInterface(m);
}

#ifdef USE_LIS
TEST(Math, CheckInterface_LisMatrix)
{
    MathLib::LisMatrix m(10);
    checkGlobalMatrixInterface(m);
}
#endif

#ifdef USE_PETSC
#include "MathLib/LinAlg/PETSc/PETScMatrix.h"
TEST(Math, CheckInterface_PETScMatrix)
{

  int mrank, msize;

   MPI_Comm_rank(PETSC_COMM_WORLD, &mrank);
   MPI_Comm_size(PETSC_COMM_WORLD, &msize);

   if(msize != 3)
   {
      PetscSynchronizedPrintf(PETSC_COMM_WORLD, "===\nThis is test of PETSc matrix. The numnber of cores must be 3 exactly");

     PetscFinalize();
     exit(EXIT_FAILURE);
   }

   int dim = 10; 

   int sparse_info[4] = {dim, dim ,dim ,dim};


    MathLib::PETScMatrix m;
    m.set_rank_size(mrank, msize);
    m.Init(dim, sparse_info);

    checkGlobalMatrixInterface(m);
}
#endif
