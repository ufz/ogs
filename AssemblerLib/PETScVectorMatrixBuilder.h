/*!
  \file PETScVectorMatrixBuilder.h
  \author Wenqing Wang
  \date   2014.11
  \brief  Interface to create PETSc matrix and vector

  \copyright
  Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
             Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license

*/

#ifndef ASSEMBLERLIB_PETSCVECTORMATRIXBUILDER_H_
#define ASSEMBLERLIB_PETSCVECTORMATRIXBUILDER_H_

#include <petscmat.h>

#include "MathLib/LinAlg/PETSc/PETScMatrixOption.h"

#include "AssemblerLib/MeshComponentMap.h"

namespace AssemblerLib
{

template <typename MatrixType_, typename VectorType_>
class PETScVectorMatrixBuilder
{
public:
    typedef VectorType_ VectorType;
    typedef MatrixType_ MatrixType;

public:
    /*!
        \param vec_size       The size of the vector, either global or local
        \param is_global_size The flag of the type of vec_size, i.e. whether it is a global size
                              or local size. The default is true.
    */
    static
    VectorType* createVector(const PetscInt vec_size, const bool is_global_size = true)
    {
        VectorType* vec = new VectorType(vec_size, is_global_size);
        return vec;
    }

    /*!
      \param nrows  The number of rows of the matrix or the local matrix.
      \param mat_op The configuration information for creating a matrix.
    */
    static
    MatrixType* createMatrix(const PetscInt nrows, 
                  const MathLib::PETScMatrixOption &mat_op = MathLib::PETScMatrixOption())
    {
        MatrixType* mat = new MatrixType(nrows, mat_op);
        return mat;
    }
    
};

}   // namespace AssemblerLib

#endif  // ASSEMBLERLIB_PETSCVECTORMATRIXBUILDER_H_
