/*!
   \file  PETScLinearSolver.h
   \brief Declaration of class PETScLinearSolver, which defines a solver object
         based on PETSc routines.

   \author Wenqing Wang
   \version
   \date Nov 2011 - Sep 2013

   \copyright
   Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/

#ifndef PETSCLINEARSOLVER_H_
#define PETSCLINEARSOLVER_H_

#include<string>

#include <petscksp.h>

#include <logog/include/logog.hpp>

#include "BaseLib/ConfigTree.h"

#include "PETScMatrix.h"
#include "PETScVector.h"

namespace MathLib
{

/*! A class of linear solver based on PETSc rountines.

 All command-line options that are not recognized by OGS are passed on to
 PETSc, i.e., they potentially affect the linear solver.
 The linear solver options in the project file take precedence over the
 command-line options, because the former are processed at a later stage.
*/
class PETScLinearSolver
{
    public:

        /*!
            Constructor
            \param prefix  Name used to distinguish the options in the command
                           line for this solver. It can be the name of the PDE
                           that owns an instance of this class.
            \param option  Petsc options, which will be inserted into the global
                           petsc options database.
        */
        PETScLinearSolver(const std::string prefix,
                          BaseLib::ConfigTree const* const option);

        ~PETScLinearSolver()
        {
            KSPDestroy(&_solver);
        }

        // TODO check if some args in LinearSolver interface can be made const&.
        bool solve(PETScMatrix& A, PETScVector &b, PETScVector &x);

        /// Get number of iterations.
        PetscInt getNumberOfIterations() const
        {
            PetscInt its = 0;
            KSPGetIterationNumber(_solver, &its);
            return its;
        }

        /// Get elapsed wall clock time.
        double getElapsedTime() const
        {
            return _elapsed_ctime;
        }

    private:
        KSP _solver; ///< Solver type.
        PC _pc;      ///< Preconditioner type.

        double _elapsed_ctime = 0.0; ///< Clock time
};

} // end namespace
#endif

