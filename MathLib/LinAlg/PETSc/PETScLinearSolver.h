/*!
   \file  PETScLinearSolver.h
   \brief Declaration of class PETScLinearSolver, which defines a solver object
         based on PETSc routines.

   \author Wenqing Wang
   \version
   \date Nov 2011 - Sep 2013

   \copyright
   Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/

#ifndef PETSCLINEARSOLVER_H_
#define PETSCLINEARSOLVER_H_

#include<string>

#include "petscksp.h"

#include "PETScMatrix.h"
#include "PETScVector.h"

namespace MathLib
{

/*!
     A class of linear solver based on PETSc rountines.

     All options for KSP and PC must given in the command line.
*/
class PETScLinearSolver
{
    public:

        /*!
            Constructor
            \param A       Matrix, cannot be constant.
            \param prefix  Name used to distinguish the options in the command
                           line for this solver. It can be the name of the PDE
                           that owns an instance of this class.
        */
        PETScLinearSolver(PETScMatrix &A, const std::string &prefix="");

        ~PETScLinearSolver()
        {
            KSPDestroy(&_solver);
        }

        /*!
            Solve a system of equations.
            \param b The right hand of the equations.
            \param x The solutions to be solved.
            \return  true: converged, false: diverged due to exceeding
                     the maximum iterations.
        */
        bool solve(const PETScVector &b, PETScVector &x);

        /*!
            \brief Get number of iterations.
            If function solve(...) returns false, the return value is
            exactly the maximum iterations.
        */
        PetscInt getNumberOfIterations() const
        {
            PetscInt its = 0;
            KSPGetIterationNumber(_solver, &its);
            return its;
        }

    private:
        /// Matrix, kept as a member only for solving successive linear equation
        /// that preconditioner matrix may vary.
        PETScMatrix &_A;

        KSP _solver; ///< Solver type.
        PC _pc;      ///< Preconditioner type.
};

} // end namespace
#endif

