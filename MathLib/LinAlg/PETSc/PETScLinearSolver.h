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

#include "PETScMatrix.h"
#include "PETScVector.h"

#include "PETScLinearSolverOption.h"

#include "KSP_Option/PETScPC_KSP_Chebyshev_Option.h"
#include "KSP_Option/PETScPC_KSP_Richards_Option.h"
#include "KSP_Option/PETScPC_KSP_GMRES_Option.h"

#include "PC_Option/PETScPC_ILU_Option.h"
#include "PC_Option/PETScPC_SOR_Option.h"
#include "PC_Option/PETScPC_LU_Option.h"
#include "PC_Option/PETScPC_ASM_Option.h"
#include "PC_Option/PETScPC_AMG_Option.h"

namespace MathLib
{

/*!
     A class of linear solver based on PETSc rountines.

*/
class PETScLinearSolver
{
    public:

        /*!
            Constructor
            \param A       Matrix, cannot be constant.
            \param option  Configuration data for solver and preconditioner.
        */
        PETScLinearSolver(PETScMatrix &A, const boost::property_tree::ptree &option);

        template <typename T_KSP_OPTION> void setKSP_Option(T_KSP_OPTION &ksp_opt)
        {
            ksp_opt.setOption(_solver);
        }

        template <typename T_PC_OPTION> void setPC_Option(T_PC_OPTION &pc_opt)
        {
            pc_opt.setOption(_pc);
        }

        ~PETScLinearSolver()
        {
            KSPDestroy(&_solver);
        }

        /*!
            Solve a system of equations.
            \param b The right hand of the equations.
            \param x The solutions to be solved.
        */
        void solve(const PETScVector &b, PETScVector &x);

    private:
        KSP _solver; ///< Slover type.
        PC _pc;      ///< Preconditioner type.
};

} // end namespace
#endif

