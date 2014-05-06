/*!
   \file  PETScLinearSolverOption.h
   \brief Define the configuration data for the PETSc linear solver.

   \author Wenqing Wang
   \date 02-2014

   \copyright
    Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/
#ifndef PETSCLINEARSOLVEROPTION_H_
#define PETSCLINEARSOLVEROPTION_H_

#include <string>
#include <boost/property_tree/ptree.hpp>

#include <petscksp.h>

namespace MathLib
{
/*!
   \brief Class to manage the basic configuration data for
          a Krylov subspace (KSP) interation method and a procondtioner (PC)
          of PETSc routines.
*/
class PETScLinearSolverOption
{
    public:

        /*!
            Default constructor
        */
        PETScLinearSolverOption() : _solver_name("bcgs"), _pc_name("bjacobi"),
                                    _preco_side(PC_LEFT),  _max_it(2000), 
                                    _rtol(1.e-5), _atol(PETSC_DEFAULT), _dtol(PETSC_DEFAULT)
        {}

        /*!
            \param ksp_option ptree of basic configuration data for solver
        */
        PETScLinearSolverOption(const boost::property_tree::ptree &ksp_option);

        /*!
            \param ksp_option ptree of basic configuration data for solver
            \param pc_option  ptree of basic configuration data for preconditioner
        */
        PETScLinearSolverOption(const boost::property_tree::ptree &ksp_option,
                                const boost::property_tree::ptree &pc_option);

        /// Set basic options for a KSP and a PC
        void setOption(KSP &ksp, PC &pc)
        {
            KSPSetType(ksp,  _solver_name.c_str());
            KSPGetPC(ksp, &pc);
            PCSetType(pc, _pc_name.c_str());
            KSPSetPCSide(ksp, _preco_side);

            KSPSetTolerances(ksp, _rtol, _atol, _dtol, _max_it);
        }

    private:
        /*!
            The name of solver, and it could be one of the following names
               "richardson"
               "chebychev"
               "cg"
               "cgne"
               "nash"
               "stcg"
               "gltr"
               "gmres"
               "fgmres"
               "lgmres"
               "dgmres"
               "tcqmr"
               "bcgs"
               "ibcgs"
               "bcgsl"
               "cgs"
               "tfqmr"
               "cr"
               "lsqr"
               "preonly"
               "qcg"
               "bicg"
               "minres"
               "symmlq"
               "lcd"
               "python"
               "broyden"
               "gcr"
               "ngmres"
               "specest"
        */
        std::string _solver_name;

        /*!
             The name of preconditioner, and it could be one of the following names
                none
                jacobi
                sor
                lu
                shell
                bjacobi
                mg
                eisenstat
                ilu
                icc
                asm
                gasm
                ksp
                composite
                redundant
                spai
                nn
                cholesky
                pbjacobi
                mat
                hypre
                parms
                fieldsplit
                tfs
                ml
                prometheus
                galerkin
                exotic
                hmpi
                supportgraph
                asa
                cp
                bfbt
                lsc
                python
                pfmg
                syspfmg
                redistribute
                sacusp
                sacusppoly
                bicgstabcusp
                svd
                ainvcusp
                gamg
        */
        std::string _pc_name;

        /// Flag for which side preconditioning.
        PCSide _preco_side;

        PetscInt _max_it; ///< Maximum iteration.

        PetscReal _rtol;  ///< Tolerance for the relative error, \f$e=|r|/|b|\f$.
        PetscReal _atol;  ///< Tolerance for the absolute error, \f$e=|r|\f$.
        PetscReal _dtol;  ///< Relative increase in the residual.
};

} // end namespace
#endif

