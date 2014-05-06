/*!
   \file  PETScPC_SOR_Option.h
   \brief Define the configuration data for PETSc SOR/SSOR preconditioner.

   \author Wenqing Wang
   \date 04-2014

   \copyright
    Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/
#ifndef PETSCPC_SOR_OPTION_H_
#define PETSCPC_SOR_OPTION_H_

#include <boost/property_tree/ptree.hpp>

#include <petscksp.h>

namespace MathLib
{
/*!
    PETSc SOR//SSOR preconditioner options

*/
struct PETScPC_SOR_Option
{
    PETScPC_SOR_Option(const boost::property_tree::ptree &option);

    /// Set SOR/SSOR option
    void setOption(PC &pc);

    /// Overloaded assign operator
    void operator = (const PETScPC_SOR_Option& opt)
    {
        omega = opt.omega;
        its = opt.its;
        lits = opt.lits;
        type = opt.type;
    }

    /// Relaxation number
    PetscReal omega;

    /// Number of parelllel iterations, each parallel iteration
    /// has 'lits' local iterations
    PetscInt its;

    /// Number of local iterations
    PetscInt lits;

    /*!
        SOR type:
        SOR_FORWARD_SWEEP
        SOR_BACKWARD_SWEEP
        SOR_SYMMETRIC_SWEEP
        SOR_LOCAL_FORWARD_SWEEP
        SOR_LOCAL_BACKWARD_SWEEP
        SOR_LOCAL_SYMMETRIC_SWEEP
    */
    MatSORType type;
};
} // end namespace
#endif

