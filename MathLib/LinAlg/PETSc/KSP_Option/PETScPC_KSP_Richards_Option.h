/*!
   \file  PETScPC_KSP_Richards_Option.h
   \brief Define the configuration data for the PETSc Richards linear solver.

   \author Wenqing Wang
   \date 04-2014

   \copyright
    Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/
#ifndef PETSCPC_KSP_RICHARDS_OPTION_H_
#define PETSCPC_KSP_RICHARDS_OPTION_H_

#include <boost/property_tree/ptree.hpp>

#include <petscksp.h>

namespace MathLib
{

/*!
    PETSc KSP Richards options

 */
struct PETScPC_KSP_Richards_Option
{
    PETScPC_KSP_Richards_Option(const boost::property_tree::ptree &option);

    /// Set Richards option
    void setOption(KSP &solver)
    {
        KSPRichardsonSetScale(solver, damping_factor_richards);
    }

    /// Overloaded assign operator
    void operator = (const PETScPC_KSP_Richards_Option& opt)
    {
        damping_factor_richards = opt.damping_factor_richards;
    }

    /// Damping factor for Richards.
    PetscReal damping_factor_richards;
};

} // end namespace
#endif

