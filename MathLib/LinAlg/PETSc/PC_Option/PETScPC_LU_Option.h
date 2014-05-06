/*!
   \file  PETScPC_LU_Option.h
   \brief Define the configuration data for LU decompoistion preconditioner.

   \author Wenqing Wang
   \date 04-2014

   \copyright
    Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/
#ifndef PETSSPC_LU_OPTION_H_
#define PETSSPC_LU_OPTION_H_

#include <boost/property_tree/ptree.hpp>

#include <petscksp.h>

namespace MathLib
{

/*!
    PETSc LU preconditioner options

*/
struct PETScPC_LU_Option
{
    PETScPC_LU_Option(const boost::property_tree::ptree &option);

    /// Set LU option
    void setOption(PC &pc)
    {
        PCFactorSetMatOrderingType(pc, mat_type);
    }

    /// Overloaded assign operator
    void operator = (const PETScPC_LU_Option& opt)
    {
        mat_type = opt.mat_type;
    }

    /// Ordering routine (to reduce fill) to be used in the LU factorization.
    MatOrderingType mat_type;
};

} // end namespace
#endif

