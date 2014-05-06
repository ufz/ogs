/*!
   \file  PETScPC_ASM_Option.h
   \brief Define the configuration data for PETSc ASM preconditioner.

   \author Wenqing Wang
   \date 04-2014

   \copyright
    Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/
#ifndef PETSCPC_ASM_OPTION_H_
#define PETSCPC_ASM_OPTION_H_

#include <boost/property_tree/ptree.hpp>

#include <petscksp.h>

namespace MathLib
{
/*!
    PETSc ASM preconditioner options

*/
struct PETScPC_ASM_Option
{
    PETScPC_ASM_Option(const boost::property_tree::ptree &option);

    /// Set ASM option
    void setOption(PC &pc)
    {
        PCASMSetType(pc, type);
    }

    /// Overloaded assign operator
    void operator = (const PETScPC_ASM_Option& opt)
    {
        type = opt.type;
    }

    /*!
        ASM type:
        PC_ASM_BASIC,
        PC_ASM_INTERPOLATE,
        PC_ASM_RESTRICT,
        PC_ASM_NONE
    */
    PCASMType type;
};
} // end namespace
#endif

