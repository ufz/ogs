/*!
   \file  PETScPC_AMG_Option.h
   \brief Define the configuration data for PETSc AMG preconditioner.

   \author Wenqing Wang
   \date 04-2014

   \copyright
    Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/
#ifndef PETSCPC_AMG_OPTION_H_
#define PETSCPC_AMG_OPTION_H_

#include <boost/property_tree/ptree.hpp>

#include <petscksp.h>

namespace MathLib
{
/*!
    PETSc AMG preconditioner options

*/
struct PETScPC_AMG_Option
{
    PETScPC_AMG_Option(const boost::property_tree::ptree &option);

    /// Set AMG option
    void setOption(PC &pc)
    {
        PCGAMGSetType(pc, type);
        if(is_agg)
        {
            PCGAMGSetNSmooths(pc, nsmooths);
        }
    }

    /// Overloaded assign operator
    void operator = (const PETScPC_AMG_Option& opt)
    {
        nsmooths = opt.nsmooths;
        type = opt.type;
    }

    /// Number of smoothing steps
    PetscInt nsmooths;

    /*!
        AMG type: PCGAMGAGG or PCGAMGGEO
    */
    PCGAMGType type;

    /// Flag of agg type
    bool is_agg;
};
} // end namespace
#endif

