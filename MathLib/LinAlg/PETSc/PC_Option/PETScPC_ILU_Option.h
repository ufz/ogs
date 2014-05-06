/*!
   \file  PETScPC_ILU_Option.h
   \brief Define the configuration data for PETSc ILU preconditioner.

   \author Wenqing Wang
   \date 04-2014

   \copyright
    Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/
#ifndef PETSCPC_ILU_OPTION_H_
#define PETSCPC_ILU_OPTION_H_

#include <boost/property_tree/ptree.hpp>

#include <petscksp.h>

namespace MathLib
{

/*!
    PETSc ILU preconditioner options

 */
struct PETScPC_ILU_Option
{
    PETScPC_ILU_Option(const boost::property_tree::ptree &option);

    /// Set ILU option
    void setOption(PC &pc);

    /// Overloaded assign operator
    void operator = (const PETScPC_ILU_Option& opt)
    {
        levels = opt.levels;
        reuse_ordering = opt.reuse_ordering;
        reuse_fill = opt.reuse_fill;
        use_in_place = opt.use_in_place;
        allow_diagonal_fill = opt.allow_diagonal_fill;
    }

    /// Number of levels of fill for ILU(k)
    int levels;

    /// Reuse ordering of factorized matrix from previous factorization
    bool reuse_ordering;

    /// Use the fill ratio computed in the initial factorization.
    bool reuse_fill;

    /*! for ILU(0) with natural ordering, reuses the space of the matrix
        for its factorization (overwrites original matrix)
    */
    bool use_in_place;

    /// fill in a zero diagonal even if levels of fill indicate it wouldn't be fill
    bool allow_diagonal_fill;
};

} // end namespace
#endif

