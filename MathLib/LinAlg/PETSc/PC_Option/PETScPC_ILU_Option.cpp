/*!
   \file  PETScPC_ILU_Option.cpp
   \brief Define the configuration data for PETSc ILU preconditioner.

   \author Wenqing Wang
   \date 04-2014

   \copyright
    Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/

#include "PETScPC_ILU_Option.h"

namespace MathLib
{

using boost::property_tree::ptree;

PETScPC_ILU_Option::
PETScPC_ILU_Option(const boost::property_tree::ptree &option)
    : levels(PETSC_DECIDE), reuse_ordering(false),
      reuse_fill(false), use_in_place(false), allow_diagonal_fill(false)
{
    auto val = option.get_optional<double>("levels");
    levels	= *val;

    auto reuse_order = option.get_optional<bool>("reuse_ordering");
    reuse_ordering = *reuse_order;

    auto rfill = option.get_optional<bool>("reuse_fill");
    reuse_fill = *rfill;

    auto inplane = option.get_optional<bool>("use_in_place");
    use_in_place = *inplane;

    auto low_diag_fill = option.get_optional<bool>("allow_diagonal_fill");
    allow_diagonal_fill = *low_diag_fill;
}

void PETScPC_ILU_Option::setOption(PC &pc)
{
    PCFactorSetLevels(pc, levels);

    if(reuse_ordering)
    {
        PCFactorSetReuseOrdering(pc, PETSC_TRUE);
    }

    if(reuse_fill)
    {
        PCFactorSetReuseFill(pc, PETSC_TRUE);
    }

    if(use_in_place)
    {
        PCFactorSetUseInPlace(pc);
    }

    if(allow_diagonal_fill)
    {
        PCFactorSetAllowDiagonalFill(pc);
    }
}

} // end namespace

