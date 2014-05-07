/*!
   \file  PETScKSP_Chebyshev_Option.cpp
   \brief Define the configuration data for the PETSc Chebyshev linear solver.

   \author Wenqing Wang
   \date 04-2014

   \copyright
    Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/

#include "PETScKSP_Chebyshev_Option.h"

namespace MathLib
{

using boost::property_tree::ptree;

PETScKSP_Chebyshev_Option::
PETScKSP_Chebyshev_Option(const boost::property_tree::ptree &option)
    : emin_chebyshev(0.01), emax_chebyshev(100.0)
{
    auto val = option.get_optional<double>("smallest_eignvalue");
    emin_chebyshev = *val;

    val = option.get_optional<double>("maximum_eignvalue");
    emax_chebyshev = *val;
}

} // end namespace

