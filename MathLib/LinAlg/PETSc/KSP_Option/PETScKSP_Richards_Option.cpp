/*!
   \file  PETScKSP_Richards_Option.cpp
   \brief Define the configuration data for the PETSc Richards linear solver.

   \author Wenqing Wang
   \date 04-2014

   \copyright
    Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/

#include "PETScKSP_Richards_Option.h"

namespace MathLib
{

using boost::property_tree::ptree;

PETScKSP_Richards_Option::
PETScKSP_Richards_Option(const boost::property_tree::ptree &option)
    : damping_factor_richards(1.0)
{
    auto damping_factor = option.get_optional<double>("damping_factor");
    damping_factor_richards = *damping_factor;
}

} // end namespace

