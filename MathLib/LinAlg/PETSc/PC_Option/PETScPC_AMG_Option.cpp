/*!
   \file  PETScPC_AMG_Option.cpp
   \brief Define the configuration data for PETSc AMG preconditioner.

   \author Wenqing Wang
   \date 04-2014

   \copyright
    Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/

#include "PETScPC_AMG_Option.h"

namespace MathLib
{

using boost::property_tree::ptree;

PETScPC_AMG_Option::
PETScPC_AMG_Option(const boost::property_tree::ptree &option)
    : nsmooths(1), type(PCGAMGAGG), is_agg(true)
{
    auto val = option.get_optional<double>("agg_nsmooths");
    nsmooths = *val;

    auto type_name = option.get_optional<std::string>("type");
    if(type_name->find("agg") != std::string::npos)
    {
        type = PCGAMGAGG;
    }
    if(type_name->find("geo") != std::string::npos)
    {
        type = PCGAMGGEO;
        is_agg = false;
    }
}

} // end namespace

