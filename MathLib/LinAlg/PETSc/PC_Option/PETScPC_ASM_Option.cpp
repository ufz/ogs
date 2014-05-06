/*!
   \file  PETScPC_ASM_Option.cpp
   \brief Define the configuration data for PETSc ASM preconditioner.

   \author Wenqing Wang
   \date 04-2014

   \copyright
    Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/

#include "PETScPC_ASM_Option.h"

namespace MathLib
{

using boost::property_tree::ptree;

PETScPC_ASM_Option::
PETScPC_ASM_Option(const boost::property_tree::ptree &option)
    : type(PC_ASM_BASIC)
{
    auto type_name = option.get_optional<std::string>("type");
    if(type_name->find("base") != std::string::npos)
    {
        type = PC_ASM_BASIC;
    }
    if(type_name->find("interpolate") != std::string::npos)
    {
        type = PC_ASM_INTERPOLATE;
    }
    if(type_name->find("restrict") != std::string::npos)
    {
        type = PC_ASM_RESTRICT;
    }
    if(type_name->find("PC_ASM_NONE") != std::string::npos)
    {
        type = PC_ASM_NONE;
    }
}

} // end namespace

