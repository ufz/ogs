/*!
   \file  PETScPC_LU_Option.cpp
   \brief Define the configuration data for LU decompoistion preconditioner.

   \author Wenqing Wang
   \date 04-2014

   \copyright
    Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/

#include "PETScPC_LU_Option.h"

namespace MathLib
{

using boost::property_tree::ptree;

PETScPC_LU_Option::
PETScPC_LU_Option(const boost::property_tree::ptree &option)
    : mat_type(MATORDERINGNATURAL)
{
    auto type_name = option.get_optional<std::string>("type");
    if(type_name->find("natural") != std::string::npos)
    {
        mat_type = MATORDERINGNATURAL;
    }
    if(type_name->find("nd") != std::string::npos)
    {
        mat_type = MATORDERINGND;
    }
    if(type_name->find("1wd") != std::string::npos)
    {
        mat_type = MATORDERING1WD;
    }
    if(type_name->find("rcm") != std::string::npos)
    {
        mat_type = MATORDERINGRCM;
    }
    if(type_name->find("qmd") != std::string::npos)
    {
        mat_type = MATORDERINGQMD;
    }
    if(type_name->find("rowlength") != std::string::npos)
    {
        mat_type = MATORDERINGROWLENGTH;
    }
}

} // end namespace

