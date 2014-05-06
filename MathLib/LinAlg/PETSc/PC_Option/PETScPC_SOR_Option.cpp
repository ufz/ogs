/*!
   \file  PETScPC_SOR_Option.cpp
   \brief Define the configuration data for PETSc SOR/SSOR preconditioner.

   \author Wenqing Wang
   \date 04-2014

   \copyright
    Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/

#include "PETScPC_SOR_Option.h"

namespace MathLib
{

using boost::property_tree::ptree;

PETScPC_SOR_Option::
PETScPC_SOR_Option(const boost::property_tree::ptree &option)
    : omega(1.), its(PETSC_DEFAULT),
      lits(PETSC_DEFAULT), type(SOR_FORWARD_SWEEP)
{
    auto val = option.get_optional<double>("omega");
    omega	= *val;

    auto n_its = option.get_optional<int>("local_iterations");
    lits = *n_its;

    n_its = option.get_optional<int>("parallel_iterations");
    its = *n_its;

    auto type_name = option.get_optional<std::string>("type");
    if(type_name->find("pc_sor_symmetric") != std::string::npos)
    {
        type = SOR_FORWARD_SWEEP;
    }
    if(type_name->find("pc_sor_backward") != std::string::npos)
    {
        type = SOR_BACKWARD_SWEEP;
    }
    if(type_name->find("pc_sor_local_forward") != std::string::npos)
    {
        type = SOR_LOCAL_FORWARD_SWEEP;
    }
    if(type_name->find("pc_sor_local_symmetric") != std::string::npos)
    {
        type = SOR_LOCAL_SYMMETRIC_SWEEP;
    }
    if(type_name->find("pc_sor_local_backward") != std::string::npos)
    {
        type = SOR_LOCAL_BACKWARD_SWEEP;
    }
}

void PETScPC_SOR_Option::setOption(PC &pc)
{
    PCSORSetOmega(pc, omega);
    PCSORSetIterations(pc, its, lits);
    PCSORSetSymmetric(pc, type);
}

} // end namespace

