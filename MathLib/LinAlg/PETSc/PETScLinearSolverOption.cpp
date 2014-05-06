/*!
   \file  PETScLinearSolverOption.cpp
   \brief Define members of PETScLinearSolverOption

   \author Wenqing Wang
   \date 02-2014

   \copyright
    Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/

#include "PETScLinearSolverOption.h"

namespace MathLib
{

using boost::property_tree::ptree;

PETScLinearSolverOption::
PETScLinearSolverOption(const boost::property_tree::ptree &ksp_option)
    : _solver_name("bcgs"), _pc_name("bjacobi"), _preco_side(PC_LEFT),
      _max_it(2000), _rtol(1.e-5), _atol(PETSC_DEFAULT), _dtol(PETSC_DEFAULT)
{
    auto solver_type = ksp_option.get_optional<std::string>("solver_type");       
    if(solver_type)
    {
        _solver_name = *solver_type;
    }

    auto max_iteration_step = ksp_option.get_optional<int>("max_it");
    if(max_iteration_step)
    {
        _max_it = *max_iteration_step;
    }

    auto error_tolerance = ksp_option.get_optional<double>("rtol");
    if(error_tolerance)
    {
        _rtol = *error_tolerance;
    }

    error_tolerance = ksp_option.get_optional<double>("atol");
    if(error_tolerance)
    {
        _atol = *error_tolerance;
    }

    error_tolerance = ksp_option.get_optional<double>("dtol");
    if(error_tolerance)
    {
        _dtol = *error_tolerance;
    }
}

PETScLinearSolverOption::
PETScLinearSolverOption(const boost::property_tree::ptree &ksp_option,
                        const boost::property_tree::ptree &pc_option)
    : _solver_name("bcgs"), _pc_name("bjacobi"), _preco_side(PC_LEFT),
      _max_it(2000), _rtol(1.e-5), _atol(PETSC_DEFAULT), _dtol(PETSC_DEFAULT)
{
    auto solver_type = ksp_option.get_optional<std::string>("solver_type");       
    if(solver_type)
    {
        _solver_name = *solver_type;
    }

    auto max_iteration_step = ksp_option.get_optional<int>("max_it");
    if(max_iteration_step)
    {
        _max_it = *max_iteration_step;
    }

    auto error_tolerance = ksp_option.get_optional<double>("rtol");
    if(error_tolerance)
    {
        _rtol = *error_tolerance;
    }

    error_tolerance = ksp_option.get_optional<double>("atol");
    if(error_tolerance)
    {
        _atol = *error_tolerance;
    }

    error_tolerance = ksp_option.get_optional<double>("dtol");
    if(error_tolerance)
    {
        _dtol = *error_tolerance;
    }

    // Preconditioners:
    auto pc_type = pc_option.get_optional<std::string>("pc_type");
    if(pc_type)
    {
        _pc_name = *pc_type;
    }
    
    auto pc_side = pc_option.get_optional<std::string>("pc_side");
    if(pc_side)
    {
        if(pc_side->find("left") != std::string::npos)
            _preco_side = PC_LEFT;
        if(pc_side->find("right") != std::string::npos)
            _preco_side = PC_RIGHT;
        if(pc_side->find("symmetric") != std::string::npos)
            _preco_side = PC_SYMMETRIC;
    }
}

} // end namespace

