/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   createLiquidFlowProcess.h
 *
 * Created on August 19, 2016, 1:30 PM
 */

#ifndef CREATELIQUIDFLOWPROCESS_H
#define CREATELIQUIDFLOWPROCESS_H

#include <memory>
#include "ProcessLib/Process.h"

namespace ProcessLib
{
namespace LiquidFlow
{
std::unique_ptr<Process> createLiquidFlowProcess(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& /*parameters*/,
    BaseLib::ConfigTree const& config);
}  // end of namespace
}  // end of namespace

#endif /* CREATELIQUIDFLOWPROCESS_H */
