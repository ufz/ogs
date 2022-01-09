/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <map>
#include <memory>
#include <string>
#include <vector>

namespace BaseLib
{
class ConfigTree;
}
namespace MaterialPropertyLib
{
class Medium;
}
namespace MathLib
{
class PiecewiseLinearInterpolation;
}
namespace MeshLib
{
class Mesh;
}
namespace ParameterLib
{
struct ParameterBase;
}
namespace ProcessLib
{
class AbstractJacobianAssembler;
}
namespace ProcessLib
{
class Process;
}
namespace ProcessLib
{
class ProcessVariable;
}

namespace ProcessLib
{
namespace TwoPhaseFlowWithPrho
{
std::unique_ptr<Process> createTwoPhaseFlowWithPrhoProcess(
    std::string name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves);
}  // namespace TwoPhaseFlowWithPrho
}  // namespace ProcessLib
