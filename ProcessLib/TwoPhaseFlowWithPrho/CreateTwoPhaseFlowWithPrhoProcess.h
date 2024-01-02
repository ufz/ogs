/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
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
    std::string const& name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config,
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media);
}  // namespace TwoPhaseFlowWithPrho
}  // namespace ProcessLib
