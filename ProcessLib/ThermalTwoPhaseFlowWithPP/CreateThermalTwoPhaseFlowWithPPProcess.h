// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <memory>

#include "MaterialLib/MPL/Medium.h"
#include "ProcessLib/Process.h"

namespace ProcessLib
{
namespace ThermalTwoPhaseFlowWithPP
{
std::unique_ptr<Process> createThermalTwoPhaseFlowWithPPProcess(
    std::string const& name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config,
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media);
}  // namespace ThermalTwoPhaseFlowWithPP
}  // namespace ProcessLib
