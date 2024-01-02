/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "CreateMFront.h"

#ifndef _WIN32
#include <dlfcn.h>
#endif

#include "CreateMFrontGeneric.h"
#include "MFront.h"

namespace MaterialLib
{
namespace Solids
{
namespace MFront
{
template <int DisplacementDim>
std::unique_ptr<MechanicsBase<DisplacementDim>> createMFront(
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    std::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    BaseLib::ConfigTree const& config)
{
    bool const library_path_is_relative_to_prj_file = true;
    auto conf = createMFrontConfig(DisplacementDim, parameters, config,
                                   library_path_is_relative_to_prj_file);

    return std::make_unique<MFront<DisplacementDim>>(
        std::move(conf.behaviour), std::move(conf.material_properties),
        std::move(conf.state_variables_initial_properties),
        local_coordinate_system);
}
}  // namespace MFront
}  // namespace Solids
}  // namespace MaterialLib

namespace MaterialLib
{
namespace Solids
{
namespace MFront
{
template std::unique_ptr<MechanicsBase<2>> createMFront<2>(
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    std::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    BaseLib::ConfigTree const& config);
template std::unique_ptr<MechanicsBase<3>> createMFront<3>(
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    std::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    BaseLib::ConfigTree const& config);
}  // namespace MFront
}  // namespace Solids
}  // namespace MaterialLib
