/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateMohrCoulomb.h"

#include "ProcessLib/Utils/ProcessUtils.h"  // required for findParameter
#include "MohrCoulomb.h"

namespace MaterialLib
{
namespace Fracture
{

template <int DisplacementDim>
std::unique_ptr<FractureModelBase<DisplacementDim>>
createMohrCoulomb(
    std::vector<std::unique_ptr<ProcessLib::ParameterBase>> const& parameters,
    BaseLib::ConfigTree const& config)
{
    config.checkConfigParameter("type", "MohrCoulomb");
    DBUG("Create MohrCoulomb material");

    auto& Kn = ProcessLib::findParameter<double>(
        config, "normal_stiffness", parameters, 1);

    auto& Ks = ProcessLib::findParameter<double>(
        config, "shear_stiffness", parameters, 1);

    auto& friction_angle = ProcessLib::findParameter<double>(
        config, "friction_angle", parameters, 1);

    auto& dilatancy_angle = ProcessLib::findParameter<double>(
        config, "dilatancy_angle", parameters, 1);

    auto& cohesion = ProcessLib::findParameter<double>(
        config, "cohesion", parameters, 1);


    typename MohrCoulomb<DisplacementDim>::MaterialProperties mp{
        Kn, Ks, friction_angle, dilatancy_angle, cohesion};

    return std::unique_ptr<MohrCoulomb<DisplacementDim>>{
        new MohrCoulomb<DisplacementDim>{mp}};
}


template
std::unique_ptr<FractureModelBase<2>>
createMohrCoulomb(
    std::vector<std::unique_ptr<ProcessLib::ParameterBase>> const& parameters,
    BaseLib::ConfigTree const& config);

template
std::unique_ptr<FractureModelBase<3>>
createMohrCoulomb(
    std::vector<std::unique_ptr<ProcessLib::ParameterBase>> const& parameters,
    BaseLib::ConfigTree const& config);

}  // namespace Fracture
}  // namespace MaterialLib
