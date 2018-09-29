/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateCohesiveZoneModeI.h"

#include "CohesiveZoneModeI.h"
#include "ProcessLib/Utils/ProcessUtils.h"  // required for findParameter

namespace MaterialLib
{
namespace Fracture
{
namespace CohesiveZoneModeI
{
template <int DisplacementDim>
std::unique_ptr<FractureModelBase<DisplacementDim>> createCohesiveZoneModeI(
    std::vector<std::unique_ptr<ProcessLib::ParameterBase>> const& parameters,
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{material__fracture_model__type}
    config.checkConfigParameter("type", "CohesiveZoneModeI");
    DBUG("Create CohesiveZoneModeI material");

    auto& Kn = ProcessLib::findParameter<double>(
        //! \ogs_file_param_special{material__fracture_model__CohesiveZoneModeI__normal_stiffness}
        config, "normal_stiffness", parameters, 1);

    auto& Ks = ProcessLib::findParameter<double>(
        //! \ogs_file_param_special{material__fracture_model__CohesiveZoneModeI__shear_stiffness}
        config, "shear_stiffness", parameters, 1);

    auto& Gc = ProcessLib::findParameter<double>(
        //! \ogs_file_param_special{material__fracture_model__CohesiveZoneModeI__fracture_toughness}
        config, "fracture_toughness", parameters, 1);

    auto& t_np = ProcessLib::findParameter<double>(
        //! \ogs_file_param_special{material__fracture_model__CohesiveZoneModeI__peak_normal_traction}
        config, "peak_normal_traction", parameters, 1);

    auto const penalty_aperture_cutoff =
        //! \ogs_file_param{material__fracture_model__CohesiveZoneModeI__penalty_aperture_cutoff}
        config.getConfigParameter<double>("penalty_aperture_cutoff");

    auto const tension_cutoff =
        //! \ogs_file_param{material__fracture_model__CohesiveZoneModeI__tension_cutoff}
        config.getConfigParameter<bool>("tension_cutoff");

    MaterialPropertiesParameters mp{Kn, Ks, Gc, t_np};

    return std::make_unique<CohesiveZoneModeI<DisplacementDim>>(
        penalty_aperture_cutoff, tension_cutoff, mp);
}

template std::unique_ptr<FractureModelBase<2>> createCohesiveZoneModeI(
    std::vector<std::unique_ptr<ProcessLib::ParameterBase>> const& parameters,
    BaseLib::ConfigTree const& config);

template std::unique_ptr<FractureModelBase<3>> createCohesiveZoneModeI(
    std::vector<std::unique_ptr<ProcessLib::ParameterBase>> const& parameters,
    BaseLib::ConfigTree const& config);

}  // namespace CohesiveZoneModeI
}  // namespace Fracture
}  // namespace MaterialLib
