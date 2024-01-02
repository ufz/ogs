/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateConstitutiveSetting.h"

#include "MaterialLib/SolidModels/CreateConstitutiveRelationsGeneric.h"
#include "MaterialLib/SolidModels/MFront/CreateMFrontGeneric.h"

namespace ProcessLib::ThermoRichardsMechanics
{
namespace ConstitutiveStressSaturation_StrainPressureTemperature
{
template <int DisplacementDim>
std::unique_ptr<SolidConstitutiveRelation<DisplacementDim>> createMFrontGeneric(
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    std::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    BaseLib::ConfigTree const& config)
{
    namespace MSM = MaterialLib::Solids::MFront;
    using namespace boost::mp11;

    bool const library_path_is_relative_to_prj_file = true;

    return MSM::createMFrontGeneric<
        DisplacementDim, mp_list<MSM::Strain, MSM::LiquidPressure>,
        mp_list<MSM::Stress, MSM::Saturation>, mp_list<MSM::Temperature>>(
        parameters, local_coordinate_system, config,
        library_path_is_relative_to_prj_file);
}

template <int DisplacementDim>
std::map<int, std::unique_ptr<SolidConstitutiveRelation<DisplacementDim>>>
CreateConstitutiveSetting<DisplacementDim>::createSolidConstitutiveRelations(
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    std::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    BaseLib::ConfigTree const& config)
{
    return MaterialLib::Solids::createConstitutiveRelationsGeneric(
        parameters, local_coordinate_system, config,
        createMFrontGeneric<DisplacementDim>);
}

template struct CreateConstitutiveSetting<2>;
template struct CreateConstitutiveSetting<3>;
}  // namespace ConstitutiveStressSaturation_StrainPressureTemperature
}  // namespace ProcessLib::ThermoRichardsMechanics
