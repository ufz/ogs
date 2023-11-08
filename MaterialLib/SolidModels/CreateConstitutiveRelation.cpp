/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  Created on July 10, 2018, 12:09 PM
 */

#include "CreateConstitutiveRelation.h"

#include "CreateCreepBGRa.h"
#include "CreateEhlers.h"
#include "CreateLinearElasticIsotropic.h"
#include "CreateLinearElasticOrthotropic.h"
#include "CreateLinearElasticTransverseIsotropic.h"
#include "CreateLubby2.h"
#ifdef OGS_USE_MFRONT
#include "MFront/CreateMFront.h"
#endif  // OGS_USE_MFRONT

#include "CreateConstitutiveRelationsGeneric.h"
#include "MechanicsBase.h"

namespace MaterialLib
{
namespace Solids
{
template <int DisplacementDim>
std::unique_ptr<MaterialLib::Solids::MechanicsBase<DisplacementDim>>
createConstitutiveRelation(
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    std::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    BaseLib::ConfigTree const& config)
{
    auto const type =
        //! \ogs_file_param{material__solid__constitutive_relation__type}
        config.peekConfigParameter<std::string>("type");

    if (type == "Ehlers")
    {
        return MaterialLib::Solids::Ehlers::createEhlers<DisplacementDim>(
            parameters, config);
    }
    if (type == "LinearElasticIsotropic")
    {
        const bool skip_type_checking = false;
        return MaterialLib::Solids::createLinearElasticIsotropic<
            DisplacementDim>(parameters, config, skip_type_checking);
    }
    if (type == "LinearElasticOrthotropic")
    {
        const bool skip_type_checking = false;
        return MaterialLib::Solids::createLinearElasticOrthotropic<
            DisplacementDim>(
            parameters, local_coordinate_system, config, skip_type_checking);
    }
    if (type == "LinearElasticTransverseIsotropic")
    {
        const bool skip_type_checking = false;
        return MaterialLib::Solids::createLinearElasticTransverseIsotropic<
            DisplacementDim>(
            parameters, local_coordinate_system, config, skip_type_checking);
    }
    if (type == "Lubby2")
    {
        return MaterialLib::Solids::Lubby2::createLubby2<DisplacementDim>(
            parameters, config);
    }
    if (type == "CreepBGRa")
    {
        return MaterialLib::Solids::Creep::createCreepBGRa<DisplacementDim>(
            parameters, config);
    }
    if (type == "MFront")
    {
#ifdef OGS_USE_MFRONT
        return MaterialLib::Solids::MFront::createMFront<DisplacementDim>(
            parameters, local_coordinate_system, config);
#else   // OGS_USE_MFRONT
        OGS_FATAL(
            "OGS is compiled without MFront support. See OGS_USE_MFRONT CMake "
            "option.");
#endif  // OGS_USE_MFRONT
    }
    OGS_FATAL("Cannot construct constitutive relation of given type '{:s}'.",
              type);
}

template std::unique_ptr<MaterialLib::Solids::MechanicsBase<2>>
createConstitutiveRelation(
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    std::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    BaseLib::ConfigTree const& config);

template std::unique_ptr<MaterialLib::Solids::MechanicsBase<3>>
createConstitutiveRelation(
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    std::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    BaseLib::ConfigTree const& config);

template <int DisplacementDim>
std::map<int,
         std::unique_ptr<MaterialLib::Solids::MechanicsBase<DisplacementDim>>>
createConstitutiveRelations(
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    std::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    BaseLib::ConfigTree const& config)
{
    return createConstitutiveRelationsGeneric<
        MaterialLib::Solids::MechanicsBase<DisplacementDim>>(
        parameters,
        local_coordinate_system,
        config,
        createConstitutiveRelation<DisplacementDim>);
}

template std::map<int, std::unique_ptr<MaterialLib::Solids::MechanicsBase<2>>>
createConstitutiveRelations<2>(
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    std::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    BaseLib::ConfigTree const& config);

template std::map<int, std::unique_ptr<MaterialLib::Solids::MechanicsBase<3>>>
createConstitutiveRelations<3>(
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    std::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    BaseLib::ConfigTree const& config);
}  // namespace Solids
}  // namespace MaterialLib
