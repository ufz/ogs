/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file   CreateConstitutiveRelation.cpp
 *  Created on July 10, 2018, 12:09 PM
 */

#include "CreateConstitutiveRelation.h"
#include "CreateEhlers.h"
#include "CreateLinearElasticIsotropic.h"
#include "CreateLubby2.h"

#include "MechanicsBase.h"

#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"

#include "ProcessLib/Parameter/Parameter.h"

namespace MaterialLib
{
namespace Solids
{
template <int DisplacementDim>
std::unique_ptr<MaterialLib::Solids::MechanicsBase<DisplacementDim>>
createConstitutiveRelation(
    std::vector<std::unique_ptr<ProcessLib::ParameterBase>> const& parameters,
    BaseLib::ConfigTree const& config)
{
    auto const constitutive_relation_config =
        //! \ogs_file_param{material__solid__constitutive_relation}
        config.getConfigSubtree("constitutive_relation");

    auto const type =
        //! \ogs_file_param{material__solid__constitutive_relation__type}
        constitutive_relation_config.peekConfigParameter<std::string>("type");

    if (type == "Ehlers")
    {
        return MaterialLib::Solids::Ehlers::createEhlers<DisplacementDim>(
            parameters, constitutive_relation_config);
    }
    else if (type == "LinearElasticIsotropic")
    {
        return MaterialLib::Solids::createLinearElasticIsotropic<
            DisplacementDim>(parameters, constitutive_relation_config);
    }
    else if (type == "Lubby2")
    {
        return MaterialLib::Solids::Lubby2::createLubby2<DisplacementDim>(
            parameters, constitutive_relation_config);
    }
    else
    {
        OGS_FATAL(
            "Cannot construct constitutive relation of given type \'%s\'.",
            type.c_str());
    }
}

template std::unique_ptr<MaterialLib::Solids::MechanicsBase<2>>
createConstitutiveRelation<2>(
    std::vector<std::unique_ptr<ProcessLib::ParameterBase>> const& parameters,
    BaseLib::ConfigTree const& config);

template std::unique_ptr<MaterialLib::Solids::MechanicsBase<3>>
createConstitutiveRelation<3>(
    std::vector<std::unique_ptr<ProcessLib::ParameterBase>> const& parameters,
    BaseLib::ConfigTree const& config);
}
}
