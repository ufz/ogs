/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on November 7, 2023, 2:57 PM
 */

#pragma once
#include <memory>

#include "BaseLib/ConfigTree.h"
#include "ParameterLib/ConstantParameter.h"
#include "ParameterLib/CoordinateSystem.h"
#include "Tests/TestTools.h"

namespace Tests
{
template <typename ConstitutiveLaw>
std::unique_ptr<ConstitutiveLaw> createTestConstitutiveRelation(
    const char xml[],
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    std::optional<ParameterLib::CoordinateSystem> const& coordinate_system,
    const bool skip_type_checking,
    std::function<std::unique_ptr<ConstitutiveLaw>(
        std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
            parameters,
        std::optional<ParameterLib::CoordinateSystem> const& coordinate_system,
        BaseLib::ConfigTree const& config, const bool skip_type_checking)>
        createConstitutiveRelation)
{
    auto ptree = Tests::readXml(xml);
    BaseLib::ConfigTree conf(std::move(ptree), "", BaseLib::ConfigTree::onerror,
                             BaseLib::ConfigTree::onwarning);

    auto const& sub_config = conf.getConfigSubtree("constitutive_relation");

    auto const type = sub_config.peekConfigParameter<std::string>("type");

    return createConstitutiveRelation(parameters, coordinate_system, sub_config,
                                      skip_type_checking);
}

template <typename ConstitutiveLaw>
std::unique_ptr<ConstitutiveLaw> createTestConstitutiveRelation(
    const char xml[],
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    const bool skip_type_checking,
    std::function<std::unique_ptr<ConstitutiveLaw>(
        std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
            parameters,
        BaseLib::ConfigTree const& config, const bool skip_type_checking)>
        createConstitutiveRelation)
{
    auto ptree = Tests::readXml(xml);
    BaseLib::ConfigTree conf(std::move(ptree), "", BaseLib::ConfigTree::onerror,
                             BaseLib::ConfigTree::onwarning);

    auto const& sub_config = conf.getConfigSubtree("constitutive_relation");

    auto const type = sub_config.peekConfigParameter<std::string>("type");

    return createConstitutiveRelation(parameters, sub_config,
                                      skip_type_checking);
}

}  // namespace Tests
