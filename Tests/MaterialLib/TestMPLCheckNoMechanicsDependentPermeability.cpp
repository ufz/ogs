// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <map>
#include <memory>
#include <stdexcept>

#include "MaterialLib/MPL/MaterialSpatialDistributionMap.h"
#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Phase.h"
#include "MaterialLib/MPL/Properties/Constant.h"
#include "MaterialLib/MPL/Properties/EmbeddedFracturePermeability.h"
#include "MaterialLib/MPL/Properties/OrthotropicEmbeddedFracturePermeability.h"
#include "MaterialLib/MPL/Properties/PermeabilityMohrCoulombFailureIndexModel.h"
#include "MaterialLib/MPL/Properties/StrainDependentPermeability.h"
#include "MaterialLib/MPL/Utils/CheckNoMechanicsDependentPermeability.h"
#include "ParameterLib/ConstantParameter.h"

namespace MPL = MaterialPropertyLib;

namespace
{
std::map<int, std::shared_ptr<MPL::Medium>> createMedia(
    std::unique_ptr<MPL::Property>&& permeability)
{
    auto properties = std::make_unique<MPL::PropertyArray>();
    (*properties)[MPL::PropertyType::permeability] = std::move(permeability);

    std::map<int, std::shared_ptr<MPL::Medium>> media;
    media[0] = std::make_shared<MPL::Medium>(
        0, std::vector<std::unique_ptr<MPL::Phase>>{}, std::move(properties));
    return media;
}

void checkMedia(std::map<int, std::shared_ptr<MPL::Medium>> const& media)
{
    MPL::MaterialSpatialDistributionMap const media_map(media, nullptr);
    MPL::checkNoMechanicsDependentPermeability(media_map, "TEST_PROCESS");
}
}  // namespace

TEST(MaterialPropertyLib, CheckNoMechanicsDependentPermeabilityAcceptsConstant)
{
    auto const media =
        createMedia(std::make_unique<MPL::Constant>("permeability", 1e-12));

    EXPECT_NO_THROW(checkMedia(media));
}

TEST(MaterialPropertyLib,
     CheckNoMechanicsDependentPermeabilityRefusesStrainDependent)
{
    ParameterLib::ConstantParameter<double> const k0("k0", 1e-20);
    auto const media =
        createMedia(std::make_unique<MPL::StrainDependentPermeability<2>>(
            "permeability", k0, 4000.0, 300.0, 50.0, 1e-22, 1e-10, nullptr));

    EXPECT_THROW(checkMedia(media), std::runtime_error);
}

TEST(MaterialPropertyLib,
     CheckNoMechanicsDependentPermeabilityRefusesStressDependent)
{
    ParameterLib::ConstantParameter<double> const k0("k0", 1e-20);
    auto const media = createMedia(
        std::make_unique<MPL::PermeabilityMohrCoulombFailureIndexModel<3>>(
            "permeability", k0, 1e-19, 3.0, 1e6, 15.0, 1e-10, 1e6, nullptr));

    EXPECT_THROW(checkMedia(media), std::runtime_error);
}

TEST(MaterialPropertyLib,
     CheckNoMechanicsDependentPermeabilityRefusesEmbeddedFracture)
{
    ParameterLib::ConstantParameter<double> const angle("angle", 0.0);
    auto const media =
        createMedia(std::make_unique<MPL::EmbeddedFracturePermeability<3>>(
            "permeability", Eigen::Vector3d{1.0, 0.0, 0.0}, true, 1e-20, 1e-5,
            1e-1, 1e-3, angle, angle, 1.0));

    EXPECT_THROW(checkMedia(media), std::runtime_error);
}

TEST(MaterialPropertyLib,
     CheckNoMechanicsDependentPermeabilityRefusesOrthotropicEmbeddedFracture)
{
    ParameterLib::ConstantParameter<double> const angle("angle", 0.0);
    ParameterLib::ConstantParameter<double> const k0("k0", 1e-20);
    auto const media = createMedia(
        std::make_unique<MPL::OrthotropicEmbeddedFracturePermeability<3>>(
            "permeability", std::vector<double>{1e-1, 1e-1, 1e-1},
            std::vector<double>{1e-3, 1e-3, 1e-3}, Eigen::Matrix3d::Identity(),
            k0, angle, angle, 1.0));

    EXPECT_THROW(checkMedia(media), std::runtime_error);
}

TEST(MaterialPropertyLib,
     CheckNoMechanicsDependentPermeabilitySkipsMediumWithoutPermeability)
{
    std::map<int, std::shared_ptr<MPL::Medium>> media;
    media[0] = std::make_shared<MPL::Medium>(
        0, std::vector<std::unique_ptr<MPL::Phase>>{},
        std::make_unique<MPL::PropertyArray>());

    EXPECT_NO_THROW(checkMedia(media));
}
