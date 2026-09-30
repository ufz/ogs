// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#ifdef OGS_USE_MFRONT

#include <gmock/gmock-matchers.h>
#include <gtest/gtest.h>

#include "BaseLib/ConfigTree.h"
#include "MaterialLib/SolidModels/MFront/CreateMFrontGeneric.h"
#include "MaterialLib/SolidModels/MFront/Variable.h"
#include "NumLib/Exceptions.h"
#include "Tests/TestTools.h"

namespace MSM = MaterialLib::Solids::MFront;

TEST(MaterialLib_MFrontGeneric, IntegrationFailureIncludesMFrontDiagnostic)
{
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> parameters;
    auto local_coordinate_system = std::nullopt;

    auto ptree = Tests::readXml(R"XML(
        <type>MFront</type>
        <behaviour>CheckIntegrationFailure</behaviour>
        <library path_is_relative_to_prj_file="false">libOgsMFrontBehaviourForUnitTests</library>
        <material_properties />
        )XML");
    BaseLib::ConfigTree config_tree(std::move(ptree), "FILENAME",
                                    &BaseLib::ConfigTree::onerror,
                                    &BaseLib::ConfigTree::onwarning);

    auto mfront_model =
        MSM::createMFrontGeneric<3, boost::mp11::mp_list<MSM::Strain>,
                                 boost::mp11::mp_list<MSM::Stress>,
                                 boost::mp11::mp_list<MSM::Temperature>>(
            parameters, local_coordinate_system, config_tree);

    namespace MPL = MaterialPropertyLib;
    using KV = MathLib::KelvinVector::KelvinVectorType<3>;

    MPL::VariableArray variable_array_prev;
    variable_array_prev.mechanical_strain = KV::Zero().eval();
    variable_array_prev.stress = KV::Zero().eval();
    variable_array_prev.temperature = 293.15;

    MPL::VariableArray variable_array;
    variable_array.mechanical_strain = KV::Zero().eval();
    variable_array.temperature = 293.15;

    auto state = mfront_model->createMaterialStateVariables();
    ParameterLib::SpatialPosition pos;

    try
    {
        static_cast<void>(mfront_model->integrateStress(
            variable_array_prev, variable_array, 1.0, pos, 1.0, *state));
        FAIL() << "MFront integration unexpectedly succeeded.";
    }
    catch (NumLib::AssemblyException const& e)
    {
        EXPECT_THAT(
            e.what(),
            testing::HasSubstr("MFront: integration failed with status -1."));
        EXPECT_THAT(e.what(),
                    testing::HasSubstr("MFront unit-test diagnostic"));
    }
}

#endif
