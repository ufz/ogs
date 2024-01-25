/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on January 25, 2024, 10:49 AM
 */

#include "CreateCoordinateSystem.h"

#include <string_view>

#include "BaseLib/ConfigTree.h"
#include "Parameter.h"
#include "Utils.h"

namespace ParameterLib
{
struct ParameterBase;
struct CoordinateSystem;

// Note: the function is used for parsing
// 1. base1 for the case of an implicit base0.
// 2. or base2.
Parameter<double> const* parseBase1OrBase2(
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    int const expected_component_number,
    std::string_view const base_tag_name)
{
    auto const base_parameter_name = config.getValue<std::string>();
    auto const& basis_vector = ParameterLib::findParameter<double>(
        base_parameter_name, parameters, 0 /* any dimension */);

    int const component_number = basis_vector.getNumberOfGlobalComponents();
    if (base_tag_name == "basis_vector_1" && component_number != 2)
    {
        OGS_FATAL(
            "The case of implicit \"basis_vector_0\" and "
            "explicit \"basis_vector_1\" is for a 2D coordinate system. The "
            "parameter for \"basis_vector_1\", {:s}, must have "
            "two components but it has {:d}. In addition, \"basis_vector_2\" "
            "should not exist in this case.",
            base_parameter_name, component_number);
    }

    if (component_number != expected_component_number)
    {
        OGS_FATAL(
            "The read parameter `{:s}' for tag {:s} has the wrong number of "
            "components ({:d} instead of {:d}).",
            base_parameter_name, base_tag_name, component_number,
            expected_component_number);
    }

    return &basis_vector;
}

void checkThirdBaseExistanceFor2D(BaseLib::ConfigTree const& config)
{
    auto const base2_config =
        //! \ogs_file_param{prj__local_coordinate_system__basis_vector_2}
        config.getConfigSubtreeOptional("basis_vector_2");
    if (base2_config)
    {
        OGS_FATAL(
            "The tag \"basis_vector_2\" is not needed for a 2D local "
            "coordinate system.");
    }
}

void confirmThirdBaseExplicit(BaseLib::ConfigTree const& config)
{
    auto const implicit_base2 =
        //! \ogs_file_attr{prj__local_coordinate_system__basis_vector_2__implicit}
        config.getConfigAttribute<bool>("implicit", false);
    if (implicit_base2)
    {
        OGS_FATAL("basis_vector_2 must be explicit.");
    }
}

std::optional<ParameterLib::CoordinateSystem>
createCoordinateSystemWithImplicitBase(
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters)
{
    //! \ogs_file_param{prj__local_coordinate_system__basis_vector_1}
    auto const& config_base1 = config.getConfigSubtree("basis_vector_1");

    // Parse the second vectors.
    auto const implicit_base1 =
        //! \ogs_file_attr{prj__local_coordinate_system__basis_vector_1__implicit}
        config_base1.getConfigAttribute<bool>("implicit", false);

    if (!implicit_base1)
    {
        // For 2D problem.
        // Only the following case is accepted:
        // <basis_vector_0 implicit = false/>
        // <basis_vector_1> [two-component parameter] </basis_vector_1>.
        int const expected_component_number = 2;
        auto const basis_vector_1 =
            parseBase1OrBase2(config_base1, parameters,
                              expected_component_number, "basis_vector_1");

        checkThirdBaseExistanceFor2D(config);

        return ParameterLib::CoordinateSystem(*basis_vector_1);
    }

    // For 3D problem.
    // Only the following case is accepted:
    // <basis_vector_0 implicit = false/>
    // <basis_vector_1 implicit = false>.
    // <basis_vector_2> [three-component parameter] </basis_vector_2>.

    // Parse the third basis vector, e2, for the three dimensional system.
    auto const config_base2 =
        //! \ogs_file_param{prj__local_coordinate_system__basis_vector_2}
        config.getConfigSubtreeOptional("basis_vector_2");
    if (!config_base2)
    {
        OGS_FATAL(
            "Both \"basis_vector_0\" and \"basis_vector_1\" are implicit but "
            "\"basis_vector_2\" does not exist. If 2D coordinate system is "
            "considered, please change \"basis_vector_1\" to explicit. If 3D "
            "coordinate system is considered, please add \"basis_vector_2\".");
    }

    confirmThirdBaseExplicit(*config_base2);

    int const expected_component_number = 3;
    auto const basis_vector_2 = parseBase1OrBase2(
        *config_base2, parameters, expected_component_number, "basis_vector_2");
    return ParameterLib::CoordinateSystem(*basis_vector_2);
}

std::optional<ParameterLib::CoordinateSystem> createCoordinateSystem(
    std::optional<BaseLib::ConfigTree> const& config,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters)
{
    if (!config)
    {
        return {};
    }

    //
    // Fetch the first basis vector; its length defines the dimension.
    //
    //! \ogs_file_param{prj__local_coordinate_system__basis_vector_0}
    auto const& config_base0 = config->getConfigSubtree("basis_vector_0");
    auto const implicit_base0 =
        //! \ogs_file_attr{prj__local_coordinate_system__basis_vector_0__implicit}
        config_base0.getConfigAttribute<bool>("implicit", false);

    // base0 and base1 can be implicit. If base0 is implicit, check whether
    // base1 is implicit.
    // If base1 is explicit, create a 2D system if its components is 2 or quit
    // if its components is not equal to 2.
    // Otherwise, read base2, create a 3D system if its components is 3 or quit
    // if its components is not equal to 3.
    if (implicit_base0)
    {
        return createCoordinateSystemWithImplicitBase(*config, parameters);
    }

    auto const base0_name = config_base0.getValue<std::string>();
    auto const& basis_vector_0 = ParameterLib::findParameter<double>(
        base0_name, parameters, 0 /* any dimension */);
    int const dimension = basis_vector_0.getNumberOfGlobalComponents();
    // check dimension
    if (dimension != 2 && dimension != 3)
    {
        OGS_FATAL(
            "Basis vector parameter '{:s}' must have two or three components, "
            "but it has {:d}.",
            base0_name, dimension);
    }

    // Fetch the second basis vector.
    //! \ogs_file_param{prj__local_coordinate_system__basis_vector_1}
    auto const& config_base1 = config->getConfigSubtree("basis_vector_1");
    auto const implicit_base1 =
        //! \ogs_file_attr{prj__local_coordinate_system__basis_vector_1__implicit}
        config_base1.getConfigAttribute<bool>("implicit", false);
    if (implicit_base1)
    {
        OGS_FATAL(
            "Since basis_vector_0 is explicitly defined, basis_vector_1"
            " must be explicit as well.");
    }
    auto const base1_name = config_base1.getValue<std::string>();
    auto const& basis_vector_1 =
        ParameterLib::findParameter<double>(base1_name, parameters, dimension);

    if (dimension == 2)
    {
        checkThirdBaseExistanceFor2D(*config);
        return ParameterLib::CoordinateSystem{basis_vector_0, basis_vector_1};
    }

    // Parse the third basis vector, e2, for the three dimensional system.
    //! \ogs_file_param{prj__local_coordinate_system__basis_vector_2}
    auto const& config_base2 = config->getConfigSubtree("basis_vector_2");
    confirmThirdBaseExplicit(config_base2);

    int const expected_component_number = 3;
    auto const basis_vector_2 = parseBase1OrBase2(
        config_base2, parameters, expected_component_number, "basis_vector_2");

    return ParameterLib::CoordinateSystem{basis_vector_0, basis_vector_1,
                                          *basis_vector_2};
}
}  // namespace ParameterLib
