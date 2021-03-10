/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "CreateMFront.h"

#include "BaseLib/FileTools.h"
#include "ParameterLib/Utils.h"

#include "MFront.h"
namespace
{
/// Prints info about MFront variables.
void varInfo(std::string const& msg,
             std::vector<mgis::behaviour::Variable> const& vars,
             mgis::behaviour::Hypothesis hypothesis)
{
    INFO("#{:s}: {:d} (array size {:d}).",
         msg,
         vars.size(),
         mgis::behaviour::getArraySize(vars, hypothesis));
    for (auto& var : vars)
    {
        INFO("  --> type `{:s}' with name `{:s}', size {:d}, offset {:d}.",
             MaterialLib::Solids::MFront::varTypeToString(var.type),
             var.name,
             mgis::behaviour::getVariableSize(var, hypothesis),
             mgis::behaviour::getVariableOffset(vars, var.name, hypothesis));
    }
}

/// Prints info about MFront parameters.
void varInfo(std::string const& msg, std::vector<std::string> const& parameters)
{
    INFO("#{:s}: {:d}.", msg, parameters.size());
    // mgis::behaviour::getArraySize(vars, hypothesis));
    for (auto const& parameter : parameters)
    {
        INFO("  --> with name `{:s}'.", parameter);
    }
}
}  // anonymous namespace

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
    INFO("### MFRONT ########################################################");

    //! \ogs_file_param{material__solid__constitutive_relation__type}
    config.checkConfigParameter("type", "MFront");

    auto const library_name =
        //! \ogs_file_param{material__solid__constitutive_relation__MFront__library}
        config.getConfigParameterOptional<std::string>("library");
    auto const lib_path =
        library_name
            ? BaseLib::joinPaths(BaseLib::getProjectDirectory(), *library_name)
            : "libOgsMFrontBehaviour";

    auto const behaviour_name =
        //! \ogs_file_param{material__solid__constitutive_relation__MFront__behaviour}
        config.getConfigParameter<std::string>("behaviour");

    static_assert(DisplacementDim == 2 || DisplacementDim == 3,
                  "Given DisplacementDim not supported.");

    mgis::behaviour::Hypothesis hypothesis;
    if (DisplacementDim == 2)
    {
        // TODO support the axial symmetry modelling hypothesis.
        WARN(
            "The model is defined in 2D. On the material level currently a "
            "plane strain setting is used. In particular it is not checked if "
            "axial symmetry or plane stress are assumed. Special material "
            "behaviour for these settings is currently not supported.");
        hypothesis = mgis::behaviour::Hypothesis::PLANESTRAIN;
    }
    else if (DisplacementDim == 3)
    {
        hypothesis = mgis::behaviour::Hypothesis::TRIDIMENSIONAL;
    }

    auto behaviour =
        mgis::behaviour::load(lib_path, behaviour_name, hypothesis);

    INFO("Behaviour:      `{:s}'.", behaviour.behaviour);
    INFO("Hypothesis:     `{:s}'.", mgis::behaviour::toString(hypothesis));
    INFO("Source:         `{:s}'.", behaviour.source);
    INFO("TFEL version:   `{:s}'.", behaviour.tfel_version);
    INFO("Behaviour type: `{:s}'.", btypeToString(behaviour.btype));
    INFO("Kinematic:      `{:s}'.", toString(behaviour.kinematic));
    INFO("Symmetry:       `{:s}'.", toString(behaviour.symmetry));

    varInfo("Mat. props.", behaviour.mps, hypothesis);
    varInfo("Gradients", behaviour.gradients, hypothesis);
    varInfo("Thdyn. forces", behaviour.thermodynamic_forces, hypothesis);
    varInfo("Int. StVars.", behaviour.isvs, hypothesis);
    varInfo("Ext. StVars.", behaviour.esvs, hypothesis);

    // TODO read parameters from prj file, not yet (2018-11-05) supported by
    // MGIS library.
    varInfo("Real-valued parameters", behaviour.params);
    varInfo("Integer parameters", behaviour.iparams);
    varInfo("Unsigned parameters", behaviour.usparams);

    std::vector<ParameterLib::Parameter<double> const*> material_properties;

    if (!behaviour.mps.empty())
    {
        std::map<std::string, std::string> map_name_to_param;

        // gather material properties from the prj file
        //! \ogs_file_param{material__solid__constitutive_relation__MFront__material_properties}
        auto const mps_config = config.getConfigSubtree("material_properties");
        for (
            auto const mp_config :
            //! \ogs_file_param{material__solid__constitutive_relation__MFront__material_properties__material_property}
            mps_config.getConfigParameterList("material_property"))
        {
            //! \ogs_file_attr{material__solid__constitutive_relation__MFront__material_properties__material_property__name}
            auto name = mp_config.getConfigAttribute<std::string>("name");
            auto const param_name =
                //! \ogs_file_attr{material__solid__constitutive_relation__MFront__material_properties__material_property__parameter}
                mp_config.getConfigAttribute<std::string>("parameter");

            map_name_to_param.emplace(std::move(name), std::move(param_name));
        }

        for (auto& mp : behaviour.mps)
        {
            auto const it = map_name_to_param.find(mp.name);
            if (it == map_name_to_param.end())
                OGS_FATAL(
                    "Material Property `{:s}' has not been configured in the "
                    "project file.",
                    mp.name);

            auto const param_name = it->second;
            auto const num_comp =
                mgis::behaviour::getVariableSize(mp, hypothesis);
            auto const* param = &ParameterLib::findParameter<double>(
                param_name, parameters, num_comp);

            INFO("Using OGS parameter `{:s}' for material property `{:s}'.",
                 param_name, mp.name);

            using V = mgis::behaviour::Variable;
            if (mp.type == V::STENSOR || mp.type == V::TENSOR)
            {
                WARN(
                    "Material property `{:s}' is a tensorial quantity. You, "
                    "the "
                    "user, have to make sure that the component order of "
                    "parameter `{:s}' matches the one required by MFront!",
                    mp.name, param_name);
            }

            material_properties.push_back(param);
            map_name_to_param.erase(it);
        }

        if (!map_name_to_param.empty())
        {
            ERR("Some material parameters that were configured are not used by "
                "the material model.");
            ERR("These parameters are:");

            for (auto& e : map_name_to_param)
            {
                ERR("  name: `{:s}', parameter: `{:s}'.", e.first, e.second);
            }

            OGS_FATAL(
                "Configuration errors occurred. Please fix the project file.");
        }
    }

    INFO("### MFRONT END ####################################################");

    return std::make_unique<MFront<DisplacementDim>>(
        std::move(behaviour), std::move(material_properties),
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
