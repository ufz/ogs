/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   CreateLiquidFlowMaterialProperties.cpp
 *
 *   Created on December 14, 2016, 1:20 PM
 */

#include "CreateLiquidFlowMaterialProperties.h"

#include "BaseLib/ConfigTree.h"
#include "BaseLib/reorderVector.h"

#include "MeshLib/PropertyVector.h"

#include "MaterialLib/Fluid/FluidProperty.h"
#include "MaterialLib/PorousMedium/Porosity/Porosity.h"
#include "MaterialLib/PorousMedium/Storage/Storage.h"
#include "MaterialLib/Fluid/FluidProperties/CreateFluidProperties.h"

#include "MaterialLib/Fluid/FluidPropertyHeaders.h"
#include "MaterialLib/PorousMedium/PorousPropertyHeaders.h"

#include "ProcessLib/Utils/ProcessUtils.h"
#include "ProcessLib/Parameter/ConstantParameter.h"

#include "LiquidFlowMaterialProperties.h"

namespace ProcessLib
{
namespace LiquidFlow
{
class LiquidFlowMaterialProperties;

std::unique_ptr<LiquidFlowMaterialProperties>
createLiquidFlowMaterialProperties(
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    bool const has_material_ids,
    MeshLib::PropertyVector<int> const& material_ids)
{
    DBUG("Reading material properties of liquid flow process.");

    //! \ogs_file_param{prj__processes__process__LIQUID_FLOW__material_property__fluid}
    auto const& fluid_config = config.getConfigSubtree("fluid");
    auto fluid_properties =
        MaterialLib::Fluid::createFluidProperties(fluid_config);

    // Get porous properties
    std::vector<Eigen::MatrixXd> intrinsic_permeability_models;
    std::vector<std::unique_ptr<MaterialLib::PorousMedium::Porosity>>
        porosity_models;
    std::vector<std::unique_ptr<MaterialLib::PorousMedium::Storage>>
        storage_models;

    std::vector<int> mat_ids;
    auto const& porous_medium_configs =
        //! \ogs_file_param{prj__processes__process__LIQUID_FLOW__material_property__porous_medium}
        config.getConfigSubtree("porous_medium");
    for (
        auto const& porous_medium_config :
        //! \ogs_file_param{prj__processes__process__LIQUID_FLOW__material_property__porous_medium__porous_medium}
        porous_medium_configs.getConfigSubtreeList("porous_medium"))
    {
        //! \ogs_file_attr{prj__processes__process__LIQUID_FLOW__material_property__porous_medium__porous_medium__id}
        auto const id = porous_medium_config.getConfigAttribute<int>("id");
        mat_ids.push_back(id);

        auto const& permeability_config =
            //! \ogs_file_param{prj__processes__process__LIQUID_FLOW__material_property__porous_medium__porous_medium__permeability}
            porous_medium_config.getConfigSubtree("permeability");
        intrinsic_permeability_models.emplace_back(
            MaterialLib::PorousMedium::createPermeabilityModel(
                permeability_config));

        auto const& porosity_config =
            //! \ogs_file_param{prj__processes__process__LIQUID_FLOW__material_property__porous_medium__porous_medium__porosity}
            porous_medium_config.getConfigSubtree("porosity");
        auto n =
            MaterialLib::PorousMedium::createPorosityModel(porosity_config);
        porosity_models.emplace_back(std::move(n));

        auto const& storage_config =
            //! \ogs_file_param{prj__processes__process__LIQUID_FLOW__material_property__porous_medium__porous_medium__storage}
            porous_medium_config.getConfigSubtree("storage");
        auto beta =
            MaterialLib::PorousMedium::createStorageModel(storage_config);
        storage_models.emplace_back(std::move(beta));
    }

    BaseLib::reorderVector(intrinsic_permeability_models, mat_ids);
    BaseLib::reorderVector(porosity_models, mat_ids);
    BaseLib::reorderVector(storage_models, mat_ids);

    //! \ogs_file_param{prj__processes__process__LIQUID_FLOW__material_property__solid}
    auto const solid_config = config.getConfigSubtreeOptional("solid");
    if (solid_config)
    {
        auto& solid_thermal_expansion = findParameter<double>(
            //! \ogs_file_param_special{prj__processes__process__LIQUID_FLOW__material_property__solid__thermal_expansion}
            *solid_config, "thermal_expansion", parameters, 1);
        DBUG("Use \'%s\' as solid thermal expansion.",
             solid_thermal_expansion.name.c_str());
        auto& biot_constant = findParameter<double>(
            //! \ogs_file_param_special{prj__processes__process__LIQUID_FLOW__material_property__solid__biot_constant}
            *solid_config, "biot_constant", parameters, 1);
        return std::unique_ptr<LiquidFlowMaterialProperties>(
            new LiquidFlowMaterialProperties(
                std::move(fluid_properties),
                std::move(intrinsic_permeability_models),
                std::move(porosity_models), std::move(storage_models),
                has_material_ids, material_ids, solid_thermal_expansion,
                biot_constant));
    }
    else
    {
        ConstantParameter<double> void_parameter("void_solid_thermal_expansion",
                                                 0.);
        return std::unique_ptr<LiquidFlowMaterialProperties>(
            new LiquidFlowMaterialProperties(
                std::move(fluid_properties),
                std::move(intrinsic_permeability_models),
                std::move(porosity_models), std::move(storage_models),
                has_material_ids, material_ids, void_parameter,
                void_parameter));
    }
}

}  // end of namespace
}  // end of namespace
