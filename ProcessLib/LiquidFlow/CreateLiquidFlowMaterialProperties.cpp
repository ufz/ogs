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

#include "LiquidFlowMaterialProperties.h"

namespace ProcessLib
{
namespace LiquidFlow
{
class LiquidFlowMaterialProperties;

std::unique_ptr<LiquidFlowMaterialProperties>
createLiquidFlowMaterialProperties(
    BaseLib::ConfigTree const& config, bool const has_material_ids,
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
    //! \ogs_file_param{prj__processes__process__LIQUID_FLOW__material_property__porous_medium}
    auto const& poro_config = config.getConfigSubtree("porous_medium");
    //! \ogs_file_param{prj__processes__process__LIQUID_FLOW__material_property__porous_medium__porous_medium}
    for (auto const& conf : poro_config.getConfigSubtreeList("porous_medium"))
    {
        //! \ogs_file_attr{prj__processes__process__LIQUID_FLOW__material_property__porous_medium__porous_medium__id}
        auto const id = conf.getConfigAttributeOptional<int>("id");
        mat_ids.push_back(*id);

        //! \ogs_file_param{prj__processes__process__LIQUID_FLOW__material_property__porous_medium__porous_medium__permeability}
        auto const& perm_conf = conf.getConfigSubtree("permeability");
        intrinsic_permeability_models.emplace_back(
            MaterialLib::PorousMedium::createPermeabilityModel(perm_conf));

        //! \ogs_file_param{prj__processes__process__LIQUID_FLOW__material_property__porous_medium__porous_medium__porosity}
        auto const& poro_conf = conf.getConfigSubtree("porosity");
        auto n = MaterialLib::PorousMedium::createPorosityModel(poro_conf);
        porosity_models.emplace_back(std::move(n));

        //! \ogs_file_param{prj__processes__process__LIQUID_FLOW__material_property__porous_medium__porous_medium__storage}
        auto const& stora_conf = conf.getConfigSubtree("storage");
        auto beta = MaterialLib::PorousMedium::createStorageModel(stora_conf);
        storage_models.emplace_back(std::move(beta));
    }

    BaseLib::reorderVector(intrinsic_permeability_models, mat_ids);
    BaseLib::reorderVector(porosity_models, mat_ids);
    BaseLib::reorderVector(storage_models, mat_ids);

    return std::unique_ptr<LiquidFlowMaterialProperties>(
        new LiquidFlowMaterialProperties(
            std::move(fluid_properties),
            std::move(intrinsic_permeability_models),
            std::move(porosity_models), std::move(storage_models),
            has_material_ids, material_ids));
}

}  // end of namespace
}  // end of namespace
