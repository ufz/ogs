/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateTwoPhaseFlowPrhoMaterialProperties.h"
#include "BaseLib/Algorithm.h"
#include "BaseLib/Logging.h"
#include "MaterialLib/Fluid/FluidProperty.h"
#include "MaterialLib/PorousMedium/Porosity/Porosity.h"
#include "MaterialLib/PorousMedium/Storage/Storage.h"
#include "MaterialLib/PorousMedium/UnsaturatedProperty/CapillaryPressure/CapillaryPressureSaturation.h"
#include "MaterialLib/PorousMedium/UnsaturatedProperty/CapillaryPressure/CreateCapillaryPressureModel.h"
#include "MaterialLib/PorousMedium/UnsaturatedProperty/RelativePermeability/CreateRelativePermeabilityModel.h"
#include "MaterialLib/PorousMedium/UnsaturatedProperty/RelativePermeability/RelativePermeability.h"
#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/PropertyVector.h"
#include "ParameterLib/Parameter.h"
#include "ParameterLib/SpatialPosition.h"
#include "TwoPhaseFlowWithPrhoMaterialProperties.h"

namespace ProcessLib
{
namespace TwoPhaseFlowWithPrho
{
std::unique_ptr<TwoPhaseFlowWithPrhoMaterialProperties>
createTwoPhaseFlowPrhoMaterialProperties(
    BaseLib::ConfigTree const& config,
    boost::optional<MeshLib::PropertyVector<int> const&>
        material_ids,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters)
{
    DBUG("Reading material properties of two-phase flow process.");

    //! \ogs_file_param{prj__processes__process__TWOPHASE_FLOW_PRHO__material_property__fluid}
    auto const& fluid_config = config.getConfigSubtree("fluid");

    // Get fluid properties
    //! \ogs_file_param{prj__processes__process__TWOPHASE_FLOW_PRHO__material_property__liquid_density}
    auto const& rho_conf = fluid_config.getConfigSubtree("liquid_density");
    auto liquid_density_ =
        MaterialLib::Fluid::createFluidDensityModel(rho_conf);
    //! \ogs_file_param{prj__processes__process__TWOPHASE_FLOW_PRHO__material_property__gas_density}
    auto const& rho_gas_conf = fluid_config.getConfigSubtree("gas_density");
    auto gas_density_ =
        MaterialLib::Fluid::createFluidDensityModel(rho_gas_conf);
    //! \ogs_file_param{prj__processes__process__TWOPHASE_FLOW_PRHO__material_property__liquid_viscosity}
    auto const& mu_conf = fluid_config.getConfigSubtree("liquid_viscosity");
    auto viscosity_ = MaterialLib::Fluid::createViscosityModel(mu_conf);
    //! \ogs_file_param{prj__processes__process__TWOPHASE_FLOW_PRHO__material_property__gas_viscosity}
    auto const& mu_gas_conf = fluid_config.getConfigSubtree("gas_viscosity");
    auto gas_viscosity_ = MaterialLib::Fluid::createViscosityModel(mu_gas_conf);

    // Get porous properties
    std::vector<int> mat_ids;
    std::vector<int> mat_krel_ids;
    std::vector<std::unique_ptr<MaterialLib::PorousMedium::Permeability>>
        intrinsic_permeability_models_;
    std::vector<std::unique_ptr<MaterialLib::PorousMedium::Porosity>>
        porosity_models_;
    std::vector<std::unique_ptr<MaterialLib::PorousMedium::Storage>>
        storage_models_;
    std::vector<
        std::unique_ptr<MaterialLib::PorousMedium::CapillaryPressureSaturation>>
        capillary_pressure_models_;
    std::vector<
        std::unique_ptr<MaterialLib::PorousMedium::RelativePermeability>>
        relative_permeability_models_;

    //! \ogs_file_param{prj__processes__process__TWOPHASE_FLOW_PRHO__material_property__porous_medium}
    auto const& poro_config = config.getConfigSubtree("porous_medium");
    //! \ogs_file_param{prj__processes__process__TWOPHASE_FLOW_PRHO__material_property__porous_medium__porous_medium}
    for (auto const& conf : poro_config.getConfigSubtreeList("porous_medium"))
    {
        //! \ogs_file_attr{prj__processes__process__TWOPHASE_FLOW_PRHO__material_property__porous_medium__porous_medium__id}
        auto const id = conf.getConfigAttributeOptional<int>("id");
        mat_ids.push_back(*id);

        //! \ogs_file_param{prj__processes__process__TWOPHASE_FLOW_PRHO__material_property__porous_medium__porous_medium__permeability}
        auto const& permeability_conf = conf.getConfigSubtree("permeability");
        intrinsic_permeability_models_.emplace_back(
            MaterialLib::PorousMedium::createPermeabilityModel(
                permeability_conf, parameters));

        //! \ogs_file_param{prj__processes__process__TWOPHASE_FLOW_PRHO__material_property__porous_medium__porous_medium__porosity}
        auto const& porosity_conf = conf.getConfigSubtree("porosity");
        auto n = MaterialLib::PorousMedium::createPorosityModel(porosity_conf,
                                                                parameters);
        porosity_models_.emplace_back(std::move(n));

        //! \ogs_file_param{prj__processes__process__TWOPHASE_FLOW_PRHO__material_property__porous_medium__porous_medium__storage}
        auto const& storage_conf = conf.getConfigSubtree("storage");
        auto beta = MaterialLib::PorousMedium::createStorageModel(storage_conf);
        storage_models_.emplace_back(std::move(beta));

        auto const& capillary_pressure_conf =
            //! \ogs_file_param{prj__processes__process__TWOPHASE_FLOW_PRHO__material_property__porous_medium__porous_medium__capillary_pressure}
            conf.getConfigSubtree("capillary_pressure");
        auto pc = MaterialLib::PorousMedium::createCapillaryPressureModel(
            capillary_pressure_conf);
        capillary_pressure_models_.emplace_back(std::move(pc));

        auto const& krel_config =
            //! \ogs_file_param{prj__processes__process__TWOPHASE_FLOW_PRHO__material_property__porous_medium__porous_medium__relative_permeability}
            conf.getConfigSubtree("relative_permeability");
        for (
            auto const& krel_conf :
            //! \ogs_file_param{prj__processes__process__TWOPHASE_FLOW_PRHO__material_property__porous_medium__porous_medium__relative_permeability__relative_permeability}
            krel_config.getConfigSubtreeList("relative_permeability"))
        {
            auto const krel_id =
                //! \ogs_file_attr{prj__processes__process__TWOPHASE_FLOW_PRHO__material_property__porous_medium__porous_medium__relative_permeability__relative_permeability__id}
                krel_conf.getConfigAttributeOptional<int>("id");
            mat_krel_ids.push_back(*krel_id);
            auto krel_n =
                MaterialLib::PorousMedium::createRelativePermeabilityModel(
                    krel_conf);
            relative_permeability_models_.emplace_back(std::move(krel_n));
        }
        BaseLib::reorderVector(relative_permeability_models_, mat_krel_ids);
    }

    BaseLib::reorderVector(intrinsic_permeability_models_, mat_ids);
    BaseLib::reorderVector(porosity_models_, mat_ids);
    BaseLib::reorderVector(storage_models_, mat_ids);

    return std::make_unique<TwoPhaseFlowWithPrhoMaterialProperties>(
        material_ids, std::move(liquid_density_), std::move(viscosity_),
        std::move(gas_density_), std::move(gas_viscosity_),
        std::move(intrinsic_permeability_models_), std::move(porosity_models_),
        std::move(storage_models_), std::move(capillary_pressure_models_),
        std::move(relative_permeability_models_));
}

}  // namespace TwoPhaseFlowWithPrho
}  // namespace ProcessLib
