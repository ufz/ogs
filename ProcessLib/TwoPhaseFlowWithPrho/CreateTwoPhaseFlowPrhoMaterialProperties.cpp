/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateTwoPhaseFlowPrhoMaterialProperties.h"
#include <logog/include/logog.hpp>
#include "BaseLib/Algorithm.h"
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
#include "ProcessLib/Parameter/Parameter.h"
#include "ProcessLib/Parameter/SpatialPosition.h"
#include "TwoPhaseFlowWithPrhoMaterialProperties.h"

namespace ProcessLib
{
namespace TwoPhaseFlowWithPrho
{
std::unique_ptr<TwoPhaseFlowWithPrhoMaterialProperties>
createTwoPhaseFlowPrhoMaterialProperties(
    BaseLib::ConfigTree const& config,
    boost::optional<MeshLib::PropertyVector<int> const&> material_ids,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters)
{
    DBUG("Reading material properties of two-phase flow process.");

    //! \ogs_file_param{prj__processes__process__TWOPHASE_FLOW_PRHO__material_property__fluid}
    auto const& fluid_config = config.getConfigSubtree("fluid");

    // Get fluid properties
    //! \ogs_file_param{prj__processes__process__TWOPHASE_FLOW_PRHO__material_property__liquid_density}
    auto const& rho_conf = fluid_config.getConfigSubtree("liquid_density");
    auto _liquid_density =
        MaterialLib::Fluid::createFluidDensityModel(rho_conf);
    //! \ogs_file_param{prj__processes__process__TWOPHASE_FLOW_PRHO__material_property__gas_density}
    auto const& rho_gas_conf = fluid_config.getConfigSubtree("gas_density");
    auto _gas_density =
        MaterialLib::Fluid::createFluidDensityModel(rho_gas_conf);
    //! \ogs_file_param{prj__processes__process__TWOPHASE_FLOW_PRHO__material_property__liquid_viscosity}
    auto const& mu_conf = fluid_config.getConfigSubtree("liquid_viscosity");
    auto _viscosity = MaterialLib::Fluid::createViscosityModel(mu_conf);
    //! \ogs_file_param{prj__processes__process__TWOPHASE_FLOW_PRHO__material_property__gas_viscosity}
    auto const& mu_gas_conf = fluid_config.getConfigSubtree("gas_viscosity");
    auto _gas_viscosity = MaterialLib::Fluid::createViscosityModel(mu_gas_conf);

    // Get porous properties
    std::vector<int> mat_ids;
    std::vector<int> mat_krel_ids;
    std::vector<std::unique_ptr<MaterialLib::PorousMedium::Permeability>>
        _intrinsic_permeability_models;
    std::vector<std::unique_ptr<MaterialLib::PorousMedium::Porosity>>
        _porosity_models;
    std::vector<std::unique_ptr<MaterialLib::PorousMedium::Storage>>
        _storage_models;
    std::vector<
        std::unique_ptr<MaterialLib::PorousMedium::CapillaryPressureSaturation>>
        _capillary_pressure_models;
    std::vector<
        std::unique_ptr<MaterialLib::PorousMedium::RelativePermeability>>
        _relative_permeability_models;

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
        _intrinsic_permeability_models.emplace_back(
            MaterialLib::PorousMedium::createPermeabilityModel(
                permeability_conf, parameters));

        //! \ogs_file_param{prj__processes__process__TWOPHASE_FLOW_PRHO__material_property__porous_medium__porous_medium__porosity}
        auto const& porosity_conf = conf.getConfigSubtree("porosity");
        auto n = MaterialLib::PorousMedium::createPorosityModel(porosity_conf,
                                                                parameters);
        _porosity_models.emplace_back(std::move(n));

        //! \ogs_file_param{prj__processes__process__TWOPHASE_FLOW_PRHO__material_property__porous_medium__porous_medium__storage}
        auto const& storage_conf = conf.getConfigSubtree("storage");
        auto beta = MaterialLib::PorousMedium::createStorageModel(storage_conf);
        _storage_models.emplace_back(std::move(beta));

        auto const& capillary_pressure_conf =
            //! \ogs_file_param{prj__processes__process__TWOPHASE_FLOW_PRHO__material_property__porous_medium__porous_medium__capillary_pressure}
            conf.getConfigSubtree("capillary_pressure");
        auto pc = MaterialLib::PorousMedium::createCapillaryPressureModel(
            capillary_pressure_conf);
        _capillary_pressure_models.emplace_back(std::move(pc));

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
            _relative_permeability_models.emplace_back(std::move(krel_n));
        }
        BaseLib::reorderVector(_relative_permeability_models, mat_krel_ids);
    }

    BaseLib::reorderVector(_intrinsic_permeability_models, mat_ids);
    BaseLib::reorderVector(_porosity_models, mat_ids);
    BaseLib::reorderVector(_storage_models, mat_ids);

    return std::make_unique<TwoPhaseFlowWithPrhoMaterialProperties>(
        material_ids, std::move(_liquid_density), std::move(_viscosity),
        std::move(_gas_density), std::move(_gas_viscosity),
        std::move(_intrinsic_permeability_models), std::move(_porosity_models),
        std::move(_storage_models), std::move(_capillary_pressure_models),
        std::move(_relative_permeability_models));
}

}  // end namespace
}  // end namespace
