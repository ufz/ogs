/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateTwoPhaseFlowMaterialProperties.h"
#include <logog/include/logog.hpp>
#include "BaseLib/reorderVector.h"
#include "MaterialLib/Fluid/FluidProperty.h"
#include "MaterialLib/PorousMedium/Porosity/Porosity.h"
#include "MaterialLib/PorousMedium/Storage/Storage.h"
#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/PropertyVector.h"
#include "ProcessLib/Parameter/Parameter.h"
#include "ProcessLib/Parameter/SpatialPosition.h"
#include "TwoPhaseFlowWithPPMaterialProperties.h"

namespace MaterialLib
{
namespace TwoPhaseFlowWithPP
{
std::unique_ptr<TwoPhaseFlowWithPPMaterialProperties>
CreateTwoPhaseFlowMaterialProperties(
    BaseLib::ConfigTree const& config,
    bool const has_material_ids,
    MeshLib::PropertyVector<int> const& material_ids)
{
    DBUG("Reading material properties of two-phase flow process.");

    //! \ogs_file_param{prj__processes__process__TWOPHASE_FLOW_PP__material_property__fluid}
    auto const& fluid_config = config.getConfigSubtree("fluid");

    // Get fluid properties
    //! \ogs_file_param{prj__processes__process__TWOPHASE_FLOW_PP__material_property__fluid__liquid_density}
    auto const& rho_conf = fluid_config.getConfigSubtree("liquid_density");
    auto _liquid_density =
        MaterialLib::Fluid::createFluidDensityModel(rho_conf);
    //! \ogs_file_param{prj__processes__process__TWOPHASE_FLOW_PP__material_property__fluid__gas_density}
    auto const& rho_gas_conf = fluid_config.getConfigSubtree("gas_density");
    auto _gas_density =
        MaterialLib::Fluid::createFluidDensityModel(rho_gas_conf);
    //! \ogs_file_param{prj__processes__process__TWOPHASE_FLOW_PP__material_property__fluid__liquid_viscosity}
    auto const& mu_conf = fluid_config.getConfigSubtree("liquid_viscosity");
    auto _viscosity = MaterialLib::Fluid::createViscosityModel(mu_conf);
    //! \ogs_file_param{prj__processes__process__TWOPHASE_FLOW_PP__material_property__fluid__gas_viscosity}
    auto const& mu_gas_conf = fluid_config.getConfigSubtree("gas_viscosity");
    auto _gas_viscosity = MaterialLib::Fluid::createViscosityModel(mu_gas_conf);

    // Get porous properties
    std::vector<int> mat_ids;
    std::vector<Eigen::MatrixXd> _intrinsic_permeability_models;
    std::vector<std::unique_ptr<MaterialLib::PorousMedium::Porosity>>
        _porosity_models;
    std::vector<std::unique_ptr<MaterialLib::PorousMedium::Storage>>
        _storage_models;
    //! \ogs_file_param{prj__processes__process__TWOPHASE_FLOW_PP__material_property__porous_medium}
    auto const& poro_config = config.getConfigSubtree("porous_medium");
    //! \ogs_file_param{prj__processes__process__TWOPHASE_FLOW_PP__material_property__porous_medium__porous_medium}
    for (auto const& conf : poro_config.getConfigSubtreeList("porous_medium"))
    {
        //! \ogs_file_attr{prj__processes__process__TWOPHASE_FLOW_PP__material_property__porous_medium__porous_medium__id}
        auto const id = conf.getConfigAttributeOptional<int>("id");
        mat_ids.push_back(*id);

        //! \ogs_file_param{prj__processes__process__TWOPHASE_FLOW_PP__material_property__porous_medium__porous_medium__permeability}
        auto const& permeability_conf = conf.getConfigSubtree("permeability");
        _intrinsic_permeability_models.emplace_back(
            MaterialLib::PorousMedium::createPermeabilityModel(permeability_conf));

        //! \ogs_file_param{prj__processes__process__TWOPHASE_FLOW_PP__material_property__porous_medium__porous_medium__porosity}
        auto const& porosity_conf = conf.getConfigSubtree("porosity");
        auto n = MaterialLib::PorousMedium::createPorosityModel(porosity_conf);
        _porosity_models.emplace_back(std::move(n));

        //! \ogs_file_param{prj__processes__process__TWOPHASE_FLOW_PP__material_property__porous_medium__porous_medium__storage}
        auto const& storage_conf = conf.getConfigSubtree("storage");
        auto beta = MaterialLib::PorousMedium::createStorageModel(storage_conf);
        _storage_models.emplace_back(std::move(beta));
    }

    BaseLib::reorderVector(_intrinsic_permeability_models, mat_ids);
    BaseLib::reorderVector(_porosity_models, mat_ids);
    BaseLib::reorderVector(_storage_models, mat_ids);

    return std::unique_ptr<TwoPhaseFlowWithPPMaterialProperties>{
        new TwoPhaseFlowWithPPMaterialProperties{
            has_material_ids, material_ids, std::move(_liquid_density),
            std::move(_viscosity), std::move(_gas_density),
            std::move(_gas_viscosity), _intrinsic_permeability_models,
            std::move(_porosity_models), std::move(_storage_models)}};
}

}  // end namespace
}  // end namespace
