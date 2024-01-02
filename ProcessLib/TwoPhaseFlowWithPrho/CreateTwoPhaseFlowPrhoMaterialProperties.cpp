/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateTwoPhaseFlowPrhoMaterialProperties.h"

#include "BaseLib/Algorithm.h"
#include "BaseLib/ConfigTree.h"
#include "BaseLib/Logging.h"
#include "MaterialLib/Fluid/FluidProperty.h"
#include "MaterialLib/PorousMedium/UnsaturatedProperty/CapillaryPressure/CapillaryPressureSaturation.h"
#include "MaterialLib/PorousMedium/UnsaturatedProperty/CapillaryPressure/CreateCapillaryPressureModel.h"
#include "MaterialLib/PorousMedium/UnsaturatedProperty/RelativePermeability/CreateRelativePermeabilityModel.h"
#include "MaterialLib/PorousMedium/UnsaturatedProperty/RelativePermeability/RelativePermeability.h"
#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/PropertyVector.h"
#include "TwoPhaseFlowWithPrhoMaterialProperties.h"

namespace ProcessLib
{
namespace TwoPhaseFlowWithPrho
{
std::unique_ptr<TwoPhaseFlowWithPrhoMaterialProperties>
createTwoPhaseFlowPrhoMaterialProperties(
    BaseLib::ConfigTree const& config,
    MeshLib::PropertyVector<int> const* const material_ids)
{
    DBUG("Reading material properties of two-phase flow process.");

    // Get porous properties
    std::vector<int> mat_ids;
    std::vector<int> mat_krel_ids;
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

    return std::make_unique<TwoPhaseFlowWithPrhoMaterialProperties>(
        material_ids,
        std::move(_capillary_pressure_models),
        std::move(_relative_permeability_models));
}

}  // namespace TwoPhaseFlowWithPrho
}  // namespace ProcessLib
