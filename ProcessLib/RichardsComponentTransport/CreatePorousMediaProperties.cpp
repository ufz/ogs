/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreatePorousMediaProperties.h"

#include "BaseLib/reorderVector.h"

#include "MaterialLib/PorousMedium/Permeability/createPermeabilityModel.h"
#include "MaterialLib/PorousMedium/Porosity/createPorosityModel.h"
#include "MaterialLib/PorousMedium/Storage/createStorageModel.h"
#include "MaterialLib/PorousMedium/UnsaturatedProperty/CapillaryPressure/CreateCapillaryPressureModel.h"
#include "MaterialLib/PorousMedium/UnsaturatedProperty/RelativePermeability/CreateRelativePermeabilityModel.h"

#include "MeshLib/Mesh.h"

namespace ProcessLib
{
namespace RichardsComponentTransport
{
PorousMediaProperties createPorousMediaProperties(
    MeshLib::Mesh& mesh, BaseLib::ConfigTree const& porous_medium_configs,
    std::vector<std::unique_ptr<ProcessLib::ParameterBase>> const& parameters)
{
    DBUG("Create PorousMediaProperties.");

    std::vector<std::unique_ptr<MaterialLib::PorousMedium::Permeability>>
        intrinsic_permeability_models;
    std::vector<std::unique_ptr<MaterialLib::PorousMedium::Porosity>>
        porosity_models;
    std::vector<std::unique_ptr<MaterialLib::PorousMedium::Storage>>
        storage_models;
    std::vector<
        std::unique_ptr<MaterialLib::PorousMedium::CapillaryPressureSaturation>>
        capillary_pressure_models;
    std::vector<
        std::unique_ptr<MaterialLib::PorousMedium::RelativePermeability>>
        relative_permeability_models;

    std::vector<int> mat_ids;
    for (auto const& porous_medium_config :
         //! \ogs_file_param{prj__processes__process__RichardsComponentTransport__porous_medium__porous_medium}
         porous_medium_configs.getConfigSubtreeList("porous_medium"))
    {
         //! \ogs_file_attr{prj__processes__process__RichardsComponentTransport__porous_medium__porous_medium__id}
        auto const id = porous_medium_config.getConfigAttribute<int>("id");
        mat_ids.push_back(id);

        auto const& porosity_config =
            //! \ogs_file_param{prj__processes__process__RichardsComponentTransport__porous_medium__porous_medium__porosity}
            porous_medium_config.getConfigSubtree("porosity");
        porosity_models.emplace_back(
            MaterialLib::PorousMedium::createPorosityModel(porosity_config,
                                                           parameters));

        // Configuration for the intrinsic permeability
        auto const& permeability_config =
            //! \ogs_file_param{prj__processes__process__RichardsComponentTransport__porous_medium__porous_medium__permeability}
            porous_medium_config.getConfigSubtree("permeability");
        intrinsic_permeability_models.emplace_back(
            MaterialLib::PorousMedium::createPermeabilityModel(
                permeability_config, parameters));

        // Configuration for the specific storage.
        auto const& storage_config =
            //! \ogs_file_param{prj__processes__process__RichardsComponentTransport__porous_medium__porous_medium__storage}
            porous_medium_config.getConfigSubtree("storage");
        storage_models.emplace_back(
            MaterialLib::PorousMedium::createStorageModel(storage_config));

        auto const& capillary_pressure_config =
            //! \ogs_file_param{prj__processes__process__RichardsComponentTransport__porous_medium__porous_medium__capillary_pressure}
            porous_medium_config.getConfigSubtree("capillary_pressure");
        auto capillary_pressure = MaterialLib::PorousMedium::createCapillaryPressureModel(
            capillary_pressure_config);
        capillary_pressure_models.emplace_back(std::move(capillary_pressure));

        auto const& krel_config =
            //! \ogs_file_param{prj__processes__process__RichardsComponentTransport__porous_medium__porous_medium__relative_permeability}
            porous_medium_config.getConfigSubtree("relative_permeability");
        auto krel = MaterialLib::PorousMedium::createRelativePermeabilityModel(
            krel_config);
        relative_permeability_models.emplace_back(std::move(krel));
    }

    BaseLib::reorderVector(intrinsic_permeability_models, mat_ids);
    BaseLib::reorderVector(porosity_models, mat_ids);
    BaseLib::reorderVector(storage_models, mat_ids);

    std::vector<int> material_ids(mesh.getNumberOfElements());
    if (mesh.getProperties().existsPropertyVector<int>("MaterialIDs"))
    {
        auto const& mesh_material_ids =
            mesh.getProperties().getPropertyVector<int>("MaterialIDs");
        material_ids.reserve(mesh_material_ids->size());
        std::copy(mesh_material_ids->cbegin(), mesh_material_ids->cend(),
                  material_ids.begin());
    }

    return PorousMediaProperties{std::move(porosity_models),
                                 std::move(intrinsic_permeability_models),
                                 std::move(storage_models),
                                 std::move(capillary_pressure_models),
                                 std::move(relative_permeability_models),
                                 std::move(material_ids)};
}

}  // namespace ComponentTransport
}  // namespace ProcessLib
