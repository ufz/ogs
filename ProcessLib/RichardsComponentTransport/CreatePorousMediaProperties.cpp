/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "CreatePorousMediaProperties.h"

#include "BaseLib/Algorithm.h"
#include "MaterialLib/PorousMedium/UnsaturatedProperty/CapillaryPressure/CreateCapillaryPressureModel.h"
#include "MeshLib/Mesh.h"

namespace ProcessLib
{
namespace RichardsComponentTransport
{
PorousMediaProperties createPorousMediaProperties(
    MeshLib::Mesh& mesh, BaseLib::ConfigTree const& porous_medium_configs)
{
    DBUG("Create PorousMediaProperties.");

    std::vector<
        std::unique_ptr<MaterialLib::PorousMedium::CapillaryPressureSaturation>>
        capillary_pressure_models;

    std::vector<int> mat_ids;
    for (
        auto const& porous_medium_config :
        //! \ogs_file_param{prj__processes__process__RichardsComponentTransport__porous_medium__porous_medium}
        porous_medium_configs.getConfigSubtreeList("porous_medium"))
    {
        //! \ogs_file_attr{prj__processes__process__RichardsComponentTransport__porous_medium__porous_medium__id}
        auto const id = porous_medium_config.getConfigAttribute<int>("id");
        mat_ids.push_back(id);

        auto const& capillary_pressure_config =
            //! \ogs_file_param{prj__processes__process__RichardsComponentTransport__porous_medium__porous_medium__capillary_pressure}
            porous_medium_config.getConfigSubtree("capillary_pressure");
        auto capillary_pressure =
            MaterialLib::PorousMedium::createCapillaryPressureModel(
                capillary_pressure_config);
        capillary_pressure_models.emplace_back(std::move(capillary_pressure));
    }

    std::vector<int> material_ids(mesh.getNumberOfElements());
    if (mesh.getProperties().existsPropertyVector<int>("MaterialIDs"))
    {
        auto const& mesh_material_ids =
            mesh.getProperties().getPropertyVector<int>("MaterialIDs");
        material_ids.reserve(mesh_material_ids->size());
        std::copy(mesh_material_ids->cbegin(), mesh_material_ids->cend(),
                  material_ids.begin());
    }

    return PorousMediaProperties{std::move(capillary_pressure_models),
                                 std::move(material_ids)};
}

}  // namespace RichardsComponentTransport
}  // namespace ProcessLib
