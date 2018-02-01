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

#include "Permeability/createPermeabilityModel.h"
#include "Porosity/createPorosityModel.h"
#include "Storage/createStorageModel.h"

#include "MeshLib/Mesh.h"

namespace MaterialLib
{
namespace PorousMedium
{
PorousMediaProperties createPorousMediaProperties(
    MeshLib::Mesh& mesh, BaseLib::ConfigTree const& configs,
    std::vector<std::unique_ptr<ProcessLib::ParameterBase>> const& parameters)
{
    DBUG("Create PorousMediaProperties.");

    auto const& porous_medium_configs =
        //! \ogs_file_param{material__porous_medium__porous_medium}
        configs.getConfigSubtree("porous_medium");

    std::vector<std::unique_ptr<MaterialLib::PorousMedium::Permeability>>
        intrinsic_permeability_models;
    std::vector<std::unique_ptr<MaterialLib::PorousMedium::Porosity>>
        porosity_models;
    std::vector<std::unique_ptr<MaterialLib::PorousMedium::Storage>>
        storage_models;

    std::vector<int> mat_ids;
    for (auto const& porous_medium_config :
         //! \ogs_file_param{material__porous_medium__porous_medium}
         porous_medium_configs.getConfigSubtreeList("porous_medium"))
    {
         //! \ogs_file_attr{material__porous_medium__porous_medium__id}
        auto const id = porous_medium_config.getConfigAttribute<int>("id");
        mat_ids.push_back(id);

        auto const& porosity_config =
            //! \ogs_file_param{material__porous_medium__porous_medium__porosity}
            porous_medium_config.getConfigSubtree("porosity");
        porosity_models.emplace_back(
            MaterialLib::PorousMedium::createPorosityModel(porosity_config,
                                                           parameters));

        // Configuration for the intrinsic permeability (only one scalar per
        // element, i.e., the isotropic case is handled at the moment)
        auto const& permeability_config =
            //! \ogs_file_param{material__porous_medium__porous_medium__permeability}
            porous_medium_config.getConfigSubtree("permeability");
        intrinsic_permeability_models.emplace_back(
            MaterialLib::PorousMedium::createPermeabilityModel(
                permeability_config, parameters));

        // Configuration for the specific storage.
        auto const& storage_config =
            //! \ogs_file_param{material__porous_medium__porous_medium__storage}
            porous_medium_config.getConfigSubtree("storage");
        storage_models.emplace_back(
            MaterialLib::PorousMedium::createStorageModel(storage_config));
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
    int const max_material_id =
        *std::max_element(material_ids.cbegin(), material_ids.cend());

    if (max_material_id > mat_ids.size() - 1)
        OGS_FATAL(
            "The maximum value of MaterialIDs in mesh is %d. As the "
            "given number of porous media definitions in the project "
            "file is %d, the maximum value of MaterialIDs in mesh must be %d "
            "(index starts with zero).",
            max_material_id, mat_ids.size(), max_material_id - 1);

    if (max_material_id < mat_ids.size() - 1)
        WARN(
            "There are %d porous medium definitions in the project file but "
            "only %d different values in the MaterialIDs vector/data_array in "
            "the mesh.",
            mat_ids.size(), max_material_id - 1);

    if (mat_ids.back() != static_cast<int>(mat_ids.size()) - 1)
        OGS_FATAL(
            "The ids in the porous media definitions in the project file have "
            "to be sequential, starting with the id zero.");

    return PorousMediaProperties{std::move(porosity_models),
                                 std::move(intrinsic_permeability_models),
                                 std::move(storage_models),
                                 std::move(material_ids)};
}

}  // namespace ComponentTransport
}  // namespace ProcessLib
