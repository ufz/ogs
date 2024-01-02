/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "SubmeshResiduumOutputConfig.h"

#include <numeric>
#include <unordered_set>

#include "CreateOutput.h"
#include "CreateOutputConfig.h"

namespace
{
std::vector<std::reference_wrapper<MeshLib::Mesh>>
filterMeshesForResiduumOutput(
    std::vector<std::unique_ptr<MeshLib::Mesh>> const& meshes,
    std::vector<std::string> const& mesh_names_for_output)
{
    std::map<std::string, std::reference_wrapper<MeshLib::Mesh>>
        map_mesh_name_to_mesh;
    for (auto const& mesh : meshes)
    {
        auto const [it, inserted] =
            map_mesh_name_to_mesh.emplace(mesh->getName(), *mesh);

        if (!inserted)
        {
            OGS_FATAL("Duplicate mesh name '{}' detected.", mesh->getName());
        }
    }

    std::vector<std::reference_wrapper<MeshLib::Mesh>> meshes_filtered;

    for (auto const& mesh_name : mesh_names_for_output)
    {
        auto const& mesh = BaseLib::getOrError(
            map_mesh_name_to_mesh, mesh_name,
            "A mesh that has been requested for output is not known to OGS.");

        meshes_filtered.push_back(mesh);
    }

    return meshes_filtered;
}

bool areElementsUnique(std::vector<std::string> const& strings)
{
    std::unordered_set<std::string_view> const strings_set(strings.begin(),
                                                           strings.end());

    return strings_set.size() == strings.size();
}

void checkBulkIDMappingsPresent(MeshLib::Mesh const& mesh)
{
    auto const& properties = mesh.getProperties();

    if (!properties.existsPropertyVector<std::size_t>(
            MeshLib::getBulkIDString(MeshLib::MeshItemType::Node),
            MeshLib::MeshItemType::Node, 1))
    {
        OGS_FATAL(
            "The required nodal property '{}' is missing in mesh '{}' or has "
            "the wrong data type or the wrong number of components",
            MeshLib::getBulkIDString(MeshLib::MeshItemType::Node),
            mesh.getName());
    }

    if (!properties.existsPropertyVector<std::size_t>(
            MeshLib::getBulkIDString(MeshLib::MeshItemType::Cell),
            MeshLib::MeshItemType::Cell, 1))
    {
        OGS_FATAL(
            "The required cell property '{}' is missing in mesh '{}' or has "
            "the wrong data type or the wrong number of components",
            MeshLib::getBulkIDString(MeshLib::MeshItemType::Cell),
            mesh.getName());
    }
}

void checkMatchingElementCounts(
    MeshLib::Mesh const& bulk_mesh,
    std::vector<std::reference_wrapper<MeshLib::Mesh>> const& submesh_refs)
{
    auto const n_elements_bulk = bulk_mesh.getNumberOfElements();
    auto const sum_elements_submeshes = [&submesh_refs]()
    {
        std::size_t n = 0;
        for (auto const& submesh_ref : submesh_refs)
        {
            n += submesh_ref.get().getNumberOfElements();
        }
        return n;
    }();

    if (n_elements_bulk != sum_elements_submeshes)
    {
        OGS_FATAL(
            "The number of bulk mesh elements does not match the sum of all "
            "submesh elements: {} != {}. Hence, the set of all submeshes "
            "cannot be a non-overlapping cover of the bulk mesh.",
            n_elements_bulk, sum_elements_submeshes);
    }
}

std::vector<bool> computeNonOverlappingBulkMeshCoverBySubmeshes(
    MeshLib::Mesh const& bulk_mesh,
    std::vector<std::reference_wrapper<MeshLib::Mesh>> const& submesh_refs)
{
    auto const n_elements_bulk = bulk_mesh.getNumberOfElements();

    std::vector<bool> bulk_element_covered(n_elements_bulk);

    for (auto const& submesh_ref : submesh_refs)
    {
        auto const& submesh = submesh_ref.get();
        auto const& bulk_element_ids =
            *submesh.getProperties().getPropertyVector<std::size_t>(
                "bulk_element_ids", MeshLib::MeshItemType::Cell, 1);

        if (bulk_element_ids.size() != submesh.getNumberOfElements())
        {
            OGS_FATAL(
                "There is something terribly wrong with the mesh '{}'. The "
                "size of 'bulk_element_ids' does not equal the number of "
                "elements in the mesh: {} != {}",
                submesh.getName(), bulk_element_ids.size(),
                submesh.getNumberOfElements());
        }

        for (auto const bulk_element_id : bulk_element_ids)
        {
            // meshes are provided as user input, so we better check the bounds
            // of the contained data
            [[unlikely]] if (bulk_element_id >= n_elements_bulk)
            {
                OGS_FATAL(
                    "Saw bulk element id {} in submesh '{}', but the bulk mesh "
                    "('{}') has only {} elements, i.e., the maximum allowed "
                    "bulk element id is {}.",
                    bulk_element_id, submesh.getName(), bulk_mesh.getName(),
                    n_elements_bulk, n_elements_bulk - 1);
            }

            [[unlikely]] if (bulk_element_covered[bulk_element_id])
            {
                OGS_FATAL(
                    "The bulk element id {} has already been covered by "
                    "another submesh. The second submesh covering this bulk "
                    "element is '{}'.",
                    bulk_element_id, submesh.getName());
            }

            bulk_element_covered[bulk_element_id] = true;
        }
    }

    return bulk_element_covered;
}

void checkNonOverlappingCover(
    MeshLib::Mesh const& bulk_mesh,
    std::vector<std::reference_wrapper<MeshLib::Mesh>> const& submesh_refs)
{
    checkMatchingElementCounts(bulk_mesh, submesh_refs);

    auto const n_elements_bulk = bulk_mesh.getNumberOfElements();

    std::vector<bool> const bulk_element_covered =
        computeNonOverlappingBulkMeshCoverBySubmeshes(bulk_mesh, submesh_refs);

    auto const n_elements_covered =
        std::accumulate(bulk_element_covered.begin(),
                        bulk_element_covered.end(), std::size_t{0});

    if (n_elements_covered == n_elements_bulk)
    {
        return;
    }

    // search first bulk element that has not been covered
    auto const non_covered_it = std::find(bulk_element_covered.begin(),
                                          bulk_element_covered.end(), false);
    auto const non_covered_index =
        std::distance(bulk_element_covered.begin(), non_covered_it);

    OGS_FATAL(
        "The bulk mesh ('{}') is not covered completely by the given "
        "submeshes. Only {} out of {} elements are covered. The first element "
        "that is not covered is #{}.",
        bulk_mesh.getName(), n_elements_covered, n_elements_bulk,
        non_covered_index);
}
}  // namespace

namespace ProcessLib
{
SubmeshResiduumOutputConfig createSubmeshResiduumOutputConfig(
    BaseLib::ConfigTree const& config, std::string const& output_directory,
    std::vector<std::unique_ptr<MeshLib::Mesh>>& meshes)
{
    auto oc = createOutputConfig(config, meshes);

    if (oc.mesh_names_for_output.empty())
    {
        OGS_FATAL(
            "You did not specify any meshes for submesh residuum output.");
    }

    if (!areElementsUnique(oc.mesh_names_for_output))
    {
        OGS_FATAL("The mesh names for submesh residuum output are not unique.");
    }

    auto const meshes_filtered =
        filterMeshesForResiduumOutput(meshes, oc.mesh_names_for_output);

    for (auto const& mesh : meshes_filtered)
    {
        checkBulkIDMappingsPresent(mesh.get());
    }

    auto const& bulk_mesh =
        *meshes.front();  // convention: the first mesh is the bulk mesh
    checkNonOverlappingCover(bulk_mesh, meshes_filtered);

    return {createOutput(std::move(oc), output_directory, meshes),
            std::move(meshes_filtered)};
}
}  // namespace ProcessLib
