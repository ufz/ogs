// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <vector>

#include "BaseLib/Logging.h"
#include "LocalAssemblerFactoryForDimGreaterEqualN.h"
#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "NumLib/Fem/Integration/IntegrationMethodProvider.h"

namespace ProcessLib
{
/// \tparam MinElementDim the minimum topological dimension of the mesh
/// elements a local assembler is built for; elements of lower dimension are
/// skipped by the underlying factory (see
/// LocalAssemblerFactoryForDimGreaterEqualN). Defaults to 1, i.e., point
/// elements are excluded, which is correct for bulk-mesh assembly; pass 0
/// (e.g. for ProcessLib::SurfaceFlux) to also build local assemblers for
/// point elements, needed for the boundary of a 1D bulk mesh.
template <int GlobalDim,
          template <typename /* shp fct */,
                    int /* global dim*/> class LocalAssemblerImplementation,
          int MinElementDim = 1, typename LocalAssemblerInterface,
          IntegrationMethodProviderOrIntegrationOrder ProviderOrOrder,
          typename... ExtraCtorArgs>
void createLocalAssemblers(
    std::vector<MeshLib::Element*> const& mesh_elements,
    NumLib::LocalToGlobalIndexMap const& dof_table,
    std::vector<std::unique_ptr<LocalAssemblerInterface>>& local_assemblers,
    ProviderOrOrder const& provider_or_order,
    ExtraCtorArgs&&... extra_ctor_args)
{
    static_assert(
        GlobalDim == 1 || GlobalDim == 2 || GlobalDim == 3,
        "Meshes with dimension greater than three are not supported.");

    DBUG("Create local assemblers.");

    auto const& integration_method_provider =
        getIntegrationMethodProvider(provider_or_order);

    using IntegrationMethodProvider =
        std::remove_cvref_t<decltype(integration_method_provider)>;
    using LocAsmFac = LocalAssemblerFactoryForDimGreaterEqualN<
        MinElementDim, LocalAssemblerInterface, LocalAssemblerImplementation,
        IntegrationMethodProvider, GlobalDim, ExtraCtorArgs...>;

    LocAsmFac factory(dof_table, integration_method_provider);
    local_assemblers.resize(mesh_elements.size());

    DBUG("Calling local assembler builder for all mesh elements.");
    GlobalExecutor::transformDereferenced(
        factory, mesh_elements, local_assemblers,
        std::forward<ExtraCtorArgs>(extra_ctor_args)...);
}

/*! Creates local assemblers for each element of the given \c mesh.
 *
 * \tparam LocalAssemblerImplementation the individual local assembler type
 * \tparam MinElementDim see the overload above.
 * \tparam LocalAssemblerInterface the general local assembler interface
 * \tparam ExtraCtorArgs types of additional constructor arguments.
 *         Those arguments will be passed to the constructor of
 *         \c LocalAssemblerImplementation.
 *
 * The first two template parameters cannot be deduced from the arguments.
 * Therefore they always have to be provided manually.
 */
template <template <typename /* shp fct */,
                    int /* global dim */> class LocalAssemblerImplementation,
          int MinElementDim = 1, typename LocalAssemblerInterface,
          IntegrationMethodProviderOrIntegrationOrder ProviderOrOrder,
          typename... ExtraCtorArgs>
void createLocalAssemblers(
    const unsigned dimension,
    std::vector<MeshLib::Element*> const& mesh_elements,
    NumLib::LocalToGlobalIndexMap const& dof_table,
    std::vector<std::unique_ptr<LocalAssemblerInterface>>& local_assemblers,
    ProviderOrOrder const& provider_or_order,
    ExtraCtorArgs&&... extra_ctor_args)
{
    DBUG("Create local assemblers.");

    switch (dimension)
    {
        case 1:
            createLocalAssemblers<1, LocalAssemblerImplementation,
                                  MinElementDim>(
                mesh_elements, dof_table, local_assemblers, provider_or_order,
                std::forward<ExtraCtorArgs>(extra_ctor_args)...);
            break;
        case 2:
            createLocalAssemblers<2, LocalAssemblerImplementation,
                                  MinElementDim>(
                mesh_elements, dof_table, local_assemblers, provider_or_order,
                std::forward<ExtraCtorArgs>(extra_ctor_args)...);
            break;
        case 3:
            createLocalAssemblers<3, LocalAssemblerImplementation,
                                  MinElementDim>(
                mesh_elements, dof_table, local_assemblers, provider_or_order,
                std::forward<ExtraCtorArgs>(extra_ctor_args)...);
            break;
        default:
            OGS_FATAL(
                "Meshes with dimension greater than three are not supported.");
    }
}

}  // namespace ProcessLib
