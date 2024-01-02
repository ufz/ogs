/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <vector>

#include "BaseLib/Logging.h"
#include "LocalAssemblerFactoryTaylorHood.h"
#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "NumLib/Fem/Integration/IntegrationMethodProvider.h"

namespace ProcessLib
{
namespace detail
{
/*! Creates local assemblers for each element of the given \c mesh.
 *
 * \tparam LocalAssemblerFactory the factory that will instantiate the local
 *         assemblers
 * \tparam LocalAssemblerImplementation the individual local assembler type
 * \tparam LocalAssemblerInterface the general local assembler interface
 * \tparam ExtraCtorArgs types of additional constructor arguments.
 *         Those arguments will be passed to the constructor of
 *         \c LocalAssemblerImplementation.
 */
template <template <typename /* LocalAssemblerInterface */,
                    template <typename /* shp fct */,
                              typename /* lower order shp fct */,
                              int /* global dim */>
                    class /* LocalAssemblerImplementation */,
                    class /* IntegrationMethodProvider */, int /* global dim */,
                    typename... /* ConstructorArgs */>
          class LocalAssemblerFactory,
          int GlobalDim,
          template <typename /* shp fct */, typename /* lower order shp fct */,
                    int /* global dim */>
          class LocalAssemblerImplementation,
          typename LocalAssemblerInterface,
          IntegrationMethodProviderOrIntegrationOrder ProviderOrOrder,
          typename... ExtraCtorArgs>
void createLocalAssemblersTaylorHood(
    std::vector<MeshLib::Element*> const& mesh_elements,
    NumLib::LocalToGlobalIndexMap const& dof_table,
    std::vector<std::unique_ptr<LocalAssemblerInterface>>& local_assemblers,
    ProviderOrOrder const& provider_or_order,
    ExtraCtorArgs&&... extra_ctor_args)
{
    DBUG("Create local assemblers.");

    auto const& integration_method_provider =
        getIntegrationMethodProvider(provider_or_order);

    using IntegrationMethodProvider =
        std::remove_cvref_t<decltype(integration_method_provider)>;
    using LocAsmFac = LocalAssemblerFactory<
        LocalAssemblerInterface, LocalAssemblerImplementation,
        IntegrationMethodProvider, GlobalDim, ExtraCtorArgs...>;

    LocAsmFac factory(dof_table, integration_method_provider);
    local_assemblers.resize(mesh_elements.size());

    DBUG("Calling local assembler builder for all mesh elements.");
    GlobalExecutor::transformDereferenced(
        factory, mesh_elements, local_assemblers,
        std::forward<ExtraCtorArgs>(extra_ctor_args)...);
}

}  // namespace detail

template <int GlobalDim,
          template <typename /* shp fct */, typename /* lower order shp fct */,
                    int /* global dim */>
          class LocalAssemblerImplementation,
          typename LocalAssemblerInterface,
          IntegrationMethodProviderOrIntegrationOrder ProviderOrOrder,
          typename... ExtraCtorArgs>
void createLocalAssemblersHM(
    std::vector<MeshLib::Element*> const& mesh_elements,
    NumLib::LocalToGlobalIndexMap const& dof_table,
    std::vector<std::unique_ptr<LocalAssemblerInterface>>& local_assemblers,
    ProviderOrOrder const& provider_or_order,
    ExtraCtorArgs&&... extra_ctor_args)
{
    detail::createLocalAssemblersTaylorHood<LocalAssemblerFactoryHM, GlobalDim,
                                            LocalAssemblerImplementation,
                                            LocalAssemblerInterface>(
        mesh_elements, dof_table, local_assemblers, provider_or_order,
        std::forward<ExtraCtorArgs>(extra_ctor_args)...);
}

template <int GlobalDim,
          template <typename /* shp fct */, typename /* lower order shp fct */,
                    int /* global dim */>
          class LocalAssemblerImplementation,
          typename LocalAssemblerInterface,
          IntegrationMethodProviderOrIntegrationOrder ProviderOrOrder,
          typename... ExtraCtorArgs>
void createLocalAssemblersStokes(
    std::vector<MeshLib::Element*> const& mesh_elements,
    NumLib::LocalToGlobalIndexMap const& dof_table,
    std::vector<std::unique_ptr<LocalAssemblerInterface>>& local_assemblers,
    ProviderOrOrder const& provider_or_order,
    ExtraCtorArgs&&... extra_ctor_args)
{
    detail::createLocalAssemblersTaylorHood<
        LocalAssemblerFactoryStokes, GlobalDim, LocalAssemblerImplementation,
        LocalAssemblerInterface>(
        mesh_elements, dof_table, local_assemblers, provider_or_order,
        std::forward<ExtraCtorArgs>(extra_ctor_args)...);
}

}  // namespace ProcessLib
