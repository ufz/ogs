/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <vector>

#include "BaseLib/Logging.h"
#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "ProcessLib/Utils/LocalAssemblerFactoryForDimGreaterEqualN.h"

namespace ProcessLib
{
namespace SmallDeformation
{
namespace detail
{
template <int GlobalDim,
          template <typename /* shp fct */, int /* global dim */>
          class LocalAssemblerImplementation,
          typename LocalAssemblerInterface,
          IntegrationMethodProviderOrIntegrationOrder ProviderOrOrder,
          typename... ExtraCtorArgs>
void createLocalAssemblers(
    NumLib::LocalToGlobalIndexMap const& dof_table,
    std::vector<MeshLib::Element*> const& mesh_elements,
    std::vector<std::unique_ptr<LocalAssemblerInterface>>& local_assemblers,
    ProviderOrOrder const& provider_or_order,
    ExtraCtorArgs&&... extra_ctor_args)
{
    DBUG("Create local assemblers.");

    auto const& integration_method_provider =
        getIntegrationMethodProvider(provider_or_order);

    using IntMethProv =
        std::remove_cvref_t<decltype(integration_method_provider)>;
    using LocAsmFactory = ProcessLib::LocalAssemblerFactorySD<
        LocalAssemblerInterface, LocalAssemblerImplementation, IntMethProv,
        GlobalDim, ExtraCtorArgs...>;

    LocAsmFactory factory(dof_table, integration_method_provider);
    local_assemblers.resize(mesh_elements.size());

    DBUG("Calling local assembler builder for all mesh elements.");
    GlobalExecutor::transformDereferenced(
        factory, mesh_elements, local_assemblers,
        std::forward<ExtraCtorArgs>(extra_ctor_args)...);
}

}  // namespace detail

/*! Creates local assemblers for each element of the given \c mesh.
 *
 * \tparam LocalAssemblerImplementation the individual local assembler type
 * \tparam LocalAssemblerInterface the general local assembler interface
 * \tparam ExtraCtorArgs types of additional constructor arguments.
 *         Those arguments will be passed to the constructor of
 *         \c LocalAssemblerImplementation.
 *
 * The first two template parameters cannot be deduced from the arguments.
 * Therefore they always have to be provided manually.
 */
template <int GlobalDim,
          template <typename /* shp fct */, int /* global dim */>
          class LocalAssemblerImplementation,
          typename LocalAssemblerInterface,
          IntegrationMethodProviderOrIntegrationOrder ProviderOrOrder,
          typename... ExtraCtorArgs>
void createLocalAssemblers(
    std::vector<MeshLib::Element*> const& mesh_elements,
    NumLib::LocalToGlobalIndexMap const& dof_table,
    std::vector<std::unique_ptr<LocalAssemblerInterface>>& local_assemblers,
    ProviderOrOrder const& provider_or_order,
    ExtraCtorArgs&&... extra_ctor_args)
{
    DBUG("Create local assemblers.");

    detail::createLocalAssemblers<GlobalDim, LocalAssemblerImplementation>(
        dof_table, mesh_elements, local_assemblers, provider_or_order,
        std::forward<ExtraCtorArgs>(extra_ctor_args)...);
}
}  // namespace SmallDeformation

}  // namespace ProcessLib
