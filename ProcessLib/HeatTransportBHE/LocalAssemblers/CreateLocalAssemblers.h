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

#include "LocalDataInitializer.h"
#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "ProcessLib/HeatTransportBHE/BHE/BHETypes.h"

namespace ProcessLib
{
namespace HeatTransportBHE
{
namespace detail
{
template <template <typename /* shp fct */>
          class LocalAssemblerSoilImplementation,
          template <typename /* shp fct */, typename /* bhe type */>
          class LocalAssemblerBHEImplementation,
          typename LocalAssemblerInterface, typename... ExtraCtorArgs>
void createLocalAssemblers(
    NumLib::LocalToGlobalIndexMap const& dof_table,
    std::vector<MeshLib::Element*> const& mesh_elements,
    std::vector<std::unique_ptr<LocalAssemblerInterface>>& local_assemblers,
    NumLib::IntegrationOrder const integration_order,
    std::unordered_map<std::size_t, BHE::BHETypes*> const& element_to_bhe_map,
    ExtraCtorArgs&&... extra_ctor_args)
{
    // Shape matrices initializer
    using LocalDataInitializer =
        LocalDataInitializer<LocalAssemblerInterface,
                             LocalAssemblerSoilImplementation,
                             LocalAssemblerBHEImplementation, ExtraCtorArgs...>;

    DBUG("Create local assemblers for the HeatTransportBHE process.");
    // Populate the vector of local assemblers.
    local_assemblers.resize(mesh_elements.size());

    LocalDataInitializer initializer(dof_table, integration_order);

    DBUG("Calling local assembler builder for all mesh elements.");
    GlobalExecutor::transformDereferenced(
        initializer, mesh_elements, local_assemblers, element_to_bhe_map,
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
template <template <typename /* shp fct */>
          class LocalAssemblerSoilImplementation,
          template <typename /* shp fct */, typename /* bhe type */>
          class LocalAssemblerBHEImplementation,
          typename LocalAssemblerInterface, typename... ExtraCtorArgs>
void createLocalAssemblers(
    std::vector<MeshLib::Element*> const& mesh_elements,
    NumLib::LocalToGlobalIndexMap const& dof_table,
    std::vector<std::unique_ptr<LocalAssemblerInterface>>& local_assemblers,
    NumLib::IntegrationOrder const integration_order,
    ExtraCtorArgs&&... extra_ctor_args)
{
    DBUG("Create local assemblers for the HeatTransportBHE process.");

    detail::createLocalAssemblers<LocalAssemblerSoilImplementation,
                                  LocalAssemblerBHEImplementation>(
        dof_table, mesh_elements, local_assemblers, integration_order,
        std::forward<ExtraCtorArgs>(extra_ctor_args)...);
}
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
