/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <vector>

#include <logog/include/logog.hpp>

#include "LocalDataInitializer.h"
#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "ProcessLib/HeatTransportBHE/BHE/BHEAbstract.h"

namespace ProcessLib
{
namespace HeatTransportBHE
{
namespace detail
{
template <
    int GlobalDim,
    template <typename, typename, int> class LocalAssemblerSoilImplementation,
    template <typename, typename, int> class LocalAssemblerBHEImplementation,
    typename LocalAssemblerInterface, typename... ExtraCtorArgs>
void createLocalAssemblers(
    NumLib::LocalToGlobalIndexMap const& dof_table,
    std::vector<MeshLib::Element*> const& mesh_elements,
    std::vector<std::unique_ptr<LocalAssemblerInterface>>& local_assemblers,
    const std::vector<std::vector<std::size_t>>& vec_ele_connected_BHE_IDs,
    const std::vector<std::unique_ptr<BHE::BHEAbstract>>& vec_BHE_property,
    ExtraCtorArgs&&... extra_ctor_args)
{
    // Shape matrices initializer
    using LocalDataInitializer = LocalDataInitializer<
        LocalAssemblerInterface, LocalAssemblerSoilImplementation,
        LocalAssemblerBHEImplementation, 3, ExtraCtorArgs...>;

    DBUG("Create local assemblers for the HeatTransportBHE process.");
    // Populate the vector of local assemblers.
    local_assemblers.resize(mesh_elements.size());

    LocalDataInitializer initializer(dof_table);

    DBUG("Calling local assembler builder for all mesh elements.");
    GlobalExecutor::transformDereferenced(
        initializer, mesh_elements, local_assemblers, vec_ele_connected_BHE_IDs,
        vec_BHE_property, std::forward<ExtraCtorArgs>(extra_ctor_args)...);
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
template <
    int GlobalDim,
    template <typename, typename, int> class LocalAssemblerSoilImplementation,
    template <typename, typename, int> class LocalAssemblerBHEImplementation,
    typename LocalAssemblerInterface, typename... ExtraCtorArgs>
void createLocalAssemblers(
    std::vector<MeshLib::Element*> const& mesh_elements,
    NumLib::LocalToGlobalIndexMap const& dof_table,
    std::vector<std::unique_ptr<LocalAssemblerInterface>>& local_assemblers,
    ExtraCtorArgs&&... extra_ctor_args)
{
    DBUG("Create local assemblers for the HeatTransportBHE process.");

    detail::createLocalAssemblers<GlobalDim, LocalAssemblerSoilImplementation,
                                  LocalAssemblerBHEImplementation>(
        dof_table, mesh_elements, local_assemblers,
        std::forward<ExtraCtorArgs>(extra_ctor_args)...);
}
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
