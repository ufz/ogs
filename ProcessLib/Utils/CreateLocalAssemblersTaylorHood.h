/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
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
template <template <typename /*LocalAssemblerInterface*/,
                    template <typename /* shp fct */,
                              typename /* lower order shp fct */,
                              typename /* int meth */, int /* global dim */>
                    class /*LocalAssemblerImplementation*/,
                    int /* global dim */, typename... /*ConstructorArgs*/>
          class LocalAssemblerFactory,
          int GlobalDim,
          template <typename /* shp fct */, typename /* lower order shp fct */,
                    typename /* int meth */, int /* global dim */>
          class LocalAssemblerImplementation,
          typename LocalAssemblerInterface, typename... ExtraCtorArgs>
void createLocalAssemblersTaylorHood(
    std::vector<MeshLib::Element*> const& mesh_elements,
    NumLib::LocalToGlobalIndexMap const& dof_table,
    std::vector<std::unique_ptr<LocalAssemblerInterface>>& local_assemblers,
    ExtraCtorArgs&&... extra_ctor_args)
{
    using LocAsmFac = LocalAssemblerFactory<LocalAssemblerInterface,
                                            LocalAssemblerImplementation,
                                            GlobalDim, ExtraCtorArgs...>;

    DBUG("Create local assemblers.");

    LocAsmFac factory(dof_table);
    local_assemblers.resize(mesh_elements.size());

    DBUG("Calling local assembler builder for all mesh elements.");
    GlobalExecutor::transformDereferenced(
        factory, mesh_elements, local_assemblers,
        std::forward<ExtraCtorArgs>(extra_ctor_args)...);
}

}  // namespace detail

template <int GlobalDim,
          template <typename /* shp fct */, typename /* lower order shp fct */,
                    typename /* int meth */, int /* global dim */>
          class LocalAssemblerImplementation,
          typename LocalAssemblerInterface, typename... ExtraCtorArgs>
void createLocalAssemblersHM(
    std::vector<MeshLib::Element*> const& mesh_elements,
    NumLib::LocalToGlobalIndexMap const& dof_table,
    std::vector<std::unique_ptr<LocalAssemblerInterface>>& local_assemblers,
    ExtraCtorArgs&&... extra_ctor_args)
{
    detail::createLocalAssemblersTaylorHood<LocalAssemblerFactoryHM, GlobalDim,
                                            LocalAssemblerImplementation,
                                            LocalAssemblerInterface>(
        mesh_elements, dof_table, local_assemblers,
        std::forward<ExtraCtorArgs>(extra_ctor_args)...);
}

template <int GlobalDim,
          template <typename /* shp fct */, typename /* lower order shp fct */,
                    typename /* int meth */, int /* global dim */>
          class LocalAssemblerImplementation,
          typename LocalAssemblerInterface, typename... ExtraCtorArgs>
void createLocalAssemblersStokes(
    std::vector<MeshLib::Element*> const& mesh_elements,
    NumLib::LocalToGlobalIndexMap const& dof_table,
    std::vector<std::unique_ptr<LocalAssemblerInterface>>& local_assemblers,
    ExtraCtorArgs&&... extra_ctor_args)
{
    detail::createLocalAssemblersTaylorHood<
        LocalAssemblerFactoryStokes, GlobalDim, LocalAssemblerImplementation,
        LocalAssemblerInterface>(
        mesh_elements, dof_table, local_assemblers,
        std::forward<ExtraCtorArgs>(extra_ctor_args)...);
}

}  // namespace ProcessLib
