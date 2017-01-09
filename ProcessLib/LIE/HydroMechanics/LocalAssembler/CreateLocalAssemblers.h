/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <vector>

#include <logog/include/logog.hpp>

#include "NumLib/DOF/LocalToGlobalIndexMap.h"

#include "LocalDataInitializer.h"

namespace ProcessLib
{
namespace LIE
{
namespace HydroMechanics
{
namespace detail
{
template <unsigned GlobalDim,
          template <typename, typename, typename, unsigned>
          class LocalAssemblerMatrixImplementation,
          template <typename, typename, typename, unsigned>
          class LocalAssemblerMatrixNearFractureImplementation,
          template <typename, typename, typename, unsigned>
          class LocalAssemblerFractureImplementation,
          typename LocalAssemblerInterface, typename... ExtraCtorArgs>
void createLocalAssemblers(
    NumLib::LocalToGlobalIndexMap const& dof_table,
    const unsigned shapefunction_order,
    std::vector<MeshLib::Element*> const& mesh_elements,
    std::vector<std::unique_ptr<LocalAssemblerInterface>>& local_assemblers,
    ExtraCtorArgs&&... extra_ctor_args)
{
    // Shape matrices initializer
    using LocalDataInitializer =
        LocalDataInitializer<LocalAssemblerInterface,
                             LocalAssemblerMatrixImplementation,
                             LocalAssemblerMatrixNearFractureImplementation,
                             LocalAssemblerFractureImplementation,
                             GlobalDim, ExtraCtorArgs...>;

    DBUG("Create local assemblers for HydroMechanics with LIE.");
    // Populate the vector of local assemblers.
    local_assemblers.resize(mesh_elements.size());

    LocalDataInitializer initializer(dof_table, shapefunction_order);

    DBUG("Calling local assembler builder for all mesh elements.");
    GlobalExecutor::transformDereferenced(
        initializer, mesh_elements, local_assemblers,
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
template <unsigned GlobalDim,
          template <typename, typename, typename, unsigned>
          class LocalAssemblerMatrixImplementation,
          template <typename, typename, typename, unsigned>
          class LocalAssemblerMatrixNearFractureImplementation,
          template <typename, typename, typename, unsigned>
          class LocalAssemblerFractureImplementation,
          typename LocalAssemblerInterface, typename... ExtraCtorArgs>
void createLocalAssemblers(
    std::vector<MeshLib::Element*> const& mesh_elements,
    NumLib::LocalToGlobalIndexMap const& dof_table,
    const unsigned shapefunction_order,
    std::vector<std::unique_ptr<LocalAssemblerInterface>>& local_assemblers,
    ExtraCtorArgs&&... extra_ctor_args)
{
    detail::createLocalAssemblers<GlobalDim,
                                  LocalAssemblerMatrixImplementation,
                                  LocalAssemblerMatrixNearFractureImplementation,
                                  LocalAssemblerFractureImplementation>(
        dof_table, shapefunction_order, mesh_elements, local_assemblers,
        std::forward<ExtraCtorArgs>(extra_ctor_args)...);
}

}  // namespace HydroMechanics
}  // namespace LIE
}  // namespace ProcessLib
