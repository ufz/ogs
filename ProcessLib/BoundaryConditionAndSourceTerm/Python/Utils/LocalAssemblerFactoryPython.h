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

#include "ProcessLib/Utils/EnabledElements.h"
#include "ProcessLib/Utils/LocalAssemblerFactoryTaylorHood.h"

namespace ProcessLib
{
namespace BoundaryConditionAndSourceTerm
{
template <typename LocalAssemblerInterface,
          template <typename /* shp */, typename /* lower order shp */,
                    typename /* int meth */, int /* global dim */>
          class LocalAssemblerImplementation,
          int GlobalDim, typename... ConstructorArgs>
using LocalAssemblerFactoryPython = LocalAssemblerFactoryTaylorHood<
    0 /* all orders (point is 0) */, 0 /* all dimensions */,
    LocalAssemblerInterface, LocalAssemblerImplementation, GlobalDim,
    ConstructorArgs...>;

}  // namespace BoundaryConditionAndSourceTerm
}  // namespace ProcessLib
