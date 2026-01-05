// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "ProcessLib/Utils/LocalAssemblerFactoryTaylorHood.h"

namespace ProcessLib
{
namespace BoundaryConditionAndSourceTerm
{
template <typename LocalAssemblerInterface,
          template <typename /* shp */, typename /* lower order shp */,
                    int /* global dim */>
          class LocalAssemblerImplementation,
          int GlobalDim, typename... ConstructorArgs>
using LocalAssemblerFactoryPython = LocalAssemblerFactoryTaylorHood<
    0 /* all orders (point is 0) */, 0 /* all dimensions */,
    LocalAssemblerInterface, LocalAssemblerImplementation,
    NumLib::DefaultIntegrationMethodProvider, GlobalDim, ConstructorArgs...>;

}  // namespace BoundaryConditionAndSourceTerm
}  // namespace ProcessLib
