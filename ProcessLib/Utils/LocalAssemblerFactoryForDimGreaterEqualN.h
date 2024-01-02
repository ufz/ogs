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

#include "EnabledElements.h"
#include "GenericLocalAssemblerFactory.h"

namespace ProcessLib
{
/**
 * Creates local assemblers for mesh elements with dimension greater than or
 * equal to \c MinElementDim.
 *
 * This local assembler factory supports a single type of shape functions, i.e.,
 * all primary variables are discretized with the same shape function.
 */
template <int MinElementDim,
          typename LocalAssemblerInterface,
          template <typename /* shp fct */, int /* global dim */>
          class LocalAssemblerImplementation,
          NumLib::IntegrationMethodProvider IntegrationMethodProvider,
          int GlobalDim,
          typename... ConstructorArgs>
class LocalAssemblerFactoryForDimGreaterEqualN final
    : public GenericLocalAssemblerFactory<LocalAssemblerInterface,
                                          IntegrationMethodProvider,
                                          ConstructorArgs...>
{
    using Base = GenericLocalAssemblerFactory<LocalAssemblerInterface,
                                              IntegrationMethodProvider,
                                              ConstructorArgs...>;

    struct IsElementEnabled
    {
        template <typename ElementTraits>
        constexpr bool operator()(ElementTraits*) const
        {
            if constexpr (GlobalDim < ElementTraits::ShapeFunction::DIM)
            {
                return false;
            }

            return ElementTraits::Element::dimension >= MinElementDim;
        }
    };

public:
    explicit LocalAssemblerFactoryForDimGreaterEqualN(
        NumLib::LocalToGlobalIndexMap const& dof_table,
        IntegrationMethodProvider const& integration_method_provider)
        : Base{dof_table, integration_method_provider}
    {
        using EnabledElementTraits =
            decltype(BaseLib::TMP::filter<EnabledElementTraitsLagrange>(
                std::declval<IsElementEnabled>()));

        BaseLib::TMP::foreach<EnabledElementTraits>(
            [this]<typename ET>(ET*)
            {
                using MeshElement = typename ET::Element;
                using ShapeFunction = typename ET::ShapeFunction;
                Base::_builders[std::type_index(typeid(MeshElement))] =
                    LocalAssemblerBuilderFactory<
                        ShapeFunction,
                        LocalAssemblerInterface,
                        LocalAssemblerImplementation,
                        IntegrationMethodProvider,
                        GlobalDim,
                        ConstructorArgs...>::template create<MeshElement>();
            });
    }
};

/// By default processes in OGS are defined in 1D, 2D and 3D.
template <typename LocalAssemblerInterface,
          template <typename /* shp fct */, int /* global dim */>
          class LocalAssemblerImplementation,
          NumLib::IntegrationMethodProvider IntegrationMethodProvider,
          int GlobalDim,
          typename... ConstructorArgs>
using LocalAssemblerFactory =
    LocalAssemblerFactoryForDimGreaterEqualN<1,
                                             LocalAssemblerInterface,
                                             LocalAssemblerImplementation,
                                             IntegrationMethodProvider,
                                             GlobalDim,
                                             ConstructorArgs...>;

/// Mechanics processes in OGS are defined in 2D and 3D only.
template <typename LocalAssemblerInterface,
          template <typename /* shp fct */, int /* global dim */>
          class LocalAssemblerImplementation,
          NumLib::IntegrationMethodProvider IntegrationMethodProvider,
          int GlobalDim,
          typename... ConstructorArgs>
using LocalAssemblerFactorySD =
    LocalAssemblerFactoryForDimGreaterEqualN<2,
                                             LocalAssemblerInterface,
                                             LocalAssemblerImplementation,
                                             IntegrationMethodProvider,
                                             GlobalDim,
                                             ConstructorArgs...>;
}  // namespace ProcessLib
