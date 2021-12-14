/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
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
template <int MinElementDim,
          typename LocalAssemblerInterface,
          template <typename, typename, int>
          class LocalAssemblerImplementation,
          int GlobalDim,
          typename... ConstructorArgs>
class LocalAssemberFactoryForDimGeN final
    : public GenericLocalAssemblerFactory<LocalAssemblerInterface,
                                          ConstructorArgs...>
{
    using Base = GenericLocalAssemblerFactory<LocalAssemblerInterface,
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
    explicit LocalAssemberFactoryForDimGeN(
        NumLib::LocalToGlobalIndexMap const& dof_table)
        : Base{dof_table}
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
                    LocalAssemblerBuilderFactory<ShapeFunction,
                                                 LocalAssemblerInterface,
                                                 LocalAssemblerImplementation,
                                                 GlobalDim,
                                                 ConstructorArgs...>::create();
            });
    }
};

/// By default processes in OGS are defined in 1D, 2D and 3D.
template <typename LocalAssemblerInterface,
          template <typename, typename, int>
          class LocalAssemblerImplementation,
          int GlobalDim,
          typename... ConstructorArgs>
using LocalAssemberFactory =
    LocalAssemberFactoryForDimGeN<1,
                                  LocalAssemblerInterface,
                                  LocalAssemblerImplementation,
                                  GlobalDim,
                                  ConstructorArgs...>;

/// Mechanics processes in OGS are defined in 2D and 3D only.
template <typename LocalAssemblerInterface,
          template <typename, typename, int>
          class LocalAssemblerImplementation,
          int GlobalDim,
          typename... ConstructorArgs>
using LocalAssemberFactorySD =
    LocalAssemberFactoryForDimGeN<2,
                                  LocalAssemblerInterface,
                                  LocalAssemblerImplementation,
                                  GlobalDim,
                                  ConstructorArgs...>;
}  // namespace ProcessLib
