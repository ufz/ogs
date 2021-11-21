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

#include "ProcessLib/Utils/EnabledElements.h"
#include "ProcessLib/Utils/GenericLocalDataInitializer.h"

namespace ProcessLib
{
namespace BoundaryConditionAndSourceTerm
{
template <typename LocalAssemblerInterface,
          template <typename, typename, int>
          class LocalAssemblerImplementation,
          int GlobalDim,
          typename... ConstructorArgs>
class LocalDataInitializer final
    : public ProcessLib::GenericLocalDataInitializer<LocalAssemblerInterface,
                                                     ConstructorArgs...>
{
    using Base =
        ProcessLib::GenericLocalDataInitializer<LocalAssemblerInterface,
                                                ConstructorArgs...>;

    template <typename ShapeFunction>
    using LocAsmBuilderFactory =
        ProcessLib::LocalAssemblerBuilderFactory<ShapeFunction,
                                                 LocalAssemblerInterface,
                                                 LocalAssemblerImplementation,
                                                 GlobalDim,
                                                 ConstructorArgs...>;

    struct HasSuitableDimension
    {
        template <typename ElementTraits>
        constexpr bool operator()(ElementTraits*) const
        {
            return GlobalDim >= ElementTraits::ShapeFunction::DIM;
        }
    };

    struct Is2ndOrderElementOfSuitableDimension
    {
        template <typename ElementTraits>
        constexpr bool operator()(ElementTraits*) const
        {
            if constexpr (GlobalDim < ElementTraits::ShapeFunction::DIM)
            {
                return false;
            }
            return ElementTraits::ShapeFunction::ORDER == 2 ||
                   // points are needed for 2nd order, too
                   std::is_same_v<MeshLib::Point,
                                  typename ElementTraits::Element>;
        }
    };

public:
    LocalDataInitializer(NumLib::LocalToGlobalIndexMap const& dof_table,
                         const unsigned shapefunction_order)
        : Base(dof_table)
    {
        if (shapefunction_order < 1 || 2 < shapefunction_order)
        {
            OGS_FATAL("The given shape function order {:d} is not supported",
                      shapefunction_order);
        }

        if (shapefunction_order == 1)
        {
            // 1st order is enabled on all elements with suitable dimension
            using EnabledElementTraits =
                decltype(BaseLib::TMP::filter<EnabledElementTraitsLagrange>(
                    std::declval<HasSuitableDimension>()));

            BaseLib::TMP::foreach<EnabledElementTraits>(
                [this]<typename ET>(ET*)
                {
                    using MeshElement = typename ET::Element;
                    // this will use linear shape functions on higher order
                    // elements and the linear shape function on linear elements
                    using LowerOrderShapeFunction =
                        typename ET::LowerOrderShapeFunction;
                    _builder[std::type_index(typeid(MeshElement))] =
                        LocAsmBuilderFactory<LowerOrderShapeFunction>::create();
                });
        }
        else if (shapefunction_order == 2)
        {
            // 2nd order only on 2nd order elements
            using EnabledElementTraits =
                decltype(BaseLib::TMP::filter<EnabledElementTraitsLagrange>(
                    std::declval<Is2ndOrderElementOfSuitableDimension>()));

            BaseLib::TMP::foreach<EnabledElementTraits>(
                [this]<typename ET>(ET*)
                {
                    using MeshElement = typename ET::Element;
                    using ShapeFunction2ndOrder = typename ET::ShapeFunction;
                    _builder[std::type_index(typeid(MeshElement))] =
                        LocAsmBuilderFactory<ShapeFunction2ndOrder>::create();
                });
        }
    }
};

}  // namespace BoundaryConditionAndSourceTerm
}  // namespace ProcessLib
