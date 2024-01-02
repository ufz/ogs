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

#include "NumLib/Fem/Integration/IntegrationMethodProvider.h"
#include "ProcessLib/Utils/EnabledElements.h"
#include "ProcessLib/Utils/GenericLocalAssemblerFactory.h"

namespace ProcessLib
{
namespace BoundaryConditionAndSourceTerm
{
template <typename LocalAssemblerInterface,
          template <typename /* shp fct */, int /* global dim */>
          class LocalAssemblerImplementation,
          int GlobalDim,
          typename... ConstructorArgs>
class LocalAssemblerFactory final
    : public ProcessLib::GenericLocalAssemblerFactory<
          LocalAssemblerInterface,
          NumLib::DefaultIntegrationMethodProvider,
          ConstructorArgs...>
{
    using Base = ProcessLib::GenericLocalAssemblerFactory<
        LocalAssemblerInterface,
        NumLib::DefaultIntegrationMethodProvider,
        ConstructorArgs...>;

    template <typename ShapeFunction>
    using LocAsmBuilderFactory = ProcessLib::LocalAssemblerBuilderFactory<
        ShapeFunction,
        LocalAssemblerInterface,
        LocalAssemblerImplementation,
        NumLib::DefaultIntegrationMethodProvider,
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
    LocalAssemblerFactory(NumLib::LocalToGlobalIndexMap const& dof_table,
                          NumLib::DefaultIntegrationMethodProvider const&
                              integration_method_provider,
                          const unsigned shapefunction_order)
        : Base{dof_table, integration_method_provider}
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
                    Base::_builders[std::type_index(typeid(MeshElement))] =
                        LocAsmBuilderFactory<LowerOrderShapeFunction>::
                            template create<MeshElement>();
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
                    Base::_builders[std::type_index(typeid(MeshElement))] =
                        LocAsmBuilderFactory<ShapeFunction2ndOrder>::
                            template create<MeshElement>();
                });
        }
    }
};

}  // namespace BoundaryConditionAndSourceTerm
}  // namespace ProcessLib
