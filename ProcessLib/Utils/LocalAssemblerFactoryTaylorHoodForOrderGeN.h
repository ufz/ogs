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
template <typename ShapeFunction,
          typename LowerOrderShapeFunction,
          typename LocalAssemblerInterface,
          template <typename /* shp fct */,
                    typename /* lower order shp fct */,
                    typename /* int meth */,
                    int /* global dim */>
          class LocalAssemblerImplementation,
          int GlobalDim,
          typename... ConstructorArgs>
class LocalAssemblerBuilderFactoryTaylorHood
{
    using GLAF = GenericLocalAssemblerFactory<LocalAssemblerInterface,
                                              ConstructorArgs...>;
    using LocAsmIntfPtr = typename GLAF::LocAsmIntfPtr;
    using LocAsmBuilder = typename GLAF::LocAsmBuilder;

    using IntegrationMethod = typename NumLib::GaussLegendreIntegrationPolicy<
        typename ShapeFunction::MeshElement>::IntegrationMethod;

    using LocAsmImpl = LocalAssemblerImplementation<ShapeFunction,
                                                    LowerOrderShapeFunction,
                                                    IntegrationMethod,
                                                    GlobalDim>;

    LocalAssemblerBuilderFactoryTaylorHood() = delete;

public:
    /// Generates a function that creates a new local assembler of type
    /// \c LocAsmImpl.
    static LocAsmBuilder create()
    {
        return [](MeshLib::Element const& e,
                  std::size_t const local_matrix_size,
                  ConstructorArgs&&... args)
        {
            return std::make_unique<LocAsmImpl>(
                e, local_matrix_size, std::forward<ConstructorArgs>(args)...);
        };
    }
};

/**
 * Local assembler factory for Taylor-Hood elements.
 *
 * Elements/shape functions must be of order greater than or equal to \c
 * MinShapeFctOrder.
 *
 * If \c MinShapeFctOrder is 1, local assemblers are instantiated also for
 * linear mesh elements. In this case we don't have Taylor-Hood elements for
 * linear mesh elements. Instead, on linear mesh elements all shape functions
 * will have the same order (namely 1).
 */
template <int MinShapeFctOrder,
          typename LocalAssemblerInterface,
          template <typename /* shp fct */,
                    typename /* lower order shp fct */,
                    typename /* int meth */,
                    int /* global dim */>
          class LocalAssemblerImplementation,
          int GlobalDim,
          typename... ConstructorArgs>
class LocalAssemblerFactoryTaylorHoodForOrderGeN final
    : public ProcessLib::GenericLocalAssemblerFactory<LocalAssemblerInterface,
                                                      ConstructorArgs...>
{
    using Base =
        ProcessLib::GenericLocalAssemblerFactory<LocalAssemblerInterface,
                                                 ConstructorArgs...>;

    template <typename ShapeFunction, typename LowerOrderShapeFunction>
    using LocAsmBuilderFactory =
        LocalAssemblerBuilderFactoryTaylorHood<ShapeFunction,
                                               LowerOrderShapeFunction,
                                               LocalAssemblerInterface,
                                               LocalAssemblerImplementation,
                                               GlobalDim,
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

            if constexpr (ElementTraits::Element::dimension < 2)
            {
                // in OGS Taylor-Hood elements are available in 2D and 3D only
                return false;
            }

            return ElementTraits::ShapeFunction::ORDER >= MinShapeFctOrder;
        }
    };

public:
    explicit LocalAssemblerFactoryTaylorHoodForOrderGeN(
        NumLib::LocalToGlobalIndexMap const& dof_table)
        : Base(dof_table)
    {
        using EnabledElementTraits =
            decltype(BaseLib::TMP::filter<EnabledElementTraitsLagrange>(
                std::declval<IsElementEnabled>()));

        BaseLib::TMP::foreach<EnabledElementTraits>(
            [this]<typename ET>(ET*)
            {
                using Elt = typename ET::Element;
                using Shp = typename ET::ShapeFunction;
                using ShpLow = typename ET::LowerOrderShapeFunction;

                Base::_builders[std::type_index(typeid(Elt))] =
                    LocalAssemblerBuilderFactoryTaylorHood<
                        Shp,
                        ShpLow,
                        LocalAssemblerInterface,
                        LocalAssemblerImplementation,
                        GlobalDim,
                        ConstructorArgs...>::create();
            });
    }
};

/// HM processes in OGS are defined for linear and higher order elements.
template <typename LocalAssemblerInterface,
          template <typename, typename, typename, int>
          class LocalAssemblerImplementation,
          int GlobalDim,
          typename... ConstructorArgs>
using LocalAssemblerFactoryHM =
    LocalAssemblerFactoryTaylorHoodForOrderGeN<1,
                                               LocalAssemblerInterface,
                                               LocalAssemblerImplementation,
                                               GlobalDim,
                                               ConstructorArgs...>;

/// Stokes flow in OGS is defined for higher order elements only.
template <typename LocalAssemblerInterface,
          template <typename, typename, typename, int>
          class LocalAssemblerImplementation,
          int GlobalDim,
          typename... ConstructorArgs>
using LocalAssemblerFactoryStokes =
    LocalAssemblerFactoryTaylorHoodForOrderGeN<2,
                                               LocalAssemblerInterface,
                                               LocalAssemblerImplementation,
                                               GlobalDim,
                                               ConstructorArgs...>;

}  // namespace ProcessLib
