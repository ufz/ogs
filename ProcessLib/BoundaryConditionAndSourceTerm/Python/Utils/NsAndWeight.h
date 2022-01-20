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

#include <Eigen/Core>

#include "MeshLib/Elements/Element.h"
#include "NumLib/Fem/InitShapeMatrices.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"

namespace ProcessLib::BoundaryConditionAndSourceTerm::Python
{

/**
 * Holds shape matrices and integration weight.
 *
 * This struct facilitates the handling of shape functions of different order
 * defined on the same mesh element, e.g., for Taylor-Hood elements.
 *
 * \note The Shape function template parameters are necessary only to
 * distinguish between different specializations of this class.
 */
template <typename ShapeFunction, typename LowerOrderShapeFunction,
          typename ShapeMatrix, typename LowerOrderShapeMatrix>
struct NsAndWeight
{
    static_assert(
        ShapeFunction::ORDER > LowerOrderShapeFunction::ORDER,
        "Shape function orders must be different from each other. For same "
        "orders the class template specialization below should be used.");

    static_assert(ShapeFunction::ORDER == 2);
    static_assert(LowerOrderShapeFunction::ORDER < 2);

    NsAndWeight(ShapeMatrix&& N1, LowerOrderShapeMatrix&& N2,
                double const weight)
        : N_higher_{std::move(N1)}, N_lower_{std::move(N2)}, weight_{weight}
    {
    }

    double weight() const { return weight_; }

    Eigen::Ref<const Eigen::RowVectorXd> N(unsigned order) const
    {
        if (order < 2)
        {
            return N_lower_;
        }

        if (order == 2)
        {
            return N_higher_;
        }

        OGS_FATAL(
            "Only shape functions up to order 2 are supported currently. You "
            "have requested order {}. There might be an error in the OGS "
            "project file.",
            order);
    }

    ShapeMatrix const& NHigherOrder() const { return N_higher_; }

private:
    ShapeMatrix N_higher_;
    LowerOrderShapeMatrix N_lower_;
    double weight_;
};

/// Specialization if shape function and lower order shape function are the
/// same.
template <typename ShapeFunction, typename ShapeMatrix>
struct NsAndWeight<ShapeFunction, ShapeFunction, ShapeMatrix, ShapeMatrix>
{
    static_assert(ShapeFunction::ORDER < 2,
                  "For higher order shape functions the general case of this "
                  "class template should be used. See above.");

    NsAndWeight(ShapeMatrix&& N, double const weight)
        : N_{std::move(N)}, weight_{weight}
    {
    }

    Eigen::Ref<const Eigen::RowVectorXd> N(unsigned order) const
    {
        if (order >= 2)
        {
            OGS_FATAL(
                "Only shape functions of order < 2 are available in the "
                "current setting. You have requested order {}. Maybe there is "
                "an error in the OGS project file.",
                order);
        }

        return N_;
    }

    ShapeMatrix const& NHigherOrder() const { return N_; }

    double weight() const { return weight_; }

private:
    ShapeMatrix N_;
    double weight_;
};

/// Collects common type aliases needed when working with NsAndWeight.
template <typename ShapeFunction, typename LowerOrderShapeFunction,
          int GlobalDim>
struct NsAndWeightsTraits
{
    using ShapeMatrixPolicy = ShapeMatrixPolicyType<ShapeFunction, GlobalDim>;
    using ShapeMatrix = typename ShapeMatrixPolicy::ShapeMatrices::ShapeType;

    using LowerOrderShapeMatrixPolicy =
        ShapeMatrixPolicyType<LowerOrderShapeFunction, GlobalDim>;
    using LowerOrderShapeMatrix =
        typename LowerOrderShapeMatrixPolicy::ShapeMatrices::ShapeType;

    using NsAndWeight =
        ProcessLib::BoundaryConditionAndSourceTerm::Python::NsAndWeight<
            ShapeFunction, LowerOrderShapeFunction, ShapeMatrix,
            LowerOrderShapeMatrix>;
};

/**
 * Computes shape matrices and integration weights for all integration points in
 * the given mesh element.
 *
 * \note This is an extension of
 * GenericNaturalBoundaryConditionLocalAssembler::initNsAndWeights().
 */
template <typename ShapeFunction, typename LowerOrderShapeFunction,
          int GlobalDim, typename IntegrationMethod>
auto computeNsAndWeights(MeshLib::Element const& element,
                         bool const is_axially_symmetric,
                         IntegrationMethod const& integration_method)
{
    using Traits =
        NsAndWeightsTraits<ShapeFunction, LowerOrderShapeFunction, GlobalDim>;
    using VecOfNsAndWeight = std::vector<typename Traits::NsAndWeight>;

    VecOfNsAndWeight nss_and_weights;
    nss_and_weights.reserve(integration_method.getNumberOfPoints());

    auto sms =
        NumLib::initShapeMatrices<ShapeFunction,
                                  typename Traits::ShapeMatrixPolicy, GlobalDim,
                                  NumLib::ShapeMatrixType::N_J>(
            element, is_axially_symmetric, integration_method);

    if constexpr (std::is_same_v<ShapeFunction, LowerOrderShapeFunction>)
    {
        static_assert(ShapeFunction::ORDER < 2,
                      "We do not expect higher order shape functions here. "
                      "Something must have gone terribly wrong.");

        for (unsigned ip = 0; ip < sms.size(); ++ip)
        {
            auto& sm = sms[ip];
            double const w =
                sm.detJ * sm.integralMeasure *
                integration_method.getWeightedPoint(ip).getWeight();

            nss_and_weights.emplace_back(std::move(sm.N), w);
        }
    }
    else
    {
        auto sms_lower = NumLib::initShapeMatrices<
            LowerOrderShapeFunction,
            typename Traits::LowerOrderShapeMatrixPolicy, GlobalDim,
            NumLib::ShapeMatrixType::N>(element, is_axially_symmetric,
                                        integration_method);

        for (unsigned ip = 0; ip < sms.size(); ++ip)
        {
            auto& sm = sms[ip];

            // Note: we use det(J) of the higher order shape function. For
            // linear geometries (i.e., higher order elements but with flat
            // surfaces) there should be no difference.
            double const w =
                sm.detJ * sm.integralMeasure *
                integration_method.getWeightedPoint(ip).getWeight();

            nss_and_weights.emplace_back(std::move(sm.N),
                                         std::move(sms_lower[ip].N), w);
        }
    }

    return nss_and_weights;
}

}  // namespace ProcessLib::Python
