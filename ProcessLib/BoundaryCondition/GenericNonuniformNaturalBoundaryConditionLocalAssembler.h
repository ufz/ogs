/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "ProcessLib/Utils/InitShapeMatrices.h"

namespace ProcessLib
{
class GenericNonuniformNaturalBoundaryConditionLocalAssemblerInterface
{
public:
    virtual ~GenericNonuniformNaturalBoundaryConditionLocalAssemblerInterface() =
        default;

    virtual void assemble(
        std::size_t const id,
        NumLib::LocalToGlobalIndexMap const& dof_table_boundary, double const t,
        const GlobalVector& x, GlobalMatrix& K, GlobalVector& b,
        GlobalMatrix* Jac) = 0;
};

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
class GenericNonuniformNaturalBoundaryConditionLocalAssembler
    : public GenericNonuniformNaturalBoundaryConditionLocalAssemblerInterface
{
protected:
    using ShapeMatricesType = ShapeMatrixPolicyType<ShapeFunction, GlobalDim>;
    using NodalMatrixType = typename ShapeMatricesType::NodalMatrixType;
    using NodalVectorType = typename ShapeMatricesType::NodalVectorType;

    struct NAndWeight
    {
        NAndWeight(typename ShapeMatricesType::ShapeMatrices::ShapeType&& N_,
                   double const weight_)
            : N(std::move(N_)), weight(weight_)
        {
        }
        typename ShapeMatricesType::ShapeMatrices::ShapeType const N;
        double const weight;
    };

private:
    static std::vector<NAndWeight, Eigen::aligned_allocator<NAndWeight>>
    initNsAndWeights(MeshLib::Element const& e, bool is_axially_symmetric,
                     unsigned const integration_order)
    {
        IntegrationMethod const integration_method(integration_order);
        std::vector<NAndWeight, Eigen::aligned_allocator<NAndWeight>>
            ns_and_weights;
        ns_and_weights.reserve(integration_method.getNumberOfPoints());

        auto sms = initShapeMatrices<ShapeFunction, ShapeMatricesType,
                                     IntegrationMethod, GlobalDim>(
            e, is_axially_symmetric, integration_method);
        for (unsigned ip = 0; ip < sms.size(); ++ip)
        {
            auto& sm = sms[ip];
            double const w =
                sm.detJ * sm.integralMeasure *
                integration_method.getWeightedPoint(ip).getWeight();

            ns_and_weights.emplace_back(std::move(sm.N), w);
        }

        return ns_and_weights;
    }

public:
    GenericNonuniformNaturalBoundaryConditionLocalAssembler(
        MeshLib::Element const& e, bool is_axially_symmetric,
        unsigned const integration_order)
        : _ns_and_weights(
              initNsAndWeights(e, is_axially_symmetric, integration_order))
    {
    }

protected:
    std::vector<NAndWeight, Eigen::aligned_allocator<NAndWeight>> const
        _ns_and_weights;
};

}  // namespace ProcessLib
