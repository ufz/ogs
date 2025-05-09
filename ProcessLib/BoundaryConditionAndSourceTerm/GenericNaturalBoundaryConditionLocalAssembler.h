/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "MeshLib/Elements/Element.h"
#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "NumLib/Fem/InitShapeMatrices.h"
#include "NumLib/Fem/Integration/GenericIntegrationMethod.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"

namespace ProcessLib
{
class GenericNaturalBoundaryConditionLocalAssemblerInterface
{
public:
    virtual ~GenericNaturalBoundaryConditionLocalAssemblerInterface() = default;

    virtual void assemble(
        std::size_t const id,
        NumLib::LocalToGlobalIndexMap const& dof_table_boundary, double const t,
        std::vector<GlobalVector*> const& x, int const process_id,
        GlobalMatrix* K, GlobalVector& b, GlobalMatrix* Jac) = 0;
};

template <typename ShapeFunction, int GlobalDim>
class GenericNaturalBoundaryConditionLocalAssembler
    : public GenericNaturalBoundaryConditionLocalAssemblerInterface
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
                     NumLib::GenericIntegrationMethod const& integration_method)
    {
        std::vector<NAndWeight, Eigen::aligned_allocator<NAndWeight>>
            ns_and_weights;
        ns_and_weights.reserve(integration_method.getNumberOfPoints());

        auto sms = NumLib::initShapeMatrices<ShapeFunction, ShapeMatricesType,
                                             GlobalDim>(e, is_axially_symmetric,
                                                        integration_method);
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
    GenericNaturalBoundaryConditionLocalAssembler(
        MeshLib::Element const& e, bool is_axially_symmetric,
        NumLib::GenericIntegrationMethod const& integration_method)
        : _integration_method(integration_method),
          _ns_and_weights(
              initNsAndWeights(e, is_axially_symmetric, integration_method)),
          _element(e)
    {
    }

protected:
    NumLib::GenericIntegrationMethod const& _integration_method;
    std::vector<NAndWeight, Eigen::aligned_allocator<NAndWeight>> const
        _ns_and_weights;
    MeshLib::Element const& _element;
};

}  // namespace ProcessLib
