/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_GROUNDWATERFLOW_FEM_H_
#define PROCESS_LIB_GROUNDWATERFLOW_FEM_H_

#include <vector>

#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"

#include "Parameter.h"
#include "ProcessUtil.h"


namespace ProcessLib
{

namespace GroundwaterFlow
{

// TODO now this interface is basically the same for all processes that assemble a
//      FirstOrderImplicitQuasiLinear ODE system.
template <typename GlobalMatrix, typename GlobalVector>
class LocalAssemblerDataInterface
{
public:
    virtual ~LocalAssemblerDataInterface() = default;

    virtual void assemble(double const t, std::vector<double> const& local_x) = 0;

    virtual void addToGlobal(AssemblerLib::LocalToGlobalIndexMap::RowColumnIndices const&,
            GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b) const = 0;
};

template <typename ShapeFunction_,
         typename IntegrationMethod_,
         typename GlobalMatrix,
         typename GlobalVector,
         unsigned GlobalDim>
class LocalAssemblerData : public LocalAssemblerDataInterface<GlobalMatrix, GlobalVector>
{
public:
    using ShapeFunction = ShapeFunction_;
    using ShapeMatricesType = ShapeMatrixPolicyType<ShapeFunction, GlobalDim>;
    using NodalMatrixType = typename ShapeMatricesType::NodalMatrixType;
    using NodalVectorType = typename ShapeMatricesType::NodalVectorType;
    using ShapeMatrices = typename ShapeMatricesType::ShapeMatrices;

    /// The hydraulic_conductivity factor is directly integrated into the local
    /// element matrix.
    LocalAssemblerData(MeshLib::Element const& e,
                       std::size_t const local_matrix_size,
                       unsigned const integration_order,
                       Parameter<double, MeshLib::Element const&> const&
                       hydraulic_conductivity)
        : _shape_matrices(
              initShapeMatrices<ShapeFunction, ShapeMatricesType, IntegrationMethod_, GlobalDim>(
                  e, integration_order))
        , _hydraulic_conductivity([&hydraulic_conductivity, &e]()
          {
              return hydraulic_conductivity(e);
          })
        // TODO narrowing conversion
        , _localA(local_matrix_size, local_matrix_size)
        , _localRhs(local_matrix_size)
        , _integration_order(integration_order)
    {}

    void assemble(double const /*t*/, std::vector<double> const& /*local_x*/) override
    {
        _localA.setZero();
        _localRhs.setZero();

        IntegrationMethod_ integration_method(_integration_order);
        unsigned const n_integration_points = integration_method.getNPoints();

        for (std::size_t ip(0); ip < n_integration_points; ip++) {
            auto const& sm = _shape_matrices[ip];
            auto const& wp = integration_method.getWeightedPoint(ip);
            _localA.noalias() += sm.dNdx.transpose() *
                                 _hydraulic_conductivity() * sm.dNdx *
                                 sm.detJ * wp.getWeight();
        }
    }

    void addToGlobal(AssemblerLib::LocalToGlobalIndexMap::RowColumnIndices const& indices,
        GlobalMatrix& /*M*/, GlobalMatrix& K, GlobalVector& b)
        const override
    {
        K.add(indices, _localA);
        b.add(indices.rows, _localRhs);
    }

private:
    std::vector<ShapeMatrices> _shape_matrices;
    std::function<double(void)> _hydraulic_conductivity;

    NodalMatrixType _localA;
    NodalVectorType _localRhs;

    unsigned const _integration_order;
};


}   // namespace GroundwaterFlow
}   // namespace ProcessLib

#endif  // PROCESS_LIB_GROUNDWATERFLOW_FEM_H_
