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
#include "ProcessLib/LocalAssemblerInterface.h"
#include "ProcessLib/LocalAssemblerTraits.h"
#include "ProcessLib/Parameter.h"
#include "ProcessLib/Utils/InitShapeMatrices.h"
#include "GroundwaterFlowProcessData.h"

namespace ProcessLib
{

namespace GroundwaterFlow
{

enum class IntegrationPointValue {
    DarcyVelocityX,
    DarcyVelocityY,
    DarcyVelocityZ
};

const unsigned NUM_NODAL_DOF = 1;

template <typename GlobalMatrix, typename GlobalVector>
class GroundwaterFlowLocalAssemblerInterface
        : public ProcessLib::LocalAssemblerInterface<GlobalMatrix, GlobalVector>
        , public NumLib::Extrapolatable<GlobalVector, IntegrationPointValue>
{};

template <typename ShapeFunction,
         typename IntegrationMethod,
         typename GlobalMatrix,
         typename GlobalVector,
         unsigned GlobalDim>
class LocalAssemblerData
        : public GroundwaterFlowLocalAssemblerInterface<GlobalMatrix, GlobalVector>
{
    using ShapeMatricesType = ShapeMatrixPolicyType<ShapeFunction, GlobalDim>;
    using ShapeMatrices = typename ShapeMatricesType::ShapeMatrices;

    using LocalAssemblerTraits = ProcessLib::LocalAssemblerTraits<
        ShapeMatricesType, ShapeFunction::NPOINTS, NUM_NODAL_DOF, GlobalDim>;

    using NodalMatrixType = typename LocalAssemblerTraits::LocalMatrix;
    using NodalVectorType = typename LocalAssemblerTraits::LocalVector;

public:
    /// The hydraulic_conductivity factor is directly integrated into the local
    /// element matrix.
    LocalAssemblerData(MeshLib::Element const& element,
                       std::size_t const local_matrix_size,
                       unsigned const integration_order,
                       GroundwaterFlowProcessData const& process_data)
        : _element(element)
        , _shape_matrices(
              initShapeMatrices<ShapeFunction, ShapeMatricesType, IntegrationMethod, GlobalDim>(
                  element, integration_order))
        , _process_data(process_data)
        , _localA(local_matrix_size, local_matrix_size) // TODO narrowing conversion
        , _localRhs(local_matrix_size)
        , _integration_order(integration_order)
    {
        // This assertion is valid only if all nodal d.o.f. use the same shape matrices.
        assert(local_matrix_size == ShapeFunction::NPOINTS * NUM_NODAL_DOF);
    }

    void assemble(double const /*t*/, std::vector<double> const& local_x) override
    {
        _localA.setZero();
        _localRhs.setZero();

        IntegrationMethod integration_method(_integration_order);
        unsigned const n_integration_points = integration_method.getNPoints();

        for (std::size_t ip(0); ip < n_integration_points; ip++)
        {
            auto const& sm = _shape_matrices[ip];
            auto const& wp = integration_method.getWeightedPoint(ip);
            auto const k = _process_data.hydraulic_conductivity(_element);

            _localA.noalias() += sm.dNdx.transpose() * k * sm.dNdx *
                                 sm.detJ * wp.getWeight();

            // Darcy velocity only computed for output.
            auto const darcy_velocity = (k * sm.dNdx *
                Eigen::Map<const NodalVectorType>(local_x.data(), ShapeFunction::NPOINTS)
                ).eval();

            for (unsigned d=0; d<GlobalDim; ++d) {
                _darcy_velocities[d][ip] = darcy_velocity[d];
            }
        }
    }

    void addToGlobal(AssemblerLib::LocalToGlobalIndexMap::RowColumnIndices const& indices,
        GlobalMatrix& /*M*/, GlobalMatrix& K, GlobalVector& b)
        const override
    {
        K.add(indices, _localA);
        b.add(indices.rows, _localRhs);
    }

    Eigen::Map<const Eigen::VectorXd>
    getShapeMatrix(const unsigned integration_point) const override
    {
        auto const& N = _shape_matrices[integration_point].N;

        // assumes N is stored contiguously in memory
        return Eigen::Map<const Eigen::VectorXd>(N.data(), N.size());
    }

    std::vector<double> const&
    getIntegrationPointValues(IntegrationPointValue const property,
                              std::vector<double>& /*cache*/) const override
    {
        switch (property)
        {
        case IntegrationPointValue::DarcyVelocityX:
            return _darcy_velocities[0];
        case IntegrationPointValue::DarcyVelocityY:
            assert(GlobalDim > 1);
            return _darcy_velocities[1];
        case IntegrationPointValue::DarcyVelocityZ:
            assert(GlobalDim > 2);
            return _darcy_velocities[2];
        }

        std::abort();
    }

private:
    MeshLib::Element const& _element;
	std::vector<ShapeMatrices, Eigen::aligned_allocator<ShapeMatrices>> _shape_matrices;
    GroundwaterFlowProcessData const& _process_data;

    NodalMatrixType _localA;
    NodalVectorType _localRhs;

    unsigned const _integration_order;

    std::vector<std::vector<double>> _darcy_velocities
        = std::vector<std::vector<double>>(
            GlobalDim, std::vector<double>(ShapeFunction::NPOINTS));
};


}   // namespace GroundwaterFlow
}   // namespace ProcessLib

#endif  // PROCESS_LIB_GROUNDWATERFLOW_FEM_H_
