/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <vector>

#include "HeatConductionProcessData.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "NumLib/Extrapolation/ExtrapolatableElement.h"
#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "ParameterLib/Parameter.h"
#include "ProcessLib/LocalAssemblerInterface.h"
#include "ProcessLib/LocalAssemblerTraits.h"
#include "ProcessLib/Utils/InitShapeMatrices.h"

namespace ProcessLib
{
namespace HeatConduction
{
const unsigned NUM_NODAL_DOF = 1;

class HeatConductionLocalAssemblerInterface
    : public ProcessLib::LocalAssemblerInterface,
      public NumLib::ExtrapolatableElement
{
public:
    virtual std::vector<double> const& getIntPtHeatFluxX(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtHeatFluxY(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtHeatFluxZ(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const = 0;
};

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
class LocalAssemblerData : public HeatConductionLocalAssemblerInterface
{
    using ShapeMatricesType = ShapeMatrixPolicyType<ShapeFunction, GlobalDim>;
    using ShapeMatrices = typename ShapeMatricesType::ShapeMatrices;

    using LocalAssemblerTraits = ProcessLib::LocalAssemblerTraits<
        ShapeMatricesType, ShapeFunction::NPOINTS, NUM_NODAL_DOF, GlobalDim>;

    using NodalMatrixType = typename LocalAssemblerTraits::LocalMatrix;
    using NodalVectorType = typename LocalAssemblerTraits::LocalVector;
    using GlobalDimVectorType = typename ShapeMatricesType::GlobalDimVectorType;

public:
    /// The thermal_conductivity factor is directly integrated into the local
    /// element matrix.
    LocalAssemblerData(MeshLib::Element const& element,
                       std::size_t const local_matrix_size,
                       bool is_axially_symmetric,
                       unsigned const integration_order,
                       HeatConductionProcessData const& process_data)
        : element_(element),
          process_data_(process_data),
          integration_method_(integration_order),
          shape_matrices_(initShapeMatrices<ShapeFunction, ShapeMatricesType,
                                            IntegrationMethod, GlobalDim>(
              element, is_axially_symmetric, integration_method_)),
          heat_fluxes_(
              GlobalDim,
              std::vector<double>(integration_method_.getNumberOfPoints()))
    {
        // This assertion is valid only if all nodal d.o.f. use the same shape
        // matrices.
        assert(local_matrix_size == ShapeFunction::NPOINTS * NUM_NODAL_DOF);
        (void)local_matrix_size;
    }

    void assemble(double const t, double const /*dt*/,
                  std::vector<double> const& local_x,
                  std::vector<double> const& /*local_xdot*/,
                  std::vector<double>& local_M_data,
                  std::vector<double>& local_K_data,
                  std::vector<double>& /*local_b_data*/) override
    {
        auto const local_matrix_size = local_x.size();
        // This assertion is valid only if all nodal d.o.f. use the same shape
        // matrices.
        assert(local_matrix_size == ShapeFunction::NPOINTS * NUM_NODAL_DOF);

        auto local_M = MathLib::createZeroedMatrix<NodalMatrixType>(
            local_M_data, local_matrix_size, local_matrix_size);
        auto local_K = MathLib::createZeroedMatrix<NodalMatrixType>(
            local_K_data, local_matrix_size, local_matrix_size);

        unsigned const n_integration_points =
            integration_method_.getNumberOfPoints();

        ParameterLib::SpatialPosition pos;
        pos.setElementID(element_.getID());

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            pos.setIntegrationPoint(ip);
            auto const& sm = shape_matrices_[ip];
            auto const& wp = integration_method_.getWeightedPoint(ip);
            auto const k = process_data_.thermal_conductivity(t, pos)[0];
            auto const heat_capacity = process_data_.heat_capacity(t, pos)[0];
            auto const density = process_data_.density(t, pos)[0];

            local_K.noalias() += sm.dNdx.transpose() * k * sm.dNdx * sm.detJ *
                                 wp.getWeight() * sm.integralMeasure;
            local_M.noalias() += sm.N.transpose() * density * heat_capacity *
                                 sm.N * sm.detJ * wp.getWeight() *
                                 sm.integralMeasure;
        }
    }

    void computeSecondaryVariableConcrete(
        double const t, double const /*dt*/, std::vector<double> const& local_x,
        std::vector<double> const& /*local_x_dot*/) override
    {
        auto const local_matrix_size = local_x.size();
        // This assertion is valid only if all nodal d.o.f. use the same shape
        // matrices.
        assert(local_matrix_size == ShapeFunction::NPOINTS * NUM_NODAL_DOF);

        unsigned const n_integration_points =
            integration_method_.getNumberOfPoints();

        ParameterLib::SpatialPosition pos;
        pos.setElementID(element_.getID());
        const auto local_x_vec =
            MathLib::toVector<NodalVectorType>(local_x, local_matrix_size);

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            pos.setIntegrationPoint(ip);
            auto const& sm = shape_matrices_[ip];
            auto const k = process_data_.thermal_conductivity(t, pos)[0];
            // heat flux only computed for output.
            GlobalDimVectorType const heat_flux = -k * sm.dNdx * local_x_vec;

            for (unsigned d = 0; d < GlobalDim; ++d)
            {
                heat_fluxes_[d][ip] = heat_flux[d];
            }
        }
    }

    Eigen::Map<const Eigen::RowVectorXd> getShapeMatrix(
        const unsigned integration_point) const override
    {
        auto const& N = shape_matrices_[integration_point].N;

        // assumes N is stored contiguously in memory
        return Eigen::Map<const Eigen::RowVectorXd>(N.data(), N.size());
    }

    std::vector<double> const& getIntPtHeatFluxX(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& /*cache*/) const override
    {
        assert(!heat_fluxes_.empty());
        return heat_fluxes_[0];
    }

    std::vector<double> const& getIntPtHeatFluxY(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& /*cache*/) const override
    {
        assert(heat_fluxes_.size() > 1);
        return heat_fluxes_[1];
    }

    std::vector<double> const& getIntPtHeatFluxZ(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& /*cache*/) const override
    {
        assert(heat_fluxes_.size() > 2);
        return heat_fluxes_[2];
    }

private:
    MeshLib::Element const& element_;
    HeatConductionProcessData const& process_data_;

    IntegrationMethod const integration_method_;
    std::vector<ShapeMatrices, Eigen::aligned_allocator<ShapeMatrices>>
        shape_matrices_;

    std::vector<std::vector<double>> heat_fluxes_;
};

}  // namespace HeatConduction
}  // namespace ProcessLib
