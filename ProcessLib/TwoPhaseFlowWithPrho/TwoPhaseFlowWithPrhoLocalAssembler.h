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
#include "MaterialLib/PhysicalConstant.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "NumLib/Extrapolation/ExtrapolatableElement.h"
#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "ParameterLib/Parameter.h"
#include "ProcessLib/LocalAssemblerInterface.h"
#include "ProcessLib/LocalAssemblerTraits.h"
#include "ProcessLib/Utils/InitShapeMatrices.h"

#include "TwoPhaseFlowWithPrhoProcessData.h"

namespace ProcessLib
{
namespace TwoPhaseFlowWithPrho
{
template <typename NodalMatrixType>
struct IntegrationPointData final
{
    explicit IntegrationPointData(
        TwoPhaseFlowWithPrhoMaterialProperties& material_property_)
        : mat_property(material_property_),
          sw(1.0),
          rho_m(0.0),
          dsw_dpg(0.0),
          dsw_drho(0.0),
          drhom_dpg(0.0),
          drhom_drho(0.0)
    {
    }
    TwoPhaseFlowWithPrhoMaterialProperties& mat_property;
    double sw;
    double rho_m;
    double dsw_dpg;
    double dsw_drho;
    double drhom_dpg;
    double drhom_drho;
    double pressure_nonwetting;

    double integration_weight;
    NodalMatrixType massOperator;
    NodalMatrixType diffusionOperator;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};
const unsigned NUM_NODAL_DOF = 2;

class TwoPhaseFlowWithPrhoLocalAssemblerInterface
    : public ProcessLib::LocalAssemblerInterface,
      public NumLib::ExtrapolatableElement
{
public:
    virtual std::vector<double> const& getIntPtSaturation(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtNonWettingPressure(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const = 0;
};

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
class TwoPhaseFlowWithPrhoLocalAssembler
    : public TwoPhaseFlowWithPrhoLocalAssemblerInterface
{
    using ShapeMatricesType = ShapeMatrixPolicyType<ShapeFunction, GlobalDim>;
    using ShapeMatrices = typename ShapeMatricesType::ShapeMatrices;

    using LocalAssemblerTraits = ProcessLib::LocalAssemblerTraits<
        ShapeMatricesType, ShapeFunction::NPOINTS, NUM_NODAL_DOF, GlobalDim>;

    using NodalMatrixType = typename ShapeMatricesType::NodalMatrixType;
    using NodalVectorType = typename ShapeMatricesType::NodalVectorType;
    using GlobalDimMatrixType = typename ShapeMatricesType::GlobalDimMatrixType;
    using GlobalDimVectorType = typename ShapeMatricesType::GlobalDimVectorType;
    using LocalMatrixType = typename LocalAssemblerTraits::LocalMatrix;
    using LocalVectorType = typename LocalAssemblerTraits::LocalVector;

public:
    TwoPhaseFlowWithPrhoLocalAssembler(
        MeshLib::Element const& element,
        std::size_t const /*local_matrix_size*/,
        bool const is_axially_symmetric,
        unsigned const integration_order,
        TwoPhaseFlowWithPrhoProcessData const& process_data)
        : element_(element),
          integration_method_(integration_order),
          shape_matrices_(initShapeMatrices<ShapeFunction, ShapeMatricesType,
                                            IntegrationMethod, GlobalDim>(
              element, is_axially_symmetric, integration_method_)),
          process_data_(process_data),
          saturation_(
              std::vector<double>(integration_method_.getNumberOfPoints())),
          pressure_nonwetting_(
              std::vector<double>(integration_method_.getNumberOfPoints()))
    {
        unsigned const n_integration_points =
            integration_method_.getNumberOfPoints();
        ip_data_.reserve(n_integration_points);
        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            ip_data_.emplace_back(*process_data_.material_);
            auto const& sm = shape_matrices_[ip];
            ip_data_[ip].integration_weight =
                sm.integralMeasure * sm.detJ *
                integration_method_.getWeightedPoint(ip).getWeight();
            ip_data_[ip].massOperator.setZero(ShapeFunction::NPOINTS,
                                              ShapeFunction::NPOINTS);
            ip_data_[ip].diffusionOperator.setZero(ShapeFunction::NPOINTS,
                                                   ShapeFunction::NPOINTS);
            ip_data_[ip].massOperator.noalias() =
                sm.N.transpose() * sm.N * ip_data_[ip].integration_weight;
            ip_data_[ip].diffusionOperator.noalias() =
                sm.dNdx.transpose() * sm.dNdx * ip_data_[ip].integration_weight;
        }
    }

    void assemble(double const t, double const /*dt*/,
                  std::vector<double> const& local_x,
                  std::vector<double> const& local_xdot,
                  std::vector<double>& local_M_data,
                  std::vector<double>& local_K_data,
                  std::vector<double>& local_b_data) override;

    Eigen::Map<const Eigen::RowVectorXd> getShapeMatrix(
        const unsigned integration_point) const override
    {
        auto const& N = shape_matrices_[integration_point].N;

        // assumes N is stored contiguously in memory
        return Eigen::Map<const Eigen::RowVectorXd>(N.data(), N.size());
    }

    std::vector<double> const& getIntPtSaturation(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& /*cache*/) const override
    {
        assert(!saturation_.empty());
        return saturation_;
    }

    std::vector<double> const& getIntPtNonWettingPressure(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& /*cache*/) const override
    {
        assert(!pressure_nonwetting_.empty());
        return pressure_nonwetting_;
    }

private:
    MeshLib::Element const& element_;

    IntegrationMethod const integration_method_;
    std::vector<ShapeMatrices, Eigen::aligned_allocator<ShapeMatrices>>
        shape_matrices_;

    TwoPhaseFlowWithPrhoProcessData const& process_data_;
    std::vector<IntegrationPointData<NodalMatrixType>,
                Eigen::aligned_allocator<IntegrationPointData<NodalMatrixType>>>
        ip_data_;

    std::vector<double> saturation_;  /// used for secondary variable output
    std::vector<double> pressure_nonwetting_;
    static const int nonwet_pressure_coeff_index = 0;
    static const int cap_pressure_coeff_index = 1;

    static const int nonwet_pressure_matrix_index = 0;
    static const int cap_pressure_matrix_index = ShapeFunction::NPOINTS;

    static const int nonwet_pressure_size = ShapeFunction::NPOINTS;
    static const int cap_pressure_size = ShapeFunction::NPOINTS;
};

}  // namespace TwoPhaseFlowWithPrho
}  // namespace ProcessLib

#include "TwoPhaseFlowWithPrhoLocalAssembler-impl.h"
