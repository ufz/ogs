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

#include <vector>

#include "MaterialLib/PhysicalConstant.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "NumLib/Extrapolation/ExtrapolatableElement.h"
#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/InitShapeMatrices.h"
#include "NumLib/Fem/Integration/GenericIntegrationMethod.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "ParameterLib/Parameter.h"
#include "ProcessLib/LocalAssemblerInterface.h"
#include "ProcessLib/LocalAssemblerTraits.h"
#include "ThermalTwoPhaseFlowWithPPProcessData.h"

namespace ProcessLib
{
namespace ThermalTwoPhaseFlowWithPP
{
template <typename NodalRowVectorType, typename GlobalDimNodalMatrixType,
          typename NodalMatrixType>
struct IntegrationPointData final
{
    explicit IntegrationPointData(NodalRowVectorType const& N_,
                                  GlobalDimNodalMatrixType const& dNdx_,
                                  double const& integration_weight_,
                                  NodalMatrixType const mass_operator_,
                                  NodalMatrixType const diffusion_operator_)
        : N(N_),
          dNdx(dNdx_),
          integration_weight(integration_weight_),
          mass_operator(mass_operator_),
          diffusion_operator(diffusion_operator_)
    {
    }
    NodalRowVectorType const N;
    GlobalDimNodalMatrixType const dNdx;
    double const integration_weight;
    NodalMatrixType const mass_operator;
    NodalMatrixType const diffusion_operator;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};
const unsigned NUM_NODAL_DOF = 4;

class ThermalTwoPhaseFlowWithPPLocalAssemblerInterface
    : public ProcessLib::LocalAssemblerInterface,
      public NumLib::ExtrapolatableElement
{
public:
    virtual std::vector<double> const& getIntPtSaturation(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtWettingPressure(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtLiquidMolFracContaminant(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtGasMolFracWater(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtGasMolFracContaminant(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const = 0;
};

template <typename ShapeFunction, int GlobalDim>
class ThermalTwoPhaseFlowWithPPLocalAssembler
    : public ThermalTwoPhaseFlowWithPPLocalAssemblerInterface
{
    using ShapeMatricesType = ShapeMatrixPolicyType<ShapeFunction, GlobalDim>;
    using ShapeMatrices = typename ShapeMatricesType::ShapeMatrices;

    using LocalAssemblerTraits = ProcessLib::LocalAssemblerTraits<
        ShapeMatricesType, ShapeFunction::NPOINTS, NUM_NODAL_DOF, GlobalDim>;
    using NodalRowVectorType = typename ShapeMatricesType::NodalRowVectorType;

    using GlobalDimNodalMatrixType =
        typename ShapeMatricesType::GlobalDimNodalMatrixType;
    using NodalMatrixType = typename ShapeMatricesType::NodalMatrixType;
    using NodalVectorType = typename ShapeMatricesType::NodalVectorType;
    using GlobalDimMatrixType = typename ShapeMatricesType::GlobalDimMatrixType;
    using GlobalDimVectorType = typename ShapeMatricesType::GlobalDimVectorType;
    using LocalMatrixType = typename LocalAssemblerTraits::LocalMatrix;
    using LocalVectorType = typename LocalAssemblerTraits::LocalVector;

public:
    ThermalTwoPhaseFlowWithPPLocalAssembler(
        MeshLib::Element const& element,
        std::size_t const /*local_matrix_size*/,
        NumLib::GenericIntegrationMethod const& integration_method,
        bool const is_axially_symmetric,
        ThermalTwoPhaseFlowWithPPProcessData const& process_data)
        : _element(element),
          _integration_method(integration_method),
          _process_data(process_data),
          _saturation(
              std::vector<double>(_integration_method.getNumberOfPoints())),
          _pressure_wetting(
              std::vector<double>(_integration_method.getNumberOfPoints())),
          _liquid_molar_fraction_contaminant(
              std::vector<double>(_integration_method.getNumberOfPoints())),
          _gas_molar_fraction_water(
              std::vector<double>(_integration_method.getNumberOfPoints())),
          _gas_molar_fraction_contaminant(
              std::vector<double>(_integration_method.getNumberOfPoints()))
    {
        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();
        _ip_data.reserve(n_integration_points);
        auto const shape_matrices =
            NumLib::initShapeMatrices<ShapeFunction, ShapeMatricesType,
                                      GlobalDim>(element, is_axially_symmetric,
                                                 _integration_method);
        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            auto const& sm = shape_matrices[ip];
            const double integration_factor = sm.integralMeasure * sm.detJ;
            _ip_data.emplace_back(
                sm.N, sm.dNdx,
                sm.integralMeasure * sm.detJ *
                    _integration_method.getWeightedPoint(ip).getWeight(),
                sm.N.transpose() * sm.N * integration_factor *
                    _integration_method.getWeightedPoint(ip).getWeight(),
                sm.dNdx.transpose() * sm.dNdx * integration_factor *
                    _integration_method.getWeightedPoint(ip).getWeight());
        }
    }

    void assemble(double const t, double const dt,
                  std::vector<double> const& local_x,
                  std::vector<double> const& local_x_prev,
                  std::vector<double>& local_M_data,
                  std::vector<double>& local_K_data,
                  std::vector<double>& local_b_data) override;

    Eigen::Map<const Eigen::RowVectorXd> getShapeMatrix(
        const unsigned integration_point) const override
    {
        auto const& N = _ip_data[integration_point].N;

        // assumes N is stored contiguously in memory
        return Eigen::Map<const Eigen::RowVectorXd>(N.data(), N.size());
    }

    std::vector<double> const& getIntPtSaturation(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& /*cache*/) const override
    {
        assert(!_saturation.empty());
        return _saturation;
    }

    std::vector<double> const& getIntPtWettingPressure(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& /*cache*/) const override
    {
        assert(!_pressure_wetting.empty());
        return _pressure_wetting;
    }

    std::vector<double> const& getIntPtLiquidMolFracContaminant(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& /*cache*/) const override
    {
        assert(!_liquid_molar_fraction_contaminant.empty());
        return _liquid_molar_fraction_contaminant;
    }

    std::vector<double> const& getIntPtGasMolFracWater(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& /*cache*/) const override
    {
        assert(!_gas_molar_fraction_water.empty());
        return _gas_molar_fraction_water;
    }

    std::vector<double> const& getIntPtGasMolFracContaminant(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& /*cache*/) const override
    {
        assert(!_gas_molar_fraction_contaminant.empty());
        return _gas_molar_fraction_contaminant;
    }

private:
    MeshLib::Element const& _element;

    NumLib::GenericIntegrationMethod const& _integration_method;

    ThermalTwoPhaseFlowWithPPProcessData const& _process_data;
    std::vector<
        IntegrationPointData<NodalRowVectorType, GlobalDimNodalMatrixType,
                             NodalMatrixType>,
        Eigen::aligned_allocator<IntegrationPointData<
            NodalRowVectorType, GlobalDimNodalMatrixType, NodalMatrixType>>>
        _ip_data;

    std::vector<double> _saturation;
    std::vector<double> _pressure_wetting;
    std::vector<double> _liquid_molar_fraction_contaminant;
    std::vector<double> _gas_molar_fraction_water;
    std::vector<double> _gas_molar_fraction_contaminant;

    static const int nonwet_pressure_matrix_index = 0;
    static const int cap_pressure_matrix_index = ShapeFunction::NPOINTS;
    static const int contaminant_matrix_index = 2 * ShapeFunction::NPOINTS;
    static const int temperature_matrix_index = 3 * ShapeFunction::NPOINTS;

    static const int nonwet_pressure_size = ShapeFunction::NPOINTS;
    static const int cap_pressure_size = ShapeFunction::NPOINTS;
    static const int contaminant_size = ShapeFunction::NPOINTS;
    static const int temperature_size = ShapeFunction::NPOINTS;
};

}  // namespace ThermalTwoPhaseFlowWithPP
}  // namespace ProcessLib

#include "ThermalTwoPhaseFlowWithPPLocalAssembler-impl.h"
