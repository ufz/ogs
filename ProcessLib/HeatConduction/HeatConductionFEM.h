/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
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
#include "ProcessLib/LocalAssemblerInterface.h"
#include "ProcessLib/LocalAssemblerTraits.h"
#include "ProcessLib/Parameter/Parameter.h"
#include "ProcessLib/Utils/InitShapeMatrices.h"

// For coupling
#include "ProcessLib/LiquidFlow/LiquidFlowProcess.h"
#include "ProcessLib/LiquidFlow/LiquidFlowMaterialProperties.h"

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
        std::vector<double>& /*cache*/) const = 0;

    virtual std::vector<double> const& getIntPtHeatFluxY(
        std::vector<double>& /*cache*/) const = 0;

    virtual std::vector<double> const& getIntPtHeatFluxZ(
        std::vector<double>& /*cache*/) const = 0;
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
        : _element(element),
          _process_data(process_data),
          _integration_method(integration_order),
          _shape_matrices(initShapeMatrices<ShapeFunction, ShapeMatricesType,
                                            IntegrationMethod, GlobalDim>(
              element, is_axially_symmetric, _integration_method)),
          _heat_fluxes(
              GlobalDim,
              std::vector<double>(_integration_method.getNumberOfPoints()))
    {
        // This assertion is valid only if all nodal d.o.f. use the same shape
        // matrices.
        assert(local_matrix_size == ShapeFunction::NPOINTS * NUM_NODAL_DOF);
        (void)local_matrix_size;
    }

    void assemble(double const t, std::vector<double> const& local_x,
                  std::vector<double>& local_M_data,
                  std::vector<double>& local_K_data,
                  std::vector<double>& local_b_data) override
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
            _integration_method.getNumberOfPoints();

        SpatialPosition pos;
        pos.setElementID(_element.getID());

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            pos.setIntegrationPoint(ip);
            auto const& sm = _shape_matrices[ip];
            auto const& wp = _integration_method.getWeightedPoint(ip);
            auto const k = _process_data.thermal_conductivity(t, pos)[0];
            auto const heat_capacity = _process_data.heat_capacity(t, pos)[0];
            auto const density = _process_data.density(t, pos)[0];

            local_K.noalias() += sm.dNdx.transpose() * k * sm.dNdx * sm.detJ *
                                 wp.getWeight() * sm.integralMeasure;
            local_M.noalias() += sm.N.transpose() * density * heat_capacity *
                                 sm.N * sm.detJ * wp.getWeight() *
                                 sm.integralMeasure;
        }
    }

    void assembleHeatTransportLiquidFlow(
        double const t, int const gravitational_axis_id,
        ProcessLib::LiquidFlow::LiquidFlowMaterialProperties const&
            liquid_flow_prop,
        std::vector<double> const& local_x, std::vector<double> const& local_p,
        std::vector<double>& local_M_data, std::vector<double>& local_K_data,
        std::vector<double>& local_b_data)
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
            _integration_method.getNumberOfPoints();

        SpatialPosition pos;
        pos.setElementID(_element.getID());
        const int material_id = liquid_flow_prop.getMaterialID(pos);

        const Eigen::MatrixXd& perm = liquid_flow_prop.getPermeability(
            material_id, t, pos, _element.getDimension());

        const double porosity_variable = 0.;
        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            pos.setIntegrationPoint(ip);
            auto const& sm = _shape_matrices[ip];
            auto const& wp = _integration_method.getWeightedPoint(ip);
            double p = 0.;
            NumLib::shapeFunctionInterpolate(local_p, sm.N, p);
            double T = 0.;
            NumLib::shapeFunctionInterpolate(local_x, sm.N, T);

            // Material parameters of solid phase
            auto const k_s = _process_data.thermal_conductivity(t, pos)[0];
            auto const cp_s = _process_data.heat_capacity(t, pos)[0];
            auto const rho_s = _process_data.density(t, pos)[0];

            // Material parameters of fluid phase
            double const cp_f = liquid_flow_prop.getHeatCapacity(p, T);
            double const k_f = liquid_flow_prop.getThermalConductivity(p, T);
            double rho_f = liquid_flow_prop.getLiquidDensity(p, T);

            // Material parameter of porosity
            double const poro =
                liquid_flow_prop.getPorosity(material_id, porosity_variable, T);

            double const effective_cp =
                (1.0 - poro) * cp_s * rho_s + poro * cp_f * rho_f;
            double const effective_K = (1.0 - poro) * k_s + poro * k_f;

            double const integration_factor =
                sm.detJ * wp.getWeight() * sm.integralMeasure;

            local_M.noalias() +=
                effective_cp * sm.N.transpose() * sm.N * integration_factor;
            local_K.noalias() += sm.dNdx.transpose() * effective_K * sm.dNdx *
                                 integration_factor;
        }
    }

    void coupling_assemble(double const t, std::vector<double> const& local_x,
                           std::vector<double>& local_M_data,
                           std::vector<double>& local_K_data,
                           std::vector<double>& local_b_data,
                           LocalCouplingTerm const& coupled_term) override
    {
        auto it = coupled_term.coupled_processes.begin();
        while (it != coupled_term.coupled_processes.end())
        {
            switch (it->first)
            {
                case ProcessLib::ProcessType::LiquidFlowProcess:
                {
                    ProcessLib::LiquidFlow::LiquidFlowProcess const& pcs =
                        static_cast<
                            ProcessLib::LiquidFlow::LiquidFlowProcess const&>(
                            it->second);
                    const auto liquid_flow_prop =
                        pcs.getLiquidFlowMaterialProperties();

                    const auto local_p = coupled_term.local_coupled_xs.at(
                        ProcessLib::ProcessType::LiquidFlowProcess);

                    int const gravitational_axis_id =
                        pcs.getGravitationalAxisID();

                    assembleHeatTransportLiquidFlow(
                        t, gravitational_axis_id, *liquid_flow_prop, local_x,
                        local_p, local_M_data, local_K_data, local_b_data);
                }
                break;
                default:
                    OGS_FATAL(
                        "This coupled process is not presented for "
                        "HeatConduction process");
            }
            it++;
        }
    }

    void computeSecondaryVariableConcrete(
        const double t, std::vector<double> const& local_x) override
    {
        auto const local_matrix_size = local_x.size();
        // This assertion is valid only if all nodal d.o.f. use the same shape
        // matrices.
        assert(local_matrix_size == ShapeFunction::NPOINTS * NUM_NODAL_DOF);

        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        SpatialPosition pos;
        pos.setElementID(_element.getID());
        const auto local_x_vec =
            MathLib::toVector<NodalVectorType>(local_x, local_matrix_size);

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            pos.setIntegrationPoint(ip);
            auto const& sm = _shape_matrices[ip];
            auto const k = _process_data.thermal_conductivity(t, pos)[0];
            // heat flux only computed for output.
            GlobalDimVectorType const heat_flux = -k * sm.dNdx * local_x_vec;

            for (unsigned d = 0; d < GlobalDim; ++d)
            {
                _heat_fluxes[d][ip] = heat_flux[d];
            }
        }
    }

    Eigen::Map<const Eigen::RowVectorXd> getShapeMatrix(
        const unsigned integration_point) const override
    {
        auto const& N = _shape_matrices[integration_point].N;

        // assumes N is stored contiguously in memory
        return Eigen::Map<const Eigen::RowVectorXd>(N.data(), N.size());
    }

    std::vector<double> const& getIntPtHeatFluxX(
        std::vector<double>& /*cache*/) const override
    {
        assert(_heat_fluxes.size() > 0);
        return _heat_fluxes[0];
    }

    std::vector<double> const& getIntPtHeatFluxY(
        std::vector<double>& /*cache*/) const override
    {
        assert(_heat_fluxes.size() > 1);
        return _heat_fluxes[1];
    }

    std::vector<double> const& getIntPtHeatFluxZ(
        std::vector<double>& /*cache*/) const override
    {
        assert(_heat_fluxes.size() > 2);
        return _heat_fluxes[2];
    }

private:
    MeshLib::Element const& _element;
    HeatConductionProcessData const& _process_data;

    IntegrationMethod const _integration_method;
    std::vector<ShapeMatrices, Eigen::aligned_allocator<ShapeMatrices>>
        _shape_matrices;

    std::vector<std::vector<double>> _heat_fluxes;
};

}  // namespace HeatConduction
}  // namespace ProcessLib
