/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_RICHARDSFLOW_FEM_H_
#define PROCESS_LIB_RICHARDSFLOW_FEM_H_

#include <vector>

#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "MeshLib/CoordinateSystem.h"
#include "NumLib/Extrapolation/ExtrapolatableElement.h"
#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "NumLib/Function/Interpolation.h"
#include "ProcessLib/LocalAssemblerInterface.h"
#include "ProcessLib/LocalAssemblerTraits.h"
#include "ProcessLib/Parameter/Parameter.h"
#include "ProcessLib/Utils/InitShapeMatrices.h"
#include "RichardsFlowProcessData.h"

namespace ProcessLib
{
namespace RichardsFlow
{
const unsigned NUM_NODAL_DOF = 1;

class RichardsFlowLocalAssemblerInterface
    : public ProcessLib::LocalAssemblerInterface,
      public NumLib::ExtrapolatableElement
{
public:
    virtual std::vector<double> const& getIntPtSaturation(
        std::vector<double>& /*cache*/) const = 0;
};

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
class LocalAssemblerData : public RichardsFlowLocalAssemblerInterface
{
    using ShapeMatricesType = ShapeMatrixPolicyType<ShapeFunction, GlobalDim>;
    using ShapeMatrices = typename ShapeMatricesType::ShapeMatrices;

    using LocalAssemblerTraits = ProcessLib::LocalAssemblerTraits<
        ShapeMatricesType, ShapeFunction::NPOINTS, NUM_NODAL_DOF, GlobalDim>;

    using NodalMatrixType = typename LocalAssemblerTraits::LocalMatrix;
    using NodalVectorType = typename LocalAssemblerTraits::LocalVector;
    using GlobalDimVectorType = typename ShapeMatricesType::GlobalDimVectorType;

public:
    LocalAssemblerData(MeshLib::Element const& element,
                       std::size_t const local_matrix_size,
                       bool is_axially_symmetric,
                       unsigned const integration_order,
                       RichardsFlowProcessData const& process_data)
        : _element(element),
          _process_data(process_data),
          _integration_method(integration_order),
          _shape_matrices(initShapeMatrices<ShapeFunction, ShapeMatricesType,
                                            IntegrationMethod, GlobalDim>(
              element, is_axially_symmetric, _integration_method)),
          _saturation(
              std::vector<double>(_integration_method.getNumberOfPoints()))
    {
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

        MathLib::PiecewiseLinearInterpolation const& interpolated_Pc =
            *_process_data.curves.at("curve_PC_S");
        MathLib::PiecewiseLinearInterpolation const& interpolated_Kr =
            *_process_data.curves.at("curve_S_Krel");

        auto local_M = MathLib::createZeroedMatrix<NodalMatrixType>(
            local_M_data, local_matrix_size, local_matrix_size);
        auto local_K = MathLib::createZeroedMatrix<NodalMatrixType>(
            local_K_data, local_matrix_size, local_matrix_size);
        auto local_b = MathLib::createZeroedVector<NodalVectorType>(
            local_b_data, local_matrix_size);

        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();
        SpatialPosition pos;
        pos.setElementID(_element.getID());
        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            pos.setIntegrationPoint(ip);
            auto const& sm = _shape_matrices[ip];
            auto const& wp = _integration_method.getWeightedPoint(ip);

            double P_int_pt = 0.0;
            NumLib::shapeFunctionInterpolate(local_x, sm.N, P_int_pt);

            auto const K = _process_data.intrinsic_permeability(t, pos)[0];
            auto const poro = _process_data.porosity(t, pos)[0];
            auto const mu = _process_data.viscosity(t, pos)[0];
            auto const storage = _process_data.storage(t, pos)[0];
            double const Pc = -P_int_pt;

            double dSwdPc = interpolated_Pc.getDerivative(Pc);
            if (Pc > interpolated_Pc.getSupportMax())
                dSwdPc = interpolated_Pc.getDerivative(
                    interpolated_Pc.getSupportMax());
            else if (Pc < interpolated_Pc.getSupportMin())
                dSwdPc = interpolated_Pc.getDerivative(
                    interpolated_Pc.getSupportMin());

            double const Sw = interpolated_Pc.getValue(Pc);
            _saturation[ip] = Sw;

            double const k_rel = interpolated_Kr.getValue(Sw);
            double const drhow_dp(0.0);

            double const mass_mat_coeff =
                storage * Sw + poro * Sw * drhow_dp - poro * dSwdPc;
            double const K_mat_coeff = K * k_rel / mu;

            local_K.noalias() += sm.dNdx.transpose() * K_mat_coeff * sm.dNdx *
                                 sm.detJ * sm.integralMeasure * wp.getWeight();

            local_M.noalias() += sm.N.transpose() * mass_mat_coeff * sm.N *
                                 sm.detJ * sm.integralMeasure * wp.getWeight();
            if (_process_data.has_gravity)
            {
                auto const rho_w = _process_data.water_density(t, pos)[0];
                auto const body_force =
                    _process_data.specific_body_force(t, pos);
                assert(body_force.size() == GlobalDim);
                auto const b = MathLib::toVector<GlobalDimVectorType>(
                    body_force, GlobalDim);
                local_b.noalias() += sm.dNdx.transpose() * K_mat_coeff * rho_w *
                                     b * sm.detJ * sm.integralMeasure *
                                     wp.getWeight();
            }
        }  // end of GP
        if (_process_data.has_mass_lumping)
        {
            for (int idx_ml = 0; idx_ml < local_M.cols(); idx_ml++)
            {
                double const mass_lump_val = local_M.col(idx_ml).sum();
                local_M.col(idx_ml).setZero();
                local_M(idx_ml, idx_ml) = mass_lump_val;
            }
        }  // end of mass lumping
    }

    Eigen::Map<const Eigen::RowVectorXd> getShapeMatrix(
        const unsigned integration_point) const override
    {
        auto const& N = _shape_matrices[integration_point].N;

        // assumes N is stored contiguously in memory
        return Eigen::Map<const Eigen::RowVectorXd>(N.data(), N.size());
    }

    std::vector<double> const& getIntPtSaturation(
        std::vector<double>& /*cache*/) const override
    {
        assert(_saturation.size() > 0);
        return _saturation;
    }

private:
    MeshLib::Element const& _element;
    RichardsFlowProcessData const& _process_data;

    IntegrationMethod const _integration_method;
    std::vector<ShapeMatrices> _shape_matrices;

    std::vector<double> _saturation;
};

}  // namespace RichardsFlow
}  // namespace ProcessLib

#endif  // PROCESS_LIB_RichardsFlow_FEM_H_
