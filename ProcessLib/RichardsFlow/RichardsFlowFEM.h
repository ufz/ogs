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
template <typename NodalRowVectorType, typename GlobalDimNodalMatrixType,
          typename NodalMatrixType>
struct IntegrationPointData final
{
    IntegrationPointData(NodalRowVectorType const& N_,
                         GlobalDimNodalMatrixType const& dNdx_,
                         double const& integration_weight_,
                         NodalMatrixType const mass_operator_)
        : N(N_),
          dNdx(dNdx_),
          integration_weight(integration_weight_),
          mass_operator(mass_operator_)
    {
    }

    NodalRowVectorType const N;
    GlobalDimNodalMatrixType const dNdx;
    double const integration_weight;
    NodalMatrixType const mass_operator;
};
const unsigned NUM_NODAL_DOF = 1;

class RichardsFlowLocalAssemblerInterface
    : public ProcessLib::LocalAssemblerInterface,
      public NumLib::ExtrapolatableElement
{
public:
    virtual std::vector<double> const& getIntPtSaturation(
        std::vector<double>& /*cache*/) const = 0;

    virtual std::vector<double> const& getIntPtDarcyVelocityX(
        std::vector<double>& /*cache*/) const = 0;

    virtual std::vector<double> const& getIntPtDarcyVelocityY(
        std::vector<double>& /*cache*/) const = 0;

    virtual std::vector<double> const& getIntPtDarcyVelocityZ(
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

    using NodalRowVectorType = typename ShapeMatricesType::NodalRowVectorType;
    using GlobalDimNodalMatrixType =
        typename ShapeMatricesType::GlobalDimNodalMatrixType;
    using NodalMatrixType = typename ShapeMatricesType::NodalMatrixType;
    using NodalVectorType = typename ShapeMatricesType::NodalVectorType;
    using GlobalDimMatrixType = typename ShapeMatricesType::GlobalDimMatrixType;
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
          _saturation(
              std::vector<double>(_integration_method.getNumberOfPoints())),
          _darcy_velocities(
              GlobalDim,
              std::vector<double>(_integration_method.getNumberOfPoints()))
    {
        // This assertion is valid only if all nodal d.o.f. use the same shape
        // matrices.
        assert(local_matrix_size == ShapeFunction::NPOINTS * NUM_NODAL_DOF);
        (void)local_matrix_size;

        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();
        _ip_data.reserve(n_integration_points);

        auto const shape_matrices =
            initShapeMatrices<ShapeFunction, ShapeMatricesType,
                              IntegrationMethod, GlobalDim>(
                element, is_axially_symmetric, _integration_method);

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            auto const& sm = shape_matrices[ip];
            const double integration_factor = sm.integralMeasure * sm.detJ;
            _ip_data.emplace_back(
                sm.N, sm.dNdx,
                _integration_method.getWeightedPoint(ip).getWeight() *
                    integration_factor,
                sm.N.transpose() * sm.N * integration_factor *
                    _integration_method.getWeightedPoint(ip).getWeight());
        }
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
        auto local_b = MathLib::createZeroedVector<NodalVectorType>(
            local_b_data, local_matrix_size);

        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();
        SpatialPosition pos;
        pos.setElementID(_element.getID());
        const int material_id =
            _process_data.material->getMaterialID(pos.getElementID().get());
        const Eigen::MatrixXd& perm = _process_data.material->getPermeability(
            material_id, t, pos, _element.getDimension());
        assert(perm.rows() == _element.getDimension() || perm.rows() == 1);
        GlobalDimMatrixType permeability = GlobalDimMatrixType::Zero(
            _element.getDimension(), _element.getDimension());
        if (perm.rows() == _element.getDimension())
            permeability = perm;
        else if (perm.rows() == 1)
            permeability.diagonal().setConstant(perm(0, 0));

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            pos.setIntegrationPoint(ip);
            double p_int_pt = 0.0;
            NumLib::shapeFunctionInterpolate(local_x, _ip_data[ip].N, p_int_pt);
            const double& temperature = _process_data.temperature(t, pos)[0];
            auto const porosity = _process_data.material->getPorosity(
                material_id, t, pos, p_int_pt, temperature, 0);

            double const pc_int_pt = -p_int_pt;

            double const Sw = _process_data.material->getSaturation(
                material_id, t, pos, p_int_pt, temperature, pc_int_pt);
            _saturation[ip] = Sw;

            double const dSw_dpc =
                _process_data.material->getSaturationDerivative(
                    material_id, t, pos, p_int_pt, temperature, Sw);

            // \TODO Extend to pressure dependent density.
            double const drhow_dp(0.0);
            auto const storage = _process_data.material->getStorage(
                    material_id, t, pos, p_int_pt,temperature,0);
            double const mass_mat_coeff =
                storage * Sw + porosity * Sw * drhow_dp - porosity * dSw_dpc;

            local_M.noalias() += mass_mat_coeff * _ip_data[ip].mass_operator;

            double const k_rel =
                _process_data.material->getRelativePermeability(
                    t, pos, p_int_pt, temperature, Sw);
            auto const mu = _process_data.material->getFluidViscosity(
                p_int_pt, temperature);
            local_K.noalias() += _ip_data[ip].dNdx.transpose() *
                permeability * _ip_data[ip].dNdx *
                _ip_data[ip].integration_weight * (k_rel / mu);

            if (_process_data.has_gravity)
            {
                auto const rho_w = _process_data.material->getFluidDensity(
                    p_int_pt, temperature);
                auto const body_force = _process_data.specific_body_force;
                assert(body_force.size() == GlobalDim);
                NodalVectorType gravity_operator =
                    _ip_data[ip].dNdx.transpose() * permeability * body_force *
                    _ip_data[ip].integration_weight;
                local_b.noalias() += (k_rel / mu) * rho_w * gravity_operator;
            }
        }
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

    void computeSecondaryVariableConcrete(
        double const t, std::vector<double> const& local_x) override
    {
        SpatialPosition pos;
        pos.setElementID(_element.getID());
        const int material_id =
            _process_data.material->getMaterialID(pos.getElementID().get());

        const Eigen::MatrixXd& perm = _process_data.material->getPermeability(
            material_id, t, pos, _element.getDimension());
        assert(perm.rows() == _element.getDimension() || perm.rows() == 1);
        GlobalDimMatrixType permeability = GlobalDimMatrixType::Zero(
            _element.getDimension(), _element.getDimension());
        if (perm.rows() == _element.getDimension())
            permeability = perm;
        else if (perm.rows() == 1)
            permeability.diagonal().setConstant(perm(0, 0));

        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        auto const p_nodal_values = Eigen::Map<const NodalVectorType>(
            &local_x[0], ShapeFunction::NPOINTS);

        for (unsigned ip = 0; ip < n_integration_points; ++ip)
        {
            double p_int_pt = 0.0;
            NumLib::shapeFunctionInterpolate(local_x, _ip_data[ip].N, p_int_pt);
            double const pc_int_pt = -p_int_pt;
            const double& temperature = _process_data.temperature(t, pos)[0];
            double const Sw = _process_data.material->getSaturation(
                material_id, t, pos, p_int_pt, temperature, pc_int_pt);
            double const k_rel =
                _process_data.material->getRelativePermeability(
                    t, pos, p_int_pt, temperature, Sw);
            auto const mu = _process_data.material->getFluidViscosity(
                p_int_pt, temperature);
            GlobalDimVectorType velocity = -permeability * (k_rel / mu) *
                                           _ip_data[ip].dNdx * p_nodal_values;
            if (_process_data.has_gravity)
            {
                auto const rho_w = _process_data.material->getFluidDensity(
                    p_int_pt, temperature);
                auto const body_force = _process_data.specific_body_force;
                assert(body_force.size() == GlobalDim);
                // here it is assumed that the vector body_force is directed
                // 'downwards'
                velocity += permeability * (k_rel / mu) * rho_w * body_force;
            }

            for (unsigned d = 0; d < GlobalDim; ++d)
            {
                _darcy_velocities[d][ip] = velocity[d];
            }
        }
    }

    Eigen::Map<const Eigen::RowVectorXd> getShapeMatrix(
        const unsigned integration_point) const override
    {
        auto const& N = _ip_data[integration_point].N;

        // assumes N is stored contiguously in memory
        return Eigen::Map<const Eigen::RowVectorXd>(N.data(), N.size());
    }

    std::vector<double> const& getIntPtSaturation(
        std::vector<double>& /*cache*/) const override
    {
        assert(!_saturation.empty());
        return _saturation;
    }

    std::vector<double> const& getIntPtDarcyVelocityX(
        std::vector<double>& /*cache*/) const override
    {
        assert(!_darcy_velocities.empty());
        return _darcy_velocities[0];
    }

    std::vector<double> const& getIntPtDarcyVelocityY(
        std::vector<double>& /*cache*/) const override
    {
        assert(_darcy_velocities.size() > 1);
        return _darcy_velocities[1];
    }

    std::vector<double> const& getIntPtDarcyVelocityZ(
        std::vector<double>& /*cache*/) const override
    {
        assert(_darcy_velocities.size() > 2);
        return _darcy_velocities[2];
    }

private:
    MeshLib::Element const& _element;
    RichardsFlowProcessData const& _process_data;

    IntegrationMethod const _integration_method;
    std::vector<ShapeMatrices, Eigen::aligned_allocator<ShapeMatrices>>
        _shape_matrices;
    std::vector<
        IntegrationPointData<NodalRowVectorType, GlobalDimNodalMatrixType,
                             NodalMatrixType>,
        Eigen::aligned_allocator<IntegrationPointData<
            NodalRowVectorType, GlobalDimNodalMatrixType, NodalMatrixType>>>
        _ip_data;
    std::vector<double> _saturation;
    std::vector<std::vector<double>> _darcy_velocities;
};

}  // namespace RichardsFlow
}  // namespace ProcessLib
