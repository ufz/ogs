/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <vector>

#include "MaterialLib/MPL/MaterialSpatialDistributionMap.h"
#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Utils/FormEigenTensor.h"
#include "MaterialLib/MPL/VariableType.h"
#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "NumLib/DOF/DOFTableUtil.h"
#include "NumLib/Extrapolation/ExtrapolatableElement.h"
#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/InitShapeMatrices.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "NumLib/Function/Interpolation.h"
#include "ParameterLib/Parameter.h"
#include "ProcessLib/LocalAssemblerInterface.h"
#include "ProcessLib/LocalAssemblerTraits.h"
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

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};
const unsigned NUM_NODAL_DOF = 1;

class RichardsFlowLocalAssemblerInterface
    : public ProcessLib::LocalAssemblerInterface,
      public NumLib::ExtrapolatableElement
{
public:
    virtual std::vector<double> const& getIntPtSaturation(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtDarcyVelocity(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const = 0;
};

template <typename ShapeFunction, typename IntegrationMethod, int GlobalDim>
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
            NumLib::initShapeMatrices<ShapeFunction, ShapeMatricesType,
                                      GlobalDim>(element, is_axially_symmetric,
                                                 _integration_method);

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

    void assemble(double const t, double const dt,
                  std::vector<double> const& local_x,
                  std::vector<double> const& /*local_xdot*/,
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
        ParameterLib::SpatialPosition pos;
        pos.setElementID(_element.getID());

        auto const& medium =
            *_process_data.media_map->getMedium(_element.getID());
        auto const& liquid_phase = medium.phase("AqueousLiquid");
        MaterialPropertyLib::VariableArray vars;
        vars[static_cast<int>(MaterialPropertyLib::Variable::temperature)] =
            medium
                .property(
                    MaterialPropertyLib::PropertyType::reference_temperature)
                .template value<double>(vars, pos, t, dt);

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            pos.setIntegrationPoint(ip);
            double p_int_pt = 0.0;
            NumLib::shapeFunctionInterpolate(local_x, _ip_data[ip].N, p_int_pt);

            vars[static_cast<int>(
                MaterialPropertyLib::Variable::phase_pressure)] = p_int_pt;
            vars[static_cast<int>(
                MaterialPropertyLib::Variable::capillary_pressure)] = -p_int_pt;

            auto const permeability =
                MaterialPropertyLib::formEigenTensor<GlobalDim>(
                    medium.property(MaterialPropertyLib::permeability)
                        .value(vars, pos, t, dt));

            auto const porosity =
                medium.property(MaterialPropertyLib::PropertyType::porosity)
                    .template value<double>(vars, pos, t, dt);

            double const Sw =
                medium
                    .property(MaterialPropertyLib::PropertyType::saturation)
                    .template value<double>(vars, pos, t, dt);
            _saturation[ip] = Sw;
            vars[static_cast<int>(
                MaterialPropertyLib::Variable::liquid_saturation)] = Sw;

            double const dSw_dpc =
                medium
                    .property(MaterialPropertyLib::PropertyType::saturation)
                    .template dValue<double>(
                        vars, MaterialPropertyLib::Variable::capillary_pressure,
                        pos, t, dt);

            auto const drhow_dp =
                liquid_phase
                    .property(MaterialPropertyLib::PropertyType::density)
                    .template dValue<double>(
                        vars, MaterialPropertyLib::Variable::phase_pressure,
                        pos, t, dt);
            auto const storage =
                medium.property(MaterialPropertyLib::PropertyType::storage)
                    .template value<double>(vars, pos, t, dt);
            double const mass_mat_coeff =
                storage * Sw + porosity * Sw * drhow_dp - porosity * dSw_dpc;

            local_M.noalias() += mass_mat_coeff * _ip_data[ip].mass_operator;

            double const k_rel =
                medium
                    .property(MaterialPropertyLib::PropertyType::
                                  relative_permeability)
                    .template value<double>(vars, pos, t, dt);

            auto const mu =
                liquid_phase.property(MaterialPropertyLib::viscosity)
                    .template value<double>(vars, pos, t, dt);
            local_K.noalias() += _ip_data[ip].dNdx.transpose() * permeability *
                                 _ip_data[ip].dNdx *
                                 _ip_data[ip].integration_weight * (k_rel / mu);

            if (_process_data.has_gravity)
            {
                auto const rho_w =
                    liquid_phase
                        .property(MaterialPropertyLib::PropertyType::density)
                        .template value<double>(vars, pos, t, dt);
                auto const& body_force = _process_data.specific_body_force;
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

    std::vector<double> const& getIntPtDarcyVelocity(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const override
    {
        // TODO (tf) Temporary value not used by current material models. Need
        // extension of secondary variable interface
        double const dt = std::numeric_limits<double>::quiet_NaN();

        constexpr int process_id = 0;  // monolithic scheme.
        auto const indices =
            NumLib::getIndices(_element.getID(), *dof_table[process_id]);
        assert(!indices.empty());
        auto const local_x = x[process_id]->get(indices);

        ParameterLib::SpatialPosition pos;
        pos.setElementID(_element.getID());

        auto const& medium =
            *_process_data.media_map->getMedium(_element.getID());
        auto const& liquid_phase = medium.phase("AqueousLiquid");

        MaterialPropertyLib::VariableArray vars;
        vars[static_cast<int>(MaterialPropertyLib::Variable::temperature)] =
            medium
                .property(
                    MaterialPropertyLib::PropertyType::reference_temperature)
                .template value<double>(vars, pos, t, dt);

        auto const p_nodal_values = Eigen::Map<const NodalVectorType>(
            &local_x[0], ShapeFunction::NPOINTS);

        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        cache.clear();
        auto cache_vec = MathLib::createZeroedMatrix<
            Eigen::Matrix<double, GlobalDim, Eigen::Dynamic, Eigen::RowMajor>>(
            cache, GlobalDim, n_integration_points);

        for (unsigned ip = 0; ip < n_integration_points; ++ip)
        {
            double p_int_pt = 0.0;
            NumLib::shapeFunctionInterpolate(local_x, _ip_data[ip].N, p_int_pt);
            vars[static_cast<int>(
                MaterialPropertyLib::Variable::phase_pressure)] = p_int_pt;
            vars[static_cast<int>(
                MaterialPropertyLib::Variable::capillary_pressure)] = -p_int_pt;

            double const Sw =
                medium
                    .property(MaterialPropertyLib::PropertyType::saturation)
                    .template value<double>(vars, pos, t, dt);
            vars[static_cast<int>(
                MaterialPropertyLib::Variable::liquid_saturation)] = Sw;

            auto const permeability =
                MaterialPropertyLib::formEigenTensor<GlobalDim>(
                    medium.property(MaterialPropertyLib::permeability)
                        .value(vars, pos, t, dt));

            double const k_rel =
                medium
                    .property(MaterialPropertyLib::PropertyType::
                                  relative_permeability)
                    .template value<double>(vars, pos, t, dt);
            auto const mu =
                liquid_phase.property(MaterialPropertyLib::viscosity)
                    .template value<double>(vars, pos, t, dt);
            auto const K_mat_coeff = permeability * (k_rel / mu);
            cache_vec.col(ip).noalias() =
                -K_mat_coeff * _ip_data[ip].dNdx * p_nodal_values;
            if (_process_data.has_gravity)
            {
                auto const rho_w =
                    liquid_phase
                        .property(MaterialPropertyLib::PropertyType::density)
                        .template value<double>(vars, pos, t, dt);
                auto const& body_force = _process_data.specific_body_force;
                assert(body_force.size() == GlobalDim);
                // here it is assumed that the vector body_force is directed
                // 'downwards'
                cache_vec.col(ip).noalias() += K_mat_coeff * rho_w * body_force;
            }
        }

        return cache;
    }

private:
    MeshLib::Element const& _element;
    RichardsFlowProcessData const& _process_data;

    IntegrationMethod const _integration_method;
    std::vector<
        IntegrationPointData<NodalRowVectorType, GlobalDimNodalMatrixType,
                             NodalMatrixType>,
        Eigen::aligned_allocator<IntegrationPointData<
            NodalRowVectorType, GlobalDimNodalMatrixType, NodalMatrixType>>>
        _ip_data;
    std::vector<double> _saturation;
};

}  // namespace RichardsFlow
}  // namespace ProcessLib
