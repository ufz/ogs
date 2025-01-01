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

#include <Eigen/Dense>
#include <numeric>
#include <vector>

#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Utils/FormEffectiveThermalConductivity.h"
#include "NumLib/DOF/DOFTableUtil.h"
#include "NumLib/Extrapolation/ExtrapolatableElement.h"
#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/InitShapeMatrices.h"
#include "NumLib/Fem/Integration/GenericIntegrationMethod.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "NumLib/NewtonRaphson.h"
#include "ParameterLib/Parameter.h"
#include "ProcessLib/Utils/SetOrGetIntegrationPointData.h"
#include "WellboreSimulatorLocalAssemblerInterface.h"
#include "WellboreSimulatorProcessData.h"

namespace ProcessLib
{
namespace WellboreSimulator
{
const unsigned NUM_NODAL_DOF = 3;

template <typename ShapeFunction, int GlobalDim>
class WellboreSimulatorFEM : public WellboreSimulatorLocalAssemblerInterface
{
    using ShapeMatricesType = ShapeMatrixPolicyType<ShapeFunction, GlobalDim>;
    using ShapeMatrices = typename ShapeMatricesType::ShapeMatrices;

    using LocalMatrixType = typename ShapeMatricesType::template MatrixType<
        NUM_NODAL_DOF * ShapeFunction::NPOINTS,
        NUM_NODAL_DOF * ShapeFunction::NPOINTS>;
    using LocalVectorType =
        typename ShapeMatricesType::template VectorType<NUM_NODAL_DOF *
                                                        ShapeFunction::NPOINTS>;

    using NodalVectorType = typename ShapeMatricesType::NodalVectorType;
    using NodalRowVectorType = typename ShapeMatricesType::NodalRowVectorType;

    using GlobalDimVectorType = typename ShapeMatricesType::GlobalDimVectorType;
    using GlobalDimNodalMatrixType =
        typename ShapeMatricesType::GlobalDimNodalMatrixType;
    using GlobalDimMatrixType = typename ShapeMatricesType::GlobalDimMatrixType;

public:
    WellboreSimulatorFEM(
        MeshLib::Element const& element,
        std::size_t const /*local_matrix_size*/,
        NumLib::GenericIntegrationMethod const& integration_method,
        bool const is_axially_symmetric,
        WellboreSimulatorProcessData const& process_data)
        : _element(element),
          _integration_method(integration_method),
          _is_axially_symmetric(is_axially_symmetric),
          _process_data(process_data)
    {
        // calculate the element direction vector
        auto const& p0 = element.getNode(0)->asEigenVector3d();
        auto const& p1 = element.getNode(1)->asEigenVector3d();

        _element_direction = (p1 - p0).normalized();

        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();
        _ip_data.reserve(n_integration_points);

        auto const shape_matrices =
            NumLib::initShapeMatrices<ShapeFunction, ShapeMatricesType,
                                      GlobalDim>(element, is_axially_symmetric,
                                                 _integration_method);

        ParameterLib::SpatialPosition pos;
        pos.setElementID(_element.getID());
        auto const& medium =
            *_process_data.media_map.getMedium(_element.getID());
        auto const& liquid_phase = medium.phase("AqueousLiquid");

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            _ip_data.emplace_back(
                shape_matrices[ip].N, shape_matrices[ip].dNdx,
                _integration_method.getWeightedPoint(ip).getWeight() *
                    shape_matrices[ip].integralMeasure *
                    shape_matrices[ip].detJ);

            pos.setIntegrationPoint(ip);
            MaterialPropertyLib::VariableArray vars;

            NodalVectorType ref_p =
                _process_data.well_ref_pressure.getNodalValuesOnElement(
                    _element, 0);
            NodalVectorType ref_h =
                _process_data.well_ref_enthalpy.getNodalValuesOnElement(
                    _element, 0);

            vars.liquid_phase_pressure = _ip_data[ip].N.dot(ref_p);
            vars.enthalpy = _ip_data[ip].N.dot(ref_h);

            //  .initialValue
            _ip_data[ip].temperature =
                liquid_phase
                    .property(MaterialPropertyLib::PropertyType::temperature)
                    .template value<double>(vars, pos, 0, 0);
            vars.temperature = _ip_data[ip].temperature;
            _ip_data[ip].mix_density =
                liquid_phase
                    .property(MaterialPropertyLib::PropertyType::density)
                    .template value<double>(vars, pos, 0, 0);
            _ip_data[ip].dryness = 0;
            _ip_data[ip].vapor_volume_fraction = 0;
            _ip_data[ip].vapor_mass_flow_rate = 0;
            _ip_data[ip].liquid_mass_flow_rate = 0;
            _ip_data[ip].pushBackState();
        }
    }

    void assemble(double const t, double const dt,
                  std::vector<double> const& local_x,
                  std::vector<double> const& local_x_prev,
                  std::vector<double>& local_M_data,
                  std::vector<double>& local_K_data,
                  std::vector<double>& local_b_data) override;

    static int const jacobian_residual_size = 1;
    using ResidualVector = Eigen::Matrix<double, jacobian_residual_size, 1>;
    using JacobianMatrix =
        Eigen::Matrix<double, jacobian_residual_size, jacobian_residual_size,
                      Eigen::RowMajor>;
    using UnknownVector = Eigen::Matrix<double, jacobian_residual_size, 1>;

    void calculateResidual(double const alpha, double const vapor_water_density,
                           double const liquid_water_density,
                           double const v_mix, double const dryness,
                           double const C_0, double const u_gu,
                           ResidualVector& res)
    {
        double const rho_mix =
            alpha * vapor_water_density + (1 - alpha) * liquid_water_density;

        res(0) =
            dryness * liquid_water_density * rho_mix * v_mix -
            alpha * C_0 * dryness * liquid_water_density * rho_mix * v_mix -
            alpha * C_0 * (1 - dryness) * vapor_water_density * rho_mix *
                v_mix -
            alpha * vapor_water_density * liquid_water_density * u_gu;
    }

    void calculateJacobian(double const alpha, double const vapor_water_density,
                           double const liquid_water_density,
                           double const v_mix, double const dryness,
                           double const C_0, double const u_gu,
                           JacobianMatrix& Jac)
    {
        Jac(0) = dryness * liquid_water_density * v_mix *
                     (vapor_water_density - liquid_water_density) -
                 (C_0 * dryness * liquid_water_density +
                  C_0 * (1 - dryness) * vapor_water_density) *
                     (2 * alpha * vapor_water_density +
                      (1 - 2 * alpha) * liquid_water_density) *
                     v_mix -
                 vapor_water_density * liquid_water_density * u_gu;
    }

    Eigen::Map<const Eigen::RowVectorXd> getShapeMatrix(
        const unsigned integration_point) const override
    {
        auto const& N = _ip_data[integration_point].N;

        // assumes N is stored contiguously in memory
        return Eigen::Map<const Eigen::RowVectorXd>(N.data(), N.size());
    }

    void computeSecondaryVariableConcrete(
        double const /*t*/,
        double const /*dt*/,
        Eigen::VectorXd const& /*local_x*/,
        Eigen::VectorXd const& /*local_x_dot*/) override
    {
        auto const n_integration_points =
            _integration_method.getNumberOfPoints();
        auto const ele_id = _element.getID();

        (*_process_data.mesh_prop_density)[ele_id] =
            std::accumulate(_ip_data.begin(), _ip_data.end(), 0.,
                            [](double const s, auto const& ip)
                            { return s + ip.mix_density; }) /
            n_integration_points;
    }

    void postTimestepConcrete(Eigen::VectorXd const& /*local_x*/,
                              Eigen::VectorXd const& /*local_x_dot*/,
                              double const /*t*/, double const /*dt*/,
                              int const /*process_id*/) override
    {
        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            _ip_data[ip].pushBackState();
        }
    }

protected:
    virtual std::vector<double> const& getIntPtVaporMassFlowRate(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        return ProcessLib::getIntegrationPointScalarData(
            _ip_data, &IpData::vapor_mass_flow_rate, cache);
    }

    virtual std::vector<double> const& getIntPtLiquidMassFlowRate(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        return ProcessLib::getIntegrationPointScalarData(
            _ip_data, &IpData::liquid_mass_flow_rate, cache);
    }

    virtual std::vector<double> const& getIntPtTemperature(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        return ProcessLib::getIntegrationPointScalarData(
            _ip_data, &IpData::temperature, cache);
    }

    virtual std::vector<double> const& getIntPtDryness(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        return ProcessLib::getIntegrationPointScalarData(
            _ip_data, &IpData::dryness, cache);
    }

    virtual std::vector<double> const& getIntPtVaporVolumeFraction(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        return ProcessLib::getIntegrationPointScalarData(
            _ip_data, &IpData::vapor_volume_fraction, cache);
    }

    MeshLib::Element const& _element;
    NumLib::GenericIntegrationMethod const& _integration_method;
    bool const _is_axially_symmetric;
    WellboreSimulatorProcessData const& _process_data;

    Eigen::Vector3d _element_direction;

    using IpData =
        IntegrationPointData<NodalRowVectorType, GlobalDimNodalMatrixType>;
    std::vector<IpData, Eigen::aligned_allocator<IpData>> _ip_data;

protected:
    static const int pressure_index = 0;
    static const int velocity_index = ShapeFunction::NPOINTS;
    static const int enthalpy_index = 2 * ShapeFunction::NPOINTS;

    static const int pressure_size = ShapeFunction::NPOINTS;
    static const int velocity_size = ShapeFunction::NPOINTS;
    static const int enthalpy_size = ShapeFunction::NPOINTS;
};

}  // namespace WellboreSimulator
}  // namespace ProcessLib

#include "WellboreSimulatorFEM-impl.h"
