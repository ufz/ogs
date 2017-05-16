/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>
#include <vector>

#include "MaterialLib/SolidModels/KelvinVector.h"
#include "MaterialLib/SolidModels/LinearElasticIsotropic.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "NumLib/Extrapolation/ExtrapolatableElement.h"
#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "NumLib/Function/Interpolation.h"
#include "ProcessLib/Deformation/BMatrixPolicy.h"
#include "ProcessLib/Deformation/LinearBMatrix.h"
#include "ProcessLib/LocalAssemblerInterface.h"
#include "ProcessLib/LocalAssemblerTraits.h"
#include "ProcessLib/Parameter/Parameter.h"
#include "ProcessLib/Utils/InitShapeMatrices.h"

#include "ThermoHydroMechanicsProcessData.h"

namespace ProcessLib
{
namespace ThermoHydroMechanics
{
template <typename BMatricesType, typename ShapeMatrixTypeDisplacement,
          typename ShapeMatricesTypePressure, int DisplacementDim, int NPoints>
struct IntegrationPointData final
{
    explicit IntegrationPointData(
        MaterialLib::Solids::MechanicsBase<DisplacementDim>& solid_material)
        : solid_material(solid_material),
          material_state_variables(
              solid_material.createMaterialStateVariables())
    {
    }

#if defined(_MSC_VER) && _MSC_VER < 1900
    // The default generated move-ctor is correctly generated for other
    // compilers.
    explicit IntegrationPointData(IntegrationPointData&& other)
        : b_matrices(std::move(other.b_matrices)),
          sigma_eff(std::move(other.sigma_eff)),
          sigma_eff_prev(std::move(other.sigma_eff_prev)),
          eps(std::move(other.eps)),
          eps_prev(std::move(other.eps_prev)),
          eps_m(std::move(other.eps_m)),
          eps_m_prev(std::move(other.eps_m_prev)),
          solid_material(other.solid_material),
          material_state_variables(std::move(other.material_state_variables)),
          C(std::move(other.C)),
          integration_weight(std::move(other.integration_weight))
    {
    }
#endif  // _MSC_VER

    typename ShapeMatrixTypeDisplacement::template MatrixType<
        DisplacementDim, NPoints * DisplacementDim>
        N_u;
    typename BMatricesType::BMatrixType b_matrices;
    typename BMatricesType::KelvinVectorType sigma_eff, sigma_eff_prev;
    typename BMatricesType::KelvinVectorType eps, eps_prev;
    typename BMatricesType::KelvinVectorType eps_m, eps_m_prev;

    typename ShapeMatricesTypePressure::NodalRowVectorType N_p;
    typename ShapeMatricesTypePressure::GlobalDimNodalMatrixType dNdx_p;

    MaterialLib::Solids::MechanicsBase<DisplacementDim>& solid_material;
    std::unique_ptr<typename MaterialLib::Solids::MechanicsBase<
        DisplacementDim>::MaterialStateVariables>
        material_state_variables;

    typename BMatricesType::KelvinMatrixType C;
    double integration_weight;

    void pushBackState()
    {
        eps_m_prev = eps_m;
        sigma_eff_prev = sigma_eff;
        material_state_variables->pushBackState();
    }

    static const int kelvin_vector_size =
        KelvinVectorDimensions<DisplacementDim>::value;
    using Invariants = MaterialLib::SolidModels::Invariants<kelvin_vector_size>;

    template <typename DisplacementVectorType>
    void updateConstitutiveRelation(double const t,
                                    SpatialPosition const& x_position,
                                    double const dt,
                                    DisplacementVectorType const& u,
                                    double const thermal_strain)
    {
        eps.noalias() = b_matrices * u;
        // assume isotropic thermal expansion
        eps_m.noalias() = eps - thermal_strain * Invariants::identity2;
        solid_material.computeConstitutiveRelation(
            t, x_position, dt, eps_m_prev, eps_m, sigma_eff_prev, sigma_eff, C,
            *material_state_variables);
    }

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

template <typename ShapeMatrixType>
struct SecondaryData
{
    std::vector<ShapeMatrixType, Eigen::aligned_allocator<ShapeMatrixType>> N_u;
};

struct ThermoHydroMechanicsLocalAssemblerInterface
    : public ProcessLib::LocalAssemblerInterface,
      public NumLib::ExtrapolatableElement
{
public:
    virtual std::vector<double> const& getIntPtSigmaXX(
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtSigmaYY(
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtSigmaZZ(
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtSigmaXY(
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtSigmaXZ(
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtSigmaYZ(
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtEpsilonXX(
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtEpsilonYY(
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtEpsilonZZ(
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtEpsilonXY(
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtEpsilonXZ(
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtEpsilonYZ(
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtDarcyVelocityX(
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtDarcyVelocityY(
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtDarcyVelocityZ(
        std::vector<double>& cache) const = 0;
};

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
class ThermoHydroMechanicsLocalAssembler
    : public ThermoHydroMechanicsLocalAssemblerInterface
{
public:
    using ShapeMatricesTypeDisplacement =
        ShapeMatrixPolicyType<ShapeFunctionDisplacement, DisplacementDim>;

    // Types for pressure.
    using ShapeMatricesTypePressure =
        ShapeMatrixPolicyType<ShapeFunctionPressure, DisplacementDim>;
    using ShapeMatrices = typename ShapeMatricesTypeDisplacement::ShapeMatrices;

    ThermoHydroMechanicsLocalAssembler(
        ThermoHydroMechanicsLocalAssembler const&) = delete;
    ThermoHydroMechanicsLocalAssembler(ThermoHydroMechanicsLocalAssembler&&) =
        delete;

    ThermoHydroMechanicsLocalAssembler(
        MeshLib::Element const& e,
        std::size_t const /*local_matrix_size*/,
        bool is_axially_symmetric,
        unsigned const integration_order,
        ThermoHydroMechanicsProcessData<DisplacementDim>& process_data)
        : _process_data(process_data),
          _integration_method(integration_order),
          _element(e),
          _darcy_velocities(
              2, std::vector<double>(_integration_method.getNumberOfPoints()))
    {
        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        _ip_data.reserve(n_integration_points);
        _secondary_data.N_u.resize(n_integration_points);

        auto const shape_matrices_u =
            initShapeMatrices<ShapeFunctionDisplacement,
                              ShapeMatricesTypeDisplacement, IntegrationMethod,
                              DisplacementDim>(e, is_axially_symmetric,
                                               _integration_method);

        auto const shape_matrices_p =
            initShapeMatrices<ShapeFunctionPressure, ShapeMatricesTypePressure,
                              IntegrationMethod, DisplacementDim>(
                e, is_axially_symmetric, _integration_method);

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            // displacement (subscript u)
            _ip_data.emplace_back(*_process_data.material);
            auto& ip_data = _ip_data[ip];
            auto const& sm = shape_matrices_u[ip];
            _ip_data[ip].integration_weight =
                _integration_method.getWeightedPoint(ip).getWeight() *
                sm.integralMeasure * shape_matrices_u[ip].detJ;
            ip_data.b_matrices.resize(
                kelvin_vector_size,
                ShapeFunctionDisplacement::NPOINTS * DisplacementDim);

            auto const x_coord =
                interpolateXCoordinate<ShapeFunctionDisplacement,
                                       ShapeMatricesTypeDisplacement>(
                    e, shape_matrices_u[ip].N);
            LinearBMatrix::computeBMatrix<DisplacementDim,
                                          ShapeFunctionDisplacement::NPOINTS>(
                shape_matrices_u[ip].dNdx, ip_data.b_matrices,
                is_axially_symmetric, shape_matrices_u[ip].N, x_coord);

            ip_data.sigma_eff.resize(kelvin_vector_size);
            ip_data.sigma_eff_prev.resize(kelvin_vector_size);
            ip_data.eps.resize(kelvin_vector_size);
            ip_data.eps_prev.resize(kelvin_vector_size);
            ip_data.eps_m.resize(kelvin_vector_size);
            ip_data.eps_m_prev.resize(kelvin_vector_size);
            ip_data.C.resize(kelvin_vector_size, kelvin_vector_size);

            ip_data.N_u = ShapeMatricesTypeDisplacement::template MatrixType<
                DisplacementDim, displacement_size>::Zero(DisplacementDim,
                                                          displacement_size);
            for (int i = 0; i < DisplacementDim; ++i)
                ip_data.N_u
                    .template block<1, displacement_size / DisplacementDim>(
                        i, i * displacement_size / DisplacementDim)
                    .noalias() = shape_matrices_u[ip].N;

            ip_data.N_p = shape_matrices_p[ip].N;
            ip_data.dNdx_p = shape_matrices_p[ip].dNdx;
            _secondary_data.N_u[ip] = shape_matrices_u[ip].N;
        }
    }

    void assemble(double const /*t*/, std::vector<double> const& /*local_x*/,
                  std::vector<double>& /*local_M_data*/,
                  std::vector<double>& /*local_K_data*/,
                  std::vector<double>& /*local_rhs_data*/) override
    {
        OGS_FATAL(
            "ThermoHydroMechanicsLocalAssembler: assembly without Jacobian is "
            "not "
            "implemented.");
    }

    void assembleWithJacobian(double const t,
                              std::vector<double> const& local_x,
                              std::vector<double> const& local_xdot,
                              const double /*dxdot_dx*/, const double /*dx_dx*/,
                              std::vector<double>& /*local_M_data*/,
                              std::vector<double>& /*local_K_data*/,
                              std::vector<double>& local_rhs_data,
                              std::vector<double>& local_Jac_data) override
    {
        assert(local_x.size() ==
               temperature_size + pressure_size + displacement_size);

        auto T =
            Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
                temperature_size> const>(local_x.data() + temperature_index,
                                         temperature_size);

        auto p =
            Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
                pressure_size> const>(local_x.data() + pressure_index,
                                      pressure_size);

        auto u = Eigen::Map<typename ShapeMatricesTypeDisplacement::
                                template VectorType<displacement_size> const>(
            local_x.data() + displacement_index, displacement_size);

        auto T_dot =
            Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
                temperature_size> const>(local_xdot.data() + temperature_index,
                                         temperature_size);

        auto p_dot =
            Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
                pressure_size> const>(local_xdot.data() + pressure_index,
                                      pressure_size);
        auto u_dot =
            Eigen::Map<typename ShapeMatricesTypeDisplacement::
                           template VectorType<displacement_size> const>(
                local_xdot.data() + displacement_index, displacement_size);

        auto local_Jac = MathLib::createZeroedMatrix<
            typename ShapeMatricesTypeDisplacement::template MatrixType<
                displacement_size + pressure_size + temperature_size,
                displacement_size + pressure_size + temperature_size>>(
            local_Jac_data,
            displacement_size + pressure_size + temperature_size,
            displacement_size + pressure_size + temperature_size);

        auto local_rhs = MathLib::createZeroedVector<
            typename ShapeMatricesTypeDisplacement::template VectorType<
                displacement_size + pressure_size + temperature_size>>(
            local_rhs_data,
            displacement_size + pressure_size + temperature_size);

        typename ShapeMatricesTypePressure::NodalMatrixType MTT;
        MTT.setZero(temperature_size, temperature_size);

        typename ShapeMatricesTypePressure::NodalMatrixType KTT_coeff;
        KTT_coeff.setZero(temperature_size, temperature_size);

        typename ShapeMatricesTypePressure::NodalMatrixType KTT;
        KTT.setZero(temperature_size, temperature_size);

        typename ShapeMatricesTypePressure::NodalMatrixType KTp;
        KTp.setZero(temperature_size, pressure_size);

        typename ShapeMatricesTypePressure::NodalMatrixType KTp_coeff;
        KTp_coeff.setZero(temperature_size, pressure_size);

        typename ShapeMatricesTypePressure::NodalMatrixType laplace_p;
        laplace_p.setZero(pressure_size, pressure_size);

        typename ShapeMatricesTypePressure::NodalMatrixType storage_p;
        storage_p.setZero(pressure_size, pressure_size);

        typename ShapeMatricesTypePressure::NodalMatrixType storage_T;
        storage_T.setZero(pressure_size, temperature_size);

        typename ShapeMatricesTypeDisplacement::template MatrixType<
            displacement_size, pressure_size>
            Kup;
        Kup.setZero(displacement_size, pressure_size);
        typename ShapeMatricesTypeDisplacement::template MatrixType<
            displacement_size, temperature_size>
            KuT;
        KuT.setZero(displacement_size, temperature_size);

        double const& dt = _process_data.dt;

        SpatialPosition x_position;
        x_position.setElementID(_element.getID());

        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();
        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            x_position.setIntegrationPoint(ip);
            auto const& w = _ip_data[ip].integration_weight;

            auto const& N_p = _ip_data[ip].N_p;
            auto const& N_u = _ip_data[ip].N_u;
            auto const& dNdx_p = _ip_data[ip].dNdx_p;

            auto const& N_T = N_p;
            auto const& dNdx_T = dNdx_p;

            auto const& B = _ip_data[ip].b_matrices;
            auto const& sigma_eff = _ip_data[ip].sigma_eff;

            auto const& C = _ip_data[ip].C;

            double const S =
                _process_data.storage_coefficient(t, x_position)[0];
            double const K_over_mu =
                _process_data.intrinsic_permeability(t, x_position)[0] /
                _process_data.fluid_viscosity(t, x_position)[0];
            double const alpha_s = _process_data.solid_linear_thermal_expansion_coefficient(t, x_position)[0];
            double const beta_f = _process_data.fluid_volumetric_thermal_expansion_coefficient(t, x_position)[0];
            double const lambda_f = _process_data.fluid_thermal_conductivity(t, x_position)[0];
            double const lambda_s = _process_data.solid_thermal_conductivity(t, x_position)[0];
            double const C_f =
                _process_data.fluid_specific_heat_capacity(t, x_position)[0];
            double const C_s =
                _process_data.solid_specific_heat_capacity(t, x_position)[0];
            double const T0 =
                _process_data.reference_temperature(t, x_position)[0];
            auto const alpha = _process_data.biot_coefficient(t, x_position)[0];
            auto rho_sr = _process_data.solid_density(t, x_position)[0];
            auto rho_fr = _process_data.fluid_density(t, x_position)[0];
            auto const porosity = _process_data.porosity(t, x_position)[0];
            auto const& b = _process_data.specific_body_force;

            // calculate linear thermal strain
            // assume isotropic thermal expansion
            auto T_int_pt = N_T * T;
            double delta_T(T_int_pt - T0);
            double const thermal_strain = alpha_s * delta_T;

            double rho_f = rho_fr * (1 - beta_f * delta_T);
            double rho_s = rho_sr * (1 - 3 * thermal_strain);

            auto velocity = (-K_over_mu * dNdx_p * p).eval();
            velocity += K_over_mu * rho_f * b;

            // mapping matrix (m)
            auto const& identity2 = MaterialLib::SolidModels::Invariants<
                kelvin_vector_size>::identity2;

            //
            // displacement equation, displacement part (K_uu)
            //
            _ip_data[ip].updateConstitutiveRelation(t, x_position, dt, u,
                                                    thermal_strain);

            local_Jac
                .template block<displacement_size, displacement_size>(
                    displacement_index, displacement_index)
                .noalias() += B.transpose() * C * B * w;

            double const rho = rho_s * (1 - porosity) + porosity * rho_f;
            local_rhs.template segment<displacement_size>(displacement_index)
                .noalias() -=
                (B.transpose() * sigma_eff - N_u.transpose() * rho * b) * w;

            //
            // displacement equation, pressure part (K_up)
            //
            Kup.noalias() += B.transpose() * alpha * identity2 * N_p * w;

            //
            // pressure equation, pressure part (K_pp and M_pp).
            //
            laplace_p.noalias() += dNdx_p.transpose() * K_over_mu * dNdx_p * w;

            storage_p.noalias() += N_p.transpose() * S * N_p * w;
            //
            // fp
            //
            local_rhs.template segment<pressure_size>(pressure_index)
                .noalias() += dNdx_p.transpose() * rho_fr * K_over_mu * b * w;
            //
            // pressure equation, temperature part (M_pT)
            //
            auto beta = porosity * beta_f + (1-porosity)* 3 * alpha_s;
            storage_T.noalias() += N_p.transpose() * beta * N_p * w;

            //
            // pressure equation, displacement part.  (M_pu)
            //
            // Reusing Kup.transpose().

            //
            // temperature equation, temperature part.
            //
            double lambda = porosity * lambda_f + (1 - porosity) * lambda_s;
            KTT.noalias() +=
                (dNdx_T.transpose() * lambda * dNdx_T +
                 dNdx_T.transpose() * velocity * N_p * rho_fr * C_f) *
                w;
            // coeff matrix using for RHS
            double heat_capacity =
                porosity * C_f * rho_fr + (1 - porosity) * C_s * rho_sr;
            MTT.noalias() += N_T.transpose() * heat_capacity * N_T * w;

            //
            // temperature equation, pressure part !!!!!!.positive or negative
            //
            KTp.noalias() += K_over_mu * rho_fr * C_f * dNdx_p.transpose() *
                             (dNdx_T * T) * N_T * w * 0;
        }
        // temperature equation, temperature part
        local_Jac
            .template block<temperature_size, temperature_size>(
                temperature_index, temperature_index)
            .noalias() += KTT + MTT / dt;

        // temperature equation, pressure part
        local_Jac
            .template block<temperature_size, pressure_size>(temperature_index,
                                                             pressure_index)
            .noalias() -= KTp;
        // displacement equation, temperature part
        local_Jac
            .template block<displacement_size, temperature_size>(
                displacement_index, temperature_index)
            .noalias() -= KuT;

        // displacement equation, pressure part
        local_Jac
            .template block<displacement_size, pressure_size>(
                displacement_index, pressure_index)
            .noalias() -= Kup;

        // pressure equation, temperature part.
        local_Jac
            .template block<pressure_size, temperature_size>(pressure_index,
                                                             temperature_index)
            .noalias() -= storage_T / dt;

        // pressure equation, pressure part.
        local_Jac
            .template block<pressure_size, pressure_size>(pressure_index,
                                                          pressure_index)
            .noalias() += laplace_p + storage_p / dt;

        // pressure equation, displacement part.
        local_Jac
            .template block<pressure_size, displacement_size>(
                pressure_index, displacement_index)
            .noalias() += Kup.transpose() / dt;

        // pressure equation (f_p)
        local_rhs.template segment<pressure_size>(pressure_index).noalias() -=
            laplace_p * p + storage_p * p_dot - storage_T * T_dot +
            Kup.transpose() * u_dot;

        // displacement equation (f_u)
        local_rhs.template segment<displacement_size>(displacement_index)
            .noalias() += Kup * p;

        // temperature equation (f_T)
        local_rhs.template segment<temperature_size>(temperature_index)
            .noalias() -= KTT * T + MTT * T_dot;
    }

    void postTimestepConcrete(std::vector<double> const& local_x) override
    {
        double const& t = _process_data.t;

        auto p =
            Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
                pressure_size> const>(local_x.data() + pressure_index,
                                      pressure_size);
        using GlobalDimVectorType =
            typename ShapeMatricesTypePressure::GlobalDimVectorType;

        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        SpatialPosition x_position;
        x_position.setElementID(_element.getID());
        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            x_position.setIntegrationPoint(ip);
            double const K_over_mu =
                _process_data.intrinsic_permeability(t, x_position)[0] /
                _process_data.fluid_viscosity(t, x_position)[0];

            auto const rho_fr = _process_data.fluid_density(t, x_position)[0];
            auto const& b = _process_data.specific_body_force;

            // Compute the velocity
            auto const& dNdx_p = _ip_data[ip].dNdx_p;
            GlobalDimVectorType const darcy_velocity =
                -K_over_mu * dNdx_p * p - K_over_mu * rho_fr * b;
            for (unsigned d = 0; d < DisplacementDim; ++d)
            {
                _darcy_velocities[d][ip] = darcy_velocity[d];
            }
        }
    }

    void preTimestepConcrete(std::vector<double> const& /*local_x*/,
                             double const /*t*/,
                             double const /*delta_t*/) override
    {
        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            _ip_data[ip].pushBackState();
        }
    }

    Eigen::Map<const Eigen::RowVectorXd> getShapeMatrix(
        const unsigned integration_point) const override
    {
        auto const& N_u = _secondary_data.N_u[integration_point];

        // assumes N is stored contiguously in memory
        return Eigen::Map<const Eigen::RowVectorXd>(N_u.data(), N_u.size());
    }

    std::vector<double> const& getIntPtSigmaXX(
        std::vector<double>& cache) const override
    {
        return getIntPtSigma(cache, 0);
    }

    std::vector<double> const& getIntPtSigmaYY(
        std::vector<double>& cache) const override
    {
        return getIntPtSigma(cache, 1);
    }

    std::vector<double> const& getIntPtSigmaZZ(
        std::vector<double>& cache) const override
    {
        return getIntPtSigma(cache, 2);
    }

    std::vector<double> const& getIntPtSigmaXY(
        std::vector<double>& cache) const override
    {
        return getIntPtSigma(cache, 3);
    }

    std::vector<double> const& getIntPtSigmaXZ(
        std::vector<double>& cache) const override
    {
        assert(DisplacementDim == 3);
        return getIntPtSigma(cache, 4);
    }

    std::vector<double> const& getIntPtSigmaYZ(
        std::vector<double>& cache) const override
    {
        assert(DisplacementDim == 3);
        return getIntPtSigma(cache, 5);
    }

    std::vector<double> const& getIntPtEpsilonXX(
        std::vector<double>& cache) const override
    {
        return getIntPtEpsilon(cache, 0);
    }

    std::vector<double> const& getIntPtEpsilonYY(
        std::vector<double>& cache) const override
    {
        return getIntPtEpsilon(cache, 1);
    }

    std::vector<double> const& getIntPtEpsilonZZ(
        std::vector<double>& cache) const override
    {
        return getIntPtEpsilon(cache, 2);
    }

    std::vector<double> const& getIntPtEpsilonXY(
        std::vector<double>& cache) const override
    {
        return getIntPtEpsilon(cache, 3);
    }

    std::vector<double> const& getIntPtEpsilonXZ(
        std::vector<double>& cache) const override
    {
        assert(DisplacementDim == 3);
        return getIntPtEpsilon(cache, 4);
    }

    std::vector<double> const& getIntPtEpsilonYZ(
        std::vector<double>& cache) const override
    {
        assert(DisplacementDim == 3);
        return getIntPtEpsilon(cache, 5);
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
    std::vector<double> const& getIntPtSigma(std::vector<double>& cache,
                                             std::size_t const component) const
    {
        cache.clear();
        cache.reserve(_ip_data.size());

        for (auto const& ip_data : _ip_data)
        {
            if (component < 3)  // xx, yy, zz components
                cache.push_back(ip_data.sigma_eff[component]);
            else  // mixed xy, yz, xz components
                cache.push_back(ip_data.sigma_eff[component] / std::sqrt(2));
        }

        return cache;
    }

    std::vector<double> const& getIntPtEpsilon(
        std::vector<double>& cache, std::size_t const component) const
    {
        cache.clear();
        cache.reserve(_ip_data.size());

        for (auto const& ip_data : _ip_data)
        {
            if (component < 3)  // xx, yy, zz components
                cache.push_back(ip_data.eps[component]);
            else  // mixed xy, yz, xz components
                cache.push_back(ip_data.eps[component] / std::sqrt(2));
        }

        return cache;
    }

    ThermoHydroMechanicsProcessData<DisplacementDim>& _process_data;

    using BMatricesType =
        BMatrixPolicyType<ShapeFunctionDisplacement, DisplacementDim>;
    using IpData =
        IntegrationPointData<BMatricesType, ShapeMatricesTypeDisplacement,
                             ShapeMatricesTypePressure, DisplacementDim,
                             ShapeFunctionDisplacement::NPOINTS>;
    std::vector<IpData, Eigen::aligned_allocator<IpData>> _ip_data;

    IntegrationMethod _integration_method;
    MeshLib::Element const& _element;
    SecondaryData<typename ShapeMatrices::ShapeType> _secondary_data;
    std::vector<std::vector<double>> _darcy_velocities =
        std::vector<std::vector<double>>(
            DisplacementDim,
            std::vector<double>(_integration_method.getNumberOfPoints()));

    static const int temperature_index = 0;
    static const int temperature_size = ShapeFunctionPressure::NPOINTS;
    static const int pressure_index = ShapeFunctionPressure::NPOINTS;
    static const int pressure_size = ShapeFunctionPressure::NPOINTS;
    static const int displacement_index = ShapeFunctionPressure::NPOINTS * 2;
    static const int displacement_size =
        ShapeFunctionDisplacement::NPOINTS * DisplacementDim;
    static const int kelvin_vector_size =
        KelvinVectorDimensions<DisplacementDim>::value;
    // std::vector<std::vector<double>> _darcy_velocities;
};

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, unsigned GlobalDim, int DisplacementDim>
class LocalAssemblerData final
    : public ThermoHydroMechanicsLocalAssembler<
          ShapeFunctionDisplacement, ShapeFunctionPressure, IntegrationMethod,
          DisplacementDim>
{
public:
    LocalAssemblerData(LocalAssemblerData const&) = delete;
    LocalAssemblerData(LocalAssemblerData&&) = delete;

    LocalAssemblerData(
        MeshLib::Element const& e,
        std::size_t const local_matrix_size,
        bool is_axially_symmetric,
        unsigned const integration_order,
        ThermoHydroMechanicsProcessData<DisplacementDim>& process_data)
        : ThermoHydroMechanicsLocalAssembler<
              ShapeFunctionDisplacement, ShapeFunctionPressure,
              IntegrationMethod, DisplacementDim>(
              e, local_matrix_size, is_axially_symmetric, integration_order,
              process_data)
    {
    }
};

}  // namespace ThermoHydroMechanics
}  // namespace ProcessLib
