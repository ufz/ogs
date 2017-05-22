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
#include <Eigen/Eigenvalues>

#include "MaterialLib/SolidModels/KelvinVector.h"
#include "MaterialLib/SolidModels/LinearElasticIsotropicPhaseField.h"
#include "MaterialLib/SolidModels/LinearElasticIsotropic.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "NumLib/Extrapolation/ExtrapolatableElement.h"
#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "ProcessLib/Deformation/BMatrixPolicy.h"
#include "ProcessLib/Deformation/LinearBMatrix.h"
#include "ProcessLib/LocalAssemblerInterface.h"
#include "ProcessLib/LocalAssemblerTraits.h"
#include "ProcessLib/Parameter/Parameter.h"
#include "ProcessLib/Utils/InitShapeMatrices.h"

#include "PhaseFieldProcessData.h"

namespace ProcessLib
{
namespace PhaseField
{
template <typename BMatricesType, typename ShapeMatrixType, int DisplacementDim>
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
          sigma(std::move(other.sigma)),
          sigma_prev(std::move(other.sigma_prev)),
          eps(std::move(other.eps)),
          eps_prev(std::move(other.eps_prev)),
          solid_material(other.solid_material),
          material_state_variables(std::move(other.material_state_variables)),
          C(std::move(other.C)),
          C_tensile(std::move(other.C_tensile)),
          C_compressive(std::move(other.C_compressive)),
          integration_weight(std::move(other.integration_weight)),
          strain_energy_tensile(std::move(other.strain_energy_tensile)),
          history_variable(std::move(other.history_variable)),
          history_variable_prev(std::move(other.history_variable_prev)),
          sigma_tensile(std::move(other.sigma_tensile)),
          sigma_compressive(std::move(other.sigma_compressive)),
          sigma_real(std::move(other.sigma_real)),
          sigma_real_prev(std::move(other.sigma_real_prev))
    {
    }
#endif  // _MSC_VER

    typename ShapeMatrixType::NodalRowVectorType N;
    typename ShapeMatrixType::GlobalDimNodalMatrixType dNdx;
    typename BMatricesType::BMatrixType b_matrices;
    typename BMatricesType::KelvinVectorType sigma, sigma_prev;
    typename BMatricesType::KelvinVectorType eps, eps_prev;

    typename BMatricesType::KelvinVectorType sigma_tensile, sigma_compressive,
        sigma_real_prev, sigma_real;
    double strain_energy_tensile;

    MaterialLib::Solids::MechanicsBase<DisplacementDim>& solid_material;
    std::unique_ptr<typename MaterialLib::Solids::MechanicsBase<
        DisplacementDim>::MaterialStateVariables> material_state_variables;

    typename BMatricesType::KelvinMatrixType C, C_tensile, C_compressive;
    double integration_weight;
    double history_variable;
    double history_variable_prev;

    void pushBackState()
    {
        eps_prev = eps;
        sigma_real_prev = sigma_real;
        material_state_variables->pushBackState();
    }

    template <typename DisplacementVectorType>
    void updateConstitutiveRelation(double const t,
                                    SpatialPosition const& x_position,
                                    double const dt,
                                    DisplacementVectorType const& u,
                                    double const degradation)
    {
        eps.noalias() = b_matrices * u;
        solid_material.computeConstitutiveRelation(t, x_position, dt, eps_prev,
                                                   eps, sigma_prev, sigma, C,
                                                   *material_state_variables);

        static_cast<MaterialLib::Solids::PhaseFieldExtension<DisplacementDim>&>(
            solid_material)
            .specialFunction(t, x_position, eps, strain_energy_tensile,
                             sigma_tensile, sigma_compressive, C_tensile,
                             C_compressive, sigma_real, degradation);
    }
};

/// Used by for extrapolation of the integration point values. It is ordered
/// (and stored) by integration points.
template <typename ShapeMatrixType>
struct SecondaryData
{
    std::vector<ShapeMatrixType> N;
};

struct PhaseFieldLocalAssemblerInterface
    : public ProcessLib::LocalAssemblerInterface,
      public NumLib::ExtrapolatableElement
{
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
};

template <typename ShapeFunction, typename IntegrationMethod,
          int DisplacementDim>
class PhaseFieldLocalAssembler : public PhaseFieldLocalAssemblerInterface
{
public:
    using ShapeMatricesType =
        ShapeMatrixPolicyType<ShapeFunction, DisplacementDim>;

    // Types for displacement.
    // (Higher order elements = ShapeFunction).
    using ShapeMatrices = typename ShapeMatricesType::ShapeMatrices;
    using BMatricesType = BMatrixPolicyType<ShapeFunction, DisplacementDim>;

    using NodalForceVectorType = typename BMatricesType::NodalForceVectorType;
    using RhsVector = typename ShapeMatricesType::template VectorType<
        ShapeFunction::NPOINTS + ShapeFunction::NPOINTS * DisplacementDim>;
    using JacobianMatrix = typename ShapeMatricesType::template MatrixType<
        ShapeFunction::NPOINTS + ShapeFunction::NPOINTS * DisplacementDim,
        ShapeFunction::NPOINTS + ShapeFunction::NPOINTS * DisplacementDim>;

    PhaseFieldLocalAssembler(PhaseFieldLocalAssembler const&) = delete;
    PhaseFieldLocalAssembler(PhaseFieldLocalAssembler&&) = delete;

    PhaseFieldLocalAssembler(
        MeshLib::Element const& e,
        std::size_t const /*local_matrix_size*/,
        bool is_axially_symmetric,
        unsigned const integration_order,
        PhaseFieldProcessData<DisplacementDim>& process_data)
        : _process_data(process_data),
          _integration_method(integration_order),
          _element(e)
    {
        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        _ip_data.reserve(n_integration_points);
        _secondary_data.N.resize(n_integration_points);

        auto const shape_matrices =
            initShapeMatrices<ShapeFunction, ShapeMatricesType,
                              IntegrationMethod, DisplacementDim>(
                e, is_axially_symmetric, _integration_method);

        SpatialPosition x_position;
        x_position.setElementID(_element.getID());

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            // displacement (subscript u)
            _ip_data.emplace_back(*_process_data.material);
            auto& ip_data = _ip_data[ip];
            ip_data.integration_weight =
                _integration_method.getWeightedPoint(ip).getWeight() *
                shape_matrices[ip].detJ;
            ip_data.b_matrices.resize(kelvin_vector_size,
                                      ShapeFunction::NPOINTS * DisplacementDim);

            auto const x_coord =
                interpolateXCoordinate<ShapeFunction, ShapeMatricesType>(
                    e, shape_matrices[ip].N);
            LinearBMatrix::computeBMatrix<DisplacementDim,
                                          ShapeFunction::NPOINTS>(
                shape_matrices[ip].dNdx, ip_data.b_matrices,
                is_axially_symmetric, shape_matrices[ip].N, x_coord);

            ip_data.sigma.resize(kelvin_vector_size);
            ip_data.sigma_prev.resize(kelvin_vector_size);
            ip_data.eps.resize(kelvin_vector_size);
            ip_data.eps_prev.resize(kelvin_vector_size);
            ip_data.C.resize(kelvin_vector_size, kelvin_vector_size);
            ip_data.C_tensile.resize(kelvin_vector_size, kelvin_vector_size);
            ip_data.C_compressive.resize(kelvin_vector_size,
                                         kelvin_vector_size);
            ip_data.sigma_tensile.resize(kelvin_vector_size);
            ip_data.sigma_compressive.resize(kelvin_vector_size);
            ip_data.strain_energy_tensile;
            ip_data.history_variable =
                process_data.history_field(0, x_position)[0];
            ip_data.history_variable_prev =
                process_data.history_field(0, x_position)[0];
            ip_data.sigma_real.resize(kelvin_vector_size);

            ip_data.N = shape_matrices[ip].N;
            ip_data.dNdx = shape_matrices[ip].dNdx;

            _secondary_data.N[ip] = shape_matrices[ip].N;
        }
    }

    void assemble(double const /*t*/, std::vector<double> const& /*local_x*/,
                  std::vector<double>& /*local_M_data*/,
                  std::vector<double>& /*local_K_data*/,
                  std::vector<double>& /*local_rhs_data*/) override
    {
        OGS_FATAL(
            "PhaseFieldLocalAssembler: assembly without jacobian is not "
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
        auto const local_matrix_size = local_x.size();
        assert(local_matrix_size == phasefield_size + displacement_size);

        auto d = Eigen::Map<typename ShapeMatricesType::template VectorType<
            phasefield_size> const>(local_x.data() + phasefield_index,
                                    phasefield_size);

        auto u = Eigen::Map<typename ShapeMatricesType::template VectorType<
            displacement_size> const>(local_x.data() + displacement_index,
                                      displacement_size);

        auto d_dot = Eigen::Map<typename ShapeMatricesType::template VectorType<
            phasefield_size> const>(local_xdot.data() + phasefield_index,
                                    phasefield_size);

        auto local_Jac = MathLib::createZeroedMatrix<JacobianMatrix>(
            local_Jac_data, local_matrix_size, local_matrix_size);

        auto local_rhs = MathLib::createZeroedVector<RhsVector>(
            local_rhs_data, local_matrix_size);

        typename ShapeMatricesType::template MatrixType<displacement_size,
                                                        phasefield_size> Kud;
        Kud.setZero(displacement_size, phasefield_size);

        typename ShapeMatricesType::template MatrixType<phasefield_size,
                                                        displacement_size> Kdu;
        Kdu.setZero(phasefield_size, displacement_size);

        typename ShapeMatricesType::NodalMatrixType Kdd;
        Kdd.setZero(phasefield_size, phasefield_size);

        typename ShapeMatricesType::NodalMatrixType Ddd;
        Ddd.setZero(phasefield_size, phasefield_size);

        double const& dt = _process_data.dt;

        SpatialPosition x_position;
        x_position.setElementID(_element.getID());

        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();
        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            x_position.setIntegrationPoint(ip);
            auto const& w = _ip_data[ip].integration_weight;

            auto const& dNdx = _ip_data[ip].dNdx;
            auto const& N = _ip_data[ip].N;

            auto const& B = _ip_data[ip].b_matrices;

            auto const& C_tensile = _ip_data[ip].C_tensile;
            auto const& C_compressive = _ip_data[ip].C_compressive;

            auto const& strain_energy_tensile =
                _ip_data[ip].strain_energy_tensile;
            auto const& sigma_tensile = _ip_data[ip].sigma_tensile;

            auto& history_variable = _ip_data[ip].history_variable;
            auto& history_variable_prev = _ip_data[ip].history_variable_prev;

            auto const& sigma_real = _ip_data[ip].sigma_real;

            // Kdd_1 defines one term which both used in Kdd and local_rhs for
            // phase field
            double const gc = _process_data.crack_resistance(t, x_position)[0];
            double const ls =
                _process_data.crack_length_scale(t, x_position)[0];
            typename ShapeMatricesType::NodalMatrixType const Kdd_1 =
                dNdx.transpose() * 2 * gc * ls * dNdx;

            //
            // displacement equation, displacement part
            //
            double const k = _process_data.residual_stiffness(t, x_position)[0];
            double const d_ip = N.dot(d);
            double const degradation = d_ip * d_ip * (1 - k) + k;
            _ip_data[ip].updateConstitutiveRelation(t, x_position, dt, u,
                                                    degradation);

            local_Jac.template block<displacement_size, displacement_size>(
                         displacement_index, displacement_index)
                .noalias() += B.transpose() *
                              (degradation * C_tensile + C_compressive) * B * w;

            typename ShapeMatricesType::template MatrixType<
                DisplacementDim, displacement_size> N_u = ShapeMatricesType::
                template MatrixType<DisplacementDim, displacement_size>::Zero(
                    DisplacementDim, displacement_size);

            for (int i = 0; i < DisplacementDim; ++i)
                N_u.template block<1, displacement_size / DisplacementDim>(
                       i, i * displacement_size / DisplacementDim)
                    .noalias() = N;

            auto const rho_sr = _process_data.solid_density(t, x_position)[0];
            auto const& b = _process_data.specific_body_force;
            local_rhs.template segment<displacement_size>(displacement_index)
                .noalias() -=
                (B.transpose() * sigma_real - N_u.transpose() * rho_sr * b) * w;

            //
            // displacement equation, phasefield part
            //
            Kud.noalias() += B.transpose() * 2 * d_ip * sigma_tensile * N * w;

            double const d_dot_ip = N.dot(d_dot);

            if (history_variable_prev < strain_energy_tensile)
            {
                history_variable = strain_energy_tensile;
                Kdu.noalias() = Kud.transpose();
            }
            else
            {
                history_variable = history_variable_prev;
            }

            //
            // phasefield equation, phasefield part.
            //
            Kdd.noalias() += (Kdd_1 + N.transpose() * 2 * history_variable * N +
                              N.transpose() * 0.5 * gc / ls * N) *
                             w;
            double const M =
                _process_data.kinetic_coefficient(t, x_position)[0];
            local_rhs.template segment<phasefield_size>(phasefield_index)
                .noalias() -= (N.transpose() * d_dot_ip / M + Kdd_1 * d +
                               N.transpose() * d_ip * 2 * history_variable -
                               N.transpose() * 0.5 * gc / ls * (1 - d_ip)) *
                              w;

            Ddd.noalias() += N.transpose() / M * N * w;
        }
        // displacement equation, phasefield part
        local_Jac.template block<displacement_size, phasefield_size>(
                     displacement_index, phasefield_index)
            .noalias() += Kud;

        // phasefield equation, phasefield part.
        local_Jac.template block<phasefield_size, phasefield_size>(
                     phasefield_index, phasefield_index)
            .noalias() += Kdd + Ddd / dt;

        // phasefield equation, displacement part.
        local_Jac.template block<phasefield_size, displacement_size>(
                     phasefield_index, displacement_index)
            .noalias() += Kdu;

        // Eigen::EigenSolver<JacobianMatrix> eigensolver(local_Jac);
        // std::cout << "eigenvalues" << eigensolver.eigenvalues() << "\n";
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
            if (_ip_data[ip].history_variable_prev <
                _ip_data[ip].history_variable)
            {
                _ip_data[ip].history_variable_prev =
                    _ip_data[ip].history_variable;
            }
        }
    }

    Eigen::Map<const Eigen::RowVectorXd> getShapeMatrix(
        const unsigned integration_point) const override
    {
        auto const& N = _secondary_data.N[integration_point];

        // assumes N is stored contiguously in memory
        return Eigen::Map<const Eigen::RowVectorXd>(N.data(), N.size());
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

    std::vector<double> const& getIntPtSigmaYZ(
        std::vector<double>& cache) const override
    {
        assert(DisplacementDim == 3);
        return getIntPtSigma(cache, 4);
    }

    std::vector<double> const& getIntPtSigmaXZ(
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

    std::vector<double> const& getIntPtEpsilonYZ(
        std::vector<double>& cache) const override
    {
        assert(DisplacementDim == 3);
        return getIntPtEpsilon(cache, 4);
    }

    std::vector<double> const& getIntPtEpsilonXZ(
        std::vector<double>& cache) const override
    {
        assert(DisplacementDim == 3);
        return getIntPtEpsilon(cache, 5);
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
                cache.push_back(ip_data.sigma_real[component]);
            else  // mixed xy, yz, xz components
                cache.push_back(ip_data.sigma_real[component] / std::sqrt(2));
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
            cache.push_back(ip_data.eps[component]);
        }

        return cache;
    }

    PhaseFieldProcessData<DisplacementDim>& _process_data;

    std::vector<IntegrationPointData<BMatricesType, ShapeMatricesType,
                                     DisplacementDim>> _ip_data;

    IntegrationMethod _integration_method;
    MeshLib::Element const& _element;
    SecondaryData<typename ShapeMatrices::ShapeType> _secondary_data;

    static const int phasefield_index = 0;
    static const int phasefield_size = ShapeFunction::NPOINTS;
    static const int displacement_index = ShapeFunction::NPOINTS;
    static const int displacement_size =
        ShapeFunction::NPOINTS * DisplacementDim;
    static const int kelvin_vector_size =
        KelvinVectorDimensions<DisplacementDim>::value;
};

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim, int DisplacementDim>
class LocalAssemblerData final
    : public PhaseFieldLocalAssembler<ShapeFunction, IntegrationMethod,
                                      DisplacementDim>
{
public:
    LocalAssemblerData(LocalAssemblerData const&) = delete;
    LocalAssemblerData(LocalAssemblerData&&) = delete;

    LocalAssemblerData(MeshLib::Element const& e,
                       std::size_t const local_matrix_size,
                       bool is_axially_symmetric,
                       unsigned const integration_order,
                       PhaseFieldProcessData<DisplacementDim>& process_data)
        : PhaseFieldLocalAssembler<ShapeFunction, IntegrationMethod,
                                   DisplacementDim>(
              e, local_matrix_size, is_axially_symmetric, integration_order,
              process_data)
    {
    }
};

}  // namespace PhaseField
}  // namespace ProcessLib
