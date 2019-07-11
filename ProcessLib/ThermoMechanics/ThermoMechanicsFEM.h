/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>
#include <vector>

#include "MaterialLib/SolidModels/SelectSolidConstitutiveRelation.h"
#include "MathLib/KelvinVector.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "NumLib/Extrapolation/ExtrapolatableElement.h"
#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "ParameterLib/Parameter.h"
#include "ProcessLib/Deformation/BMatrixPolicy.h"
#include "ProcessLib/Deformation/LinearBMatrix.h"
#include "ProcessLib/Utils/InitShapeMatrices.h"

#include "LocalAssemblerInterface.h"
#include "ThermoMechanicsProcessData.h"

namespace ProcessLib
{
namespace ThermoMechanics
{
template <typename BMatricesType, typename ShapeMatricesType,
          int DisplacementDim>
struct IntegrationPointData final
{
    explicit IntegrationPointData(
        MaterialLib::Solids::MechanicsBase<DisplacementDim> const&
            solid_material)
        : solid_material(solid_material),
          material_state_variables(
              solid_material.createMaterialStateVariables())
    {
    }

    // total stress
    typename BMatricesType::KelvinVectorType sigma, sigma_prev;

    // total strain
    typename BMatricesType::KelvinVectorType eps, eps_prev;

    // mechanical strain
    typename BMatricesType::KelvinVectorType eps_m, eps_m_prev;

    // Non-equilibrium initial stress.
    typename BMatricesType::KelvinVectorType sigma_neq =
        BMatricesType::KelvinVectorType::Zero(
            MathLib::KelvinVector::KelvinVectorDimensions<
                DisplacementDim>::value);

    MaterialLib::Solids::MechanicsBase<DisplacementDim> const& solid_material;
    std::unique_ptr<typename MaterialLib::Solids::MechanicsBase<
        DisplacementDim>::MaterialStateVariables>
        material_state_variables;
    double solid_density;
    double solid_density_prev;

    double integration_weight;
    typename ShapeMatricesType::NodalRowVectorType N;
    typename ShapeMatricesType::GlobalDimNodalMatrixType dNdx;

    void pushBackState()
    {
        eps_prev = eps;
        eps_m_prev = eps_m;
        sigma_prev = sigma;
        solid_density_prev = solid_density;
        material_state_variables->pushBackState();
    }

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

/// Used for the extrapolation of the integration point values. It is ordered
/// (and stored) by integration points.
template <typename ShapeMatrixType>
struct SecondaryData
{
    std::vector<ShapeMatrixType, Eigen::aligned_allocator<ShapeMatrixType>> N;
};

template <typename ShapeFunction, typename IntegrationMethod,
          int DisplacementDim>
class ThermoMechanicsLocalAssembler
    : public ThermoMechanicsLocalAssemblerInterface
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

    ThermoMechanicsLocalAssembler(ThermoMechanicsLocalAssembler const&) =
        delete;
    ThermoMechanicsLocalAssembler(ThermoMechanicsLocalAssembler&&) = delete;

    ThermoMechanicsLocalAssembler(
        MeshLib::Element const& e,
        std::size_t const /*local_matrix_size*/,
        bool const is_axially_symmetric,
        unsigned const integration_order,
        ThermoMechanicsProcessData<DisplacementDim>& process_data)
        : _process_data(process_data),
          _integration_method(integration_order),
          _element(e),
          _is_axially_symmetric(is_axially_symmetric)
    {
        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        _ip_data.reserve(n_integration_points);
        _secondary_data.N.resize(n_integration_points);

        auto const shape_matrices =
            initShapeMatrices<ShapeFunction, ShapeMatricesType,
                              IntegrationMethod, DisplacementDim>(
                e, is_axially_symmetric, _integration_method);

        auto& solid_material =
            MaterialLib::Solids::selectSolidConstitutiveRelation(
                _process_data.solid_materials,
                _process_data.material_ids,
                e.getID());

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            _ip_data.emplace_back(solid_material);
            auto& ip_data = _ip_data[ip];
            ip_data.integration_weight =
                _integration_method.getWeightedPoint(ip).getWeight() *
                shape_matrices[ip].integralMeasure * shape_matrices[ip].detJ;

            static const int kelvin_vector_size =
                MathLib::KelvinVector::KelvinVectorDimensions<
                    DisplacementDim>::value;
            ip_data.sigma.setZero(kelvin_vector_size);
            ip_data.sigma_prev.setZero(kelvin_vector_size);
            ip_data.eps.setZero(kelvin_vector_size);
            ip_data.eps_prev.setZero(kelvin_vector_size);
            ip_data.eps_m.setZero(kelvin_vector_size);
            ip_data.eps_m_prev.setZero(kelvin_vector_size);

            ParameterLib::SpatialPosition x_position;
            x_position.setElementID(_element.getID());
            ip_data.solid_density =
                _process_data.reference_solid_density(0, x_position)[0];
            ip_data.solid_density_prev = ip_data.solid_density;

            ip_data.N = shape_matrices[ip].N;
            ip_data.dNdx = shape_matrices[ip].dNdx;

            _secondary_data.N[ip] = shape_matrices[ip].N;

            // Initialize current time step values
            if (_process_data.nonequilibrium_stress)
            {
                // Computation of non-equilibrium stress.
                x_position.setCoordinates(MathLib::Point3d(
                    interpolateCoordinates<ShapeFunction, ShapeMatricesType>(
                        e, ip_data.N)));
                std::vector<double> sigma_neq_data =
                    (*_process_data.nonequilibrium_stress)(
                        std::numeric_limits<
                            double>::quiet_NaN() /* time independent */,
                        x_position);
                ip_data.sigma_neq =
                    Eigen::Map<typename BMatricesType::KelvinVectorType>(
                        sigma_neq_data.data(),
                        MathLib::KelvinVector::KelvinVectorDimensions<
                            DisplacementDim>::value,
                        1);
            }
            // Initialization from non-equilibrium sigma, which is zero by
            // default, or is set to some value.
            ip_data.sigma = ip_data.sigma_neq;
        }
    }

    /// Returns number of read integration points.
    std::size_t setIPDataInitialConditions(std::string const& name,
                                           double const* values,
                                           int const integration_order) override
    {
        if (integration_order !=
            static_cast<int>(_integration_method.getIntegrationOrder()))
        {
            OGS_FATAL(
                "Setting integration point initial conditions; The integration "
                "order of the local assembler for element %d is different from "
                "the integration order in the initial condition.",
                _element.getID());
        }

        if (name == "sigma_ip")
        {
            return setSigma(values);
        }
        if (name == "epsilon_ip")
        {
            return setEpsilon(values);
        }
        if (name == "epsilon_m_ip")
        {
            return setEpsilonMechanical(values);
        }

        return 0;
    }

    void assemble(double const /*t*/, std::vector<double> const& /*local_x*/,
                  std::vector<double>& /*local_M_data*/,
                  std::vector<double>& /*local_K_data*/,
                  std::vector<double>& /*local_rhs_data*/) override
    {
        OGS_FATAL(
            "ThermoMechanicsLocalAssembler: assembly without jacobian is not "
            "implemented.");
    }

    void assembleWithJacobianForStaggeredScheme(
        double const t, std::vector<double> const& local_xdot,
        const double dxdot_dx, const double dx_dx,
        std::vector<double>& local_M_data, std::vector<double>& local_K_data,
        std::vector<double>& local_b_data, std::vector<double>& local_Jac_data,
        LocalCoupledSolutions const& local_coupled_solutions) override;

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
        assert(local_matrix_size == temperature_size + displacement_size);

        auto T = Eigen::Map<typename ShapeMatricesType::template VectorType<
            temperature_size> const>(local_x.data() + temperature_index,
                                     temperature_size);

        auto u = Eigen::Map<typename ShapeMatricesType::template VectorType<
            displacement_size> const>(local_x.data() + displacement_index,
                                      displacement_size);

        auto T_dot = Eigen::Map<typename ShapeMatricesType::template VectorType<
            temperature_size> const>(local_xdot.data() + temperature_index,
                                     temperature_size);

        auto local_Jac = MathLib::createZeroedMatrix<JacobianMatrix>(
            local_Jac_data, local_matrix_size, local_matrix_size);

        auto local_rhs = MathLib::createZeroedVector<RhsVector>(
            local_rhs_data, local_matrix_size);

        typename ShapeMatricesType::template MatrixType<displacement_size,
                                                        temperature_size>
            KuT;
        KuT.setZero(displacement_size, temperature_size);

        typename ShapeMatricesType::NodalMatrixType KTT;
        KTT.setZero(temperature_size, temperature_size);

        typename ShapeMatricesType::NodalMatrixType DTT;
        DTT.setZero(temperature_size, temperature_size);

        double const& dt = _process_data.dt;

        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        ParameterLib::SpatialPosition x_position;
        x_position.setElementID(_element.getID());

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            x_position.setIntegrationPoint(ip);
            auto const& w = _ip_data[ip].integration_weight;
            auto const& N = _ip_data[ip].N;
            auto const& dNdx = _ip_data[ip].dNdx;

            auto const x_coord =
                interpolateXCoordinate<ShapeFunction, ShapeMatricesType>(
                    _element, N);
            auto const& B = LinearBMatrix::computeBMatrix<
                DisplacementDim, ShapeFunction::NPOINTS,
                typename BMatricesType::BMatrixType>(dNdx, N, x_coord,
                                                     _is_axially_symmetric);

            auto& sigma = _ip_data[ip].sigma;
            auto const& sigma_prev = _ip_data[ip].sigma_prev;
            auto const& sigma_neq = _ip_data[ip].sigma_neq;

            auto& eps = _ip_data[ip].eps;
            auto const& eps_prev = _ip_data[ip].eps_prev;

            auto& eps_m = _ip_data[ip].eps_m;
            auto const& eps_m_prev = _ip_data[ip].eps_m_prev;

            auto& state = _ip_data[ip].material_state_variables;

            double const dT = N.dot(T_dot) * dt;
            // calculate thermally induced strain
            // assume isotropic thermal expansion
            auto const alpha =
                _process_data.linear_thermal_expansion_coefficient(
                    t, x_position)[0];
            double const linear_thermal_strain_increment = alpha * dT;

            //
            // displacement equation, displacement part
            //
            eps.noalias() = B * u;

            using Invariants = MathLib::KelvinVector::Invariants<
                MathLib::KelvinVector::KelvinVectorDimensions<
                    DisplacementDim>::value>;

            // assume isotropic thermal expansion
            const double T_ip = N.dot(T);  // T at integration point
            eps_m.noalias() =
                eps_m_prev + eps - eps_prev -
                linear_thermal_strain_increment * Invariants::identity2;

            auto&& solution = _ip_data[ip].solid_material.integrateStress(
                t, x_position, dt, eps_m_prev, eps_m, sigma_prev, *state, T_ip);

            if (!solution)
            {
                OGS_FATAL("Computation of local constitutive relation failed.");
            }

            MathLib::KelvinVector::KelvinMatrixType<DisplacementDim> C;
            std::tie(sigma, state, C) = std::move(*solution);

            local_Jac
                .template block<displacement_size, displacement_size>(
                    displacement_index, displacement_index)
                .noalias() += B.transpose() * C * B * w;

            typename ShapeMatricesType::template MatrixType<DisplacementDim,
                                                            displacement_size>
                N_u = ShapeMatricesType::template MatrixType<
                    DisplacementDim,
                    displacement_size>::Zero(DisplacementDim,
                                             displacement_size);

            for (int i = 0; i < DisplacementDim; ++i)
            {
                N_u.template block<1, displacement_size / DisplacementDim>(
                       i, i * displacement_size / DisplacementDim)
                    .noalias() = N;
            }

            // calculate real density
            // rho_s_{n+1} * (V_{n} + dV) = rho_s_n * V_n
            // dV = 3 * alpha * dT * V_0
            // rho_s_{n+1} = rho_s_n / (1 + 3 * alpha * dT )
            // see reference solid density description for details.
            auto& rho_s = _ip_data[ip].solid_density;
            rho_s = _ip_data[ip].solid_density_prev /
                    (1 + 3 * linear_thermal_strain_increment);

            auto const& b = _process_data.specific_body_force;
            local_rhs
                .template block<displacement_size, 1>(displacement_index, 0)
                .noalias() -= (B.transpose() * (sigma - sigma_neq) -
                               N_u.transpose() * rho_s * b) *
                              w;

            //
            // displacement equation, temperature part
            //
            KuT.noalias() +=
                B.transpose() * C * alpha * Invariants::identity2 * N * w;
            if (_ip_data[ip].solid_material.getConstitutiveModel() ==
                MaterialLib::Solids::ConstitutiveModel::CreepBGRa)
            {
                auto const s =
                    Invariants::deviatoric_projection * (sigma - sigma_neq);
                double const norm_s = Invariants::FrobeniusNorm(s);
                const double creep_coefficient =
                    _ip_data[ip]
                        .solid_material.getTemperatureRelatedCoefficient(
                            t, dt, x_position, T_ip, norm_s);
                KuT.noalias() += creep_coefficient * B.transpose() * s * N * w;
            }

            //
            // temperature equation, temperature part;
            //
            auto const lambda =
                _process_data.thermal_conductivity(t, x_position)[0];
            KTT.noalias() += dNdx.transpose() * lambda * dNdx * w;

            auto const c =
                _process_data.specific_heat_capacity(t, x_position)[0];
            DTT.noalias() += N.transpose() * rho_s * c * N * w;
        }

        // temperature equation, temperature part
        local_Jac
            .template block<temperature_size, temperature_size>(
                temperature_index, temperature_index)
            .noalias() += KTT + DTT / dt;

        // displacement equation, temperature part
        local_Jac
            .template block<displacement_size, temperature_size>(
                displacement_index, temperature_index)
            .noalias() -= KuT;

        local_rhs.template block<temperature_size, 1>(temperature_index, 0)
            .noalias() -= KTT * T + DTT * T_dot;
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
        auto const& N = _secondary_data.N[integration_point];

        // assumes N is stored contiguously in memory
        return Eigen::Map<const Eigen::RowVectorXd>(N.data(), N.size());
    }

private:
    /**
     * Assemble local matrices and vectors arise from the linearized discretized
     * weak form of the residual of the momentum balance equation,
     *      \f[
     *            \nabla (\sigma - \mathbf{D} \alpha_T (T-T_0) \mathrm{I}) = f
     *      \f]
     * where \f$ \sigma\f$ is the effective stress tensor, \f$\mathbf{D}\f$ is
     * the tangential operator, \f$T\f$ is the  temperature, \f$T_0\f$ is the
     * initial temperature, \f$\alpha_T\f$ is the linear thermal expansion,
     * \f$\mathrm{I}\f$ is the identity tensor, and \f$f\f$ is the body force.
     *
     * @param t               Time
     * @param local_xdot      Nodal values of \f$\dot{x}\f$ of an element.
     * @param dxdot_dx        Value of \f$\dot{x} \cdot dx\f$.
     * @param dx_dx           Value of \f$ x \cdot dx\f$.
     * @param local_M_data    Mass matrix of an element, which takes the form of
     *                        \f$ \int N^T N\mathrm{d}\Omega\f$. Not used.
     * @param local_K_data    Laplacian matrix of an element, which takes the
     *         form of \f$ \int (\nabla N)^T K \nabla N\mathrm{d}\Omega\f$.
     *                        Not used.
     * @param local_b_data    Right hand side vector of an element.
     * @param local_Jac_data  Element Jacobian matrix for the Newton-Raphson
     *                        method.
     * @param local_coupled_solutions Nodal values of solutions of the coupled
     *                                temperature of an element.
     */

    void assembleWithJacobianForDeformationEquations(
        double const t, std::vector<double> const& local_xdot,
        const double dxdot_dx, const double dx_dx,
        std::vector<double>& local_M_data, std::vector<double>& local_K_data,
        std::vector<double>& local_b_data, std::vector<double>& local_Jac_data,
        LocalCoupledSolutions const& local_coupled_solutions);

    /**
     * Assemble local matrices and vectors arise from the linearized discretized
     * weak form of the residual of the energy balance equation,
     *      \f[
     *          \rho c_p \cdot{T} - \nabla (\mathbf{K} (\nabla T) = Q_T
     *      \f]
     * where \f$ rho\f$ is the solid density, \f$ c_p\f$ is the spefific heat
     * capacity, \f$ \mathbf{K} \f$ is the thermal conductivity, and \f$ Q_T\f$
     * is the source/sink term.
     *
     * @param t               Time
     * @param local_xdot      Nodal values of \f$\dot{x}\f$ of an element.
     * @param dxdot_dx        Value of \f$\dot{x} \cdot dx\f$.
     * @param dx_dx           Value of \f$ x \cdot dx\f$.
     * @param local_M_data    Mass matrix of an element, which takes the form of
     *                        \f$ \int N^T N\mathrm{d}\Omega\f$. Not used.
     * @param local_K_data    Laplacian matrix of an element, which takes the
     *         form of \f$ \int (\nabla N)^T K \nabla N\mathrm{d}\Omega\f$.
     *                        Not used.
     * @param local_b_data    Right hand side vector of an element.
     * @param local_Jac_data  Element Jacobian matrix for the Newton-Raphson
     *                        method.
     * @param local_coupled_solutions Nodal values of solutions of the coupled
     *                                displacement of an element.
     */
    void assembleWithJacobianForHeatConductionEquations(
        double const t, std::vector<double> const& local_xdot,
        const double dxdot_dx, const double dx_dx,
        std::vector<double>& local_M_data, std::vector<double>& local_K_data,
        std::vector<double>& local_b_data, std::vector<double>& local_Jac_data,
        LocalCoupledSolutions const& local_coupled_solutions);

    std::size_t setSigma(double const* values)
    {
        auto const kelvin_vector_size =
            MathLib::KelvinVector::KelvinVectorDimensions<
                DisplacementDim>::value;
        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        auto sigma_values =
            Eigen::Map<Eigen::Matrix<double, kelvin_vector_size, Eigen::Dynamic,
                                     Eigen::ColMajor> const>(
                values, kelvin_vector_size, n_integration_points);

        for (unsigned ip = 0; ip < n_integration_points; ++ip)
        {
            _ip_data[ip].sigma =
                MathLib::KelvinVector::symmetricTensorToKelvinVector(
                    sigma_values.col(ip));
        }

        return n_integration_points;
    }

    // TODO (naumov) This method is same as getIntPtSigma but for arguments and
    // the ordering of the cache_mat.
    // There should be only one.
    std::vector<double> getSigma() const override
    {
        auto const kelvin_vector_size =
            MathLib::KelvinVector::KelvinVectorDimensions<
                DisplacementDim>::value;
        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        std::vector<double> ip_sigma_values;
        auto cache_mat = MathLib::createZeroedMatrix<Eigen::Matrix<
            double, Eigen::Dynamic, kelvin_vector_size, Eigen::RowMajor>>(
            ip_sigma_values, n_integration_points, kelvin_vector_size);

        for (unsigned ip = 0; ip < n_integration_points; ++ip)
        {
            auto const& sigma = _ip_data[ip].sigma;
            cache_mat.row(ip) =
                MathLib::KelvinVector::kelvinVectorToSymmetricTensor(sigma);
        }

        return ip_sigma_values;
    }

    std::vector<double> const& getIntPtSigma(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        static const int kelvin_vector_size =
            MathLib::KelvinVector::KelvinVectorDimensions<
                DisplacementDim>::value;
        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        cache.clear();
        auto cache_mat = MathLib::createZeroedMatrix<Eigen::Matrix<
            double, kelvin_vector_size, Eigen::Dynamic, Eigen::RowMajor>>(
            cache, kelvin_vector_size, n_integration_points);

        for (unsigned ip = 0; ip < n_integration_points; ++ip)
        {
            auto const& sigma = _ip_data[ip].sigma;
            cache_mat.col(ip) =
                MathLib::KelvinVector::kelvinVectorToSymmetricTensor(sigma);
        }

        return cache;
    }

    std::size_t setEpsilon(double const* values)
    {
        auto const kelvin_vector_size =
            MathLib::KelvinVector::KelvinVectorDimensions<
                DisplacementDim>::value;
        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        auto epsilon_values =
            Eigen::Map<Eigen::Matrix<double, kelvin_vector_size, Eigen::Dynamic,
                                     Eigen::ColMajor> const>(
                values, kelvin_vector_size, n_integration_points);

        for (unsigned ip = 0; ip < n_integration_points; ++ip)
        {
            _ip_data[ip].eps =
                MathLib::KelvinVector::symmetricTensorToKelvinVector(
                    epsilon_values.col(ip));
        }

        return n_integration_points;
    }

    std::vector<double> getEpsilon() const override
    {
        auto const kelvin_vector_size =
            MathLib::KelvinVector::KelvinVectorDimensions<
                DisplacementDim>::value;
        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        std::vector<double> ip_epsilon_values;
        auto cache_mat = MathLib::createZeroedMatrix<Eigen::Matrix<
            double, Eigen::Dynamic, kelvin_vector_size, Eigen::RowMajor>>(
            ip_epsilon_values, n_integration_points, kelvin_vector_size);

        for (unsigned ip = 0; ip < n_integration_points; ++ip)
        {
            auto const& eps = _ip_data[ip].eps;
            cache_mat.row(ip) =
                MathLib::KelvinVector::kelvinVectorToSymmetricTensor(eps);
        }

        return ip_epsilon_values;
    }

    std::vector<double> const& getIntPtEpsilon(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        auto const kelvin_vector_size =
            MathLib::KelvinVector::KelvinVectorDimensions<
                DisplacementDim>::value;
        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        cache.clear();
        auto cache_mat = MathLib::createZeroedMatrix<Eigen::Matrix<
            double, kelvin_vector_size, Eigen::Dynamic, Eigen::RowMajor>>(
            cache, kelvin_vector_size, n_integration_points);

        for (unsigned ip = 0; ip < n_integration_points; ++ip)
        {
            auto const& eps = _ip_data[ip].eps;
            cache_mat.col(ip) =
                MathLib::KelvinVector::kelvinVectorToSymmetricTensor(eps);
        }

        return cache;
    }

    std::size_t setEpsilonMechanical(double const* values)
    {
        auto const kelvin_vector_size =
            MathLib::KelvinVector::KelvinVectorDimensions<
                DisplacementDim>::value;
        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        auto epsilon_m_values =
            Eigen::Map<Eigen::Matrix<double, kelvin_vector_size, Eigen::Dynamic,
                                     Eigen::ColMajor> const>(
                values, kelvin_vector_size, n_integration_points);

        for (unsigned ip = 0; ip < n_integration_points; ++ip)
        {
            _ip_data[ip].eps_m =
                MathLib::KelvinVector::symmetricTensorToKelvinVector(
                    epsilon_m_values.col(ip));
        }

        return n_integration_points;
    }

    std::vector<double> getEpsilonMechanical() const override
    {
        auto const kelvin_vector_size =
            MathLib::KelvinVector::KelvinVectorDimensions<
                DisplacementDim>::value;
        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        std::vector<double> ip_epsilon_m_values;
        auto cache_mat = MathLib::createZeroedMatrix<Eigen::Matrix<
            double, Eigen::Dynamic, kelvin_vector_size, Eigen::RowMajor>>(
            ip_epsilon_m_values, n_integration_points, kelvin_vector_size);

        for (unsigned ip = 0; ip < n_integration_points; ++ip)
        {
            auto const& eps_m = _ip_data[ip].eps_m;
            cache_mat.row(ip) =
                MathLib::KelvinVector::kelvinVectorToSymmetricTensor(eps_m);
        }

        return ip_epsilon_m_values;
    }

    ThermoMechanicsProcessData<DisplacementDim>& _process_data;

    std::vector<
        IntegrationPointData<BMatricesType, ShapeMatricesType, DisplacementDim>,
        Eigen::aligned_allocator<IntegrationPointData<
            BMatricesType, ShapeMatricesType, DisplacementDim>>>
        _ip_data;

    IntegrationMethod _integration_method;
    MeshLib::Element const& _element;
    SecondaryData<typename ShapeMatrices::ShapeType> _secondary_data;
    bool const _is_axially_symmetric;

    static const int temperature_index = 0;
    static const int temperature_size = ShapeFunction::NPOINTS;
    static const int displacement_index = ShapeFunction::NPOINTS;
    static const int displacement_size =
        ShapeFunction::NPOINTS * DisplacementDim;
};

}  // namespace ThermoMechanics
}  // namespace ProcessLib

#include "ThermoMechanicsFEM-impl.h"
