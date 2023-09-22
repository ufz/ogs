/**
 * \file
 * \copyright
 * Copyright (c) 2012-2023, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <Eigen/Eigenvalues>
#include <limits>
#include <memory>
#include <vector>

#include "LocalAssemblerInterface.h"
#include "MaterialLib/PhysicalConstant.h"
#include "MathLib/EigenBlockMatrixView.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "NumLib/DOF/LocalDOF.h"
#include "NumLib/Extrapolation/ExtrapolatableElement.h"
#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/InitShapeMatrices.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "ParameterLib/Parameter.h"
#include "ProcessLib/Deformation/BMatrixPolicy.h"
#include "ProcessLib/Deformation/GMatrixPolicy.h"
#include "ProcessLib/Deformation/LinearBMatrix.h"
#include "ProcessLib/LocalAssemblerInterface.h"
#include "ProcessLib/LocalAssemblerTraits.h"

namespace ProcessLib
{
namespace SmallDeformation
{
namespace MPL = MaterialPropertyLib;

template <typename BMatricesType, typename ShapeMatricesType,
          int DisplacementDim>
struct IntegrationPointData final
{
    typename BMatricesType::KelvinVectorType sigma, sigma_prev;
    typename BMatricesType::KelvinVectorType eps;
    double free_energy_density = 0;

    double integration_weight;
    typename ShapeMatricesType::NodalRowVectorType N;
    typename ShapeMatricesType::GlobalDimNodalMatrixType dNdx;

    void pushBackState() { sigma_prev = sigma; }

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

/// Used for the extrapolation of the integration point values. It is ordered
/// (and stored) by integration points.
template <typename ShapeMatrixType>
struct SecondaryData
{
    std::vector<ShapeMatrixType, Eigen::aligned_allocator<ShapeMatrixType>> N;
};

template <typename ShapeFunction, int DisplacementDim>
class SmallDeformationLocalAssembler
    : public SmallDeformationLocalAssemblerInterface<DisplacementDim>
{
public:
    using ShapeMatricesType =
        ShapeMatrixPolicyType<ShapeFunction, DisplacementDim>;
    using NodalMatrixType = typename ShapeMatricesType::NodalMatrixType;
    using NodalVectorType = typename ShapeMatricesType::NodalVectorType;
    using ShapeMatrices = typename ShapeMatricesType::ShapeMatrices;
    using BMatricesType = BMatrixPolicyType<ShapeFunction, DisplacementDim>;

    using BMatrixType = typename BMatricesType::BMatrixType;
    using StiffnessMatrixType = typename BMatricesType::StiffnessMatrixType;
    using NodalForceVectorType = typename BMatricesType::NodalForceVectorType;
    using NodalDisplacementVectorType =
        typename BMatricesType::NodalForceVectorType;

    using GMatricesType = GMatrixPolicyType<ShapeFunction, DisplacementDim>;
    using GradientVectorType = typename GMatricesType::GradientVectorType;
    using GradientMatrixType = typename GMatricesType::GradientMatrixType;
    using IpData =
        IntegrationPointData<BMatricesType, ShapeMatricesType, DisplacementDim>;

    static constexpr auto& N_u_op = MathLib::eigenBlockMatrixView<
        DisplacementDim, typename ShapeMatricesType::NodalRowVectorType>;

    SmallDeformationLocalAssembler(SmallDeformationLocalAssembler const&) =
        delete;
    SmallDeformationLocalAssembler(SmallDeformationLocalAssembler&&) = delete;

    SmallDeformationLocalAssembler(
        MeshLib::Element const& e,
        std::size_t const /*local_matrix_size*/,
        NumLib::GenericIntegrationMethod const& integration_method,
        bool const is_axially_symmetric,
        SmallDeformationProcessData<DisplacementDim>& process_data)
        : SmallDeformationLocalAssemblerInterface<DisplacementDim>(
              e, integration_method, is_axially_symmetric, process_data)
    {
        unsigned const n_integration_points =
            this->integration_method_.getNumberOfPoints();

        _ip_data.resize(n_integration_points);
        _secondary_data.N.resize(n_integration_points);

        auto const shape_matrices =
            NumLib::initShapeMatrices<ShapeFunction, ShapeMatricesType,
                                      DisplacementDim>(
                e, is_axially_symmetric, this->integration_method_);

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            auto& ip_data = _ip_data[ip];
            auto const& sm = shape_matrices[ip];
            _ip_data[ip].integration_weight =
                this->integration_method_.getWeightedPoint(ip).getWeight() *
                sm.integralMeasure * sm.detJ;

            ip_data.N = sm.N;
            ip_data.dNdx = sm.dNdx;

            static const int kelvin_vector_size =
                MathLib::KelvinVector::kelvin_vector_dimensions(
                    DisplacementDim);
            // Initialize current time step values
            ip_data.sigma.setZero(kelvin_vector_size);
            ip_data.eps.setZero(kelvin_vector_size);

            // Previous time step values are not initialized and are set later.
            ip_data.sigma_prev.resize(kelvin_vector_size);

            _secondary_data.N[ip] = shape_matrices[ip].N;
        }
    }

    /// Returns number of read integration points.
    std::size_t setIPDataInitialConditions(std::string const& name,
                                           double const* values,
                                           int const integration_order) override
    {
        if (integration_order !=
            static_cast<int>(this->integration_method_.getIntegrationOrder()))
        {
            OGS_FATAL(
                "Setting integration point initial conditions; The integration "
                "order of the local assembler for element {:d} is different "
                "from the integration order in the initial condition.",
                this->element_.getID());
        }

        if (name == "sigma_ip")
        {
            if (this->process_data_.initial_stress != nullptr)
            {
                OGS_FATAL(
                    "Setting initial conditions for stress from integration "
                    "point data and from a parameter '{:s}' is not possible "
                    "simultaneously.",
                    this->process_data_.initial_stress->name);
            }
            return setSigma(values);
        }
        if (name.starts_with("material_state_variable_") &&
            name.ends_with("_ip"))
        {
            std::string const variable_name =
                name.substr(24, name.size() - 24 - 3);
            DBUG("Setting material state variable '{:s}'", variable_name);

            auto const& internal_variables =
                this->solid_material_.getInternalVariables();
            if (auto const iv = std::find_if(
                    begin(internal_variables), end(internal_variables),
                    [&variable_name](auto const& iv)
                    { return iv.name == variable_name; });
                iv != end(internal_variables))
            {
                return ProcessLib::
                    setIntegrationPointDataMaterialStateVariables(
                        values, this->material_states_,
                        &MaterialStateData<
                            DisplacementDim>::material_state_variables,
                        iv->reference);
            }

            WARN(
                "Could not find variable {:s} in solid material model's "
                "internal variables.",
                variable_name);
        }

        return 0;
    }

    void initializeConcrete() override
    {
        unsigned const n_integration_points =
            this->integration_method_.getNumberOfPoints();
        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            auto& ip_data = _ip_data[ip];

            ParameterLib::SpatialPosition const x_position{
                std::nullopt, this->element_.getID(), ip,
                MathLib::Point3d(
                    NumLib::interpolateCoordinates<ShapeFunction,
                                                   ShapeMatricesType>(
                        this->element_, ip_data.N))};

            /// Set initial stress from parameter.
            if (this->process_data_.initial_stress != nullptr)
            {
                ip_data.sigma =
                    MathLib::KelvinVector::symmetricTensorToKelvinVector<
                        DisplacementDim>((*this->process_data_.initial_stress)(
                        std::numeric_limits<
                            double>::quiet_NaN() /* time independent */,
                        x_position));
            }

            double const t = 0;  // TODO (naumov) pass t from top
            auto& material_state = this->material_states_[ip];
            this->solid_material_.initializeInternalStateVariables(
                t, x_position, *material_state.material_state_variables);

            ip_data.pushBackState();
            material_state.pushBackState();
        }
    }

    /// Updates sigma, eps, and state through passed integration point data.
    /// \returns tangent stiffness and density.
    std::tuple<MathLib::KelvinVector::KelvinMatrixType<DisplacementDim>, double>
    updateConstitutiveRelations(
        Eigen::Ref<Eigen::VectorXd const> const& u,
        Eigen::Ref<Eigen::VectorXd const> const& u_prev,
        ParameterLib::SpatialPosition const& x_position, double const t,
        double const dt,
        IntegrationPointData<BMatricesType, ShapeMatricesType, DisplacementDim>&
            ip_data,
        MaterialStateData<DisplacementDim>& material_state) const
    {
        auto const& solid_phase =
            this->process_data_.media_map.getMedium(this->element_.getID())
                ->phase("Solid");

        MPL::VariableArray variables_prev;
        MPL::VariableArray variables;

        auto const& sigma_prev = ip_data.sigma_prev;

        auto& eps = ip_data.eps;
        auto& sigma = ip_data.sigma;

        auto const& N = ip_data.N;
        auto const& dNdx = ip_data.dNdx;
        auto const x_coord =
            NumLib::interpolateXCoordinate<ShapeFunction, ShapeMatricesType>(
                this->element_, N);
        auto const B =
            LinearBMatrix::computeBMatrix<DisplacementDim,
                                          ShapeFunction::NPOINTS,
                                          typename BMatricesType::BMatrixType>(
                dNdx, N, x_coord, this->is_axially_symmetric_);

        eps.noalias() = B * u;

        variables_prev.stress
            .emplace<MathLib::KelvinVector::KelvinVectorType<DisplacementDim>>(
                sigma_prev);
        variables_prev.mechanical_strain
            .emplace<MathLib::KelvinVector::KelvinVectorType<DisplacementDim>>(
                B * u_prev);

        double const T_ref =
            this->process_data_.reference_temperature
                ? (*this->process_data_.reference_temperature)(t, x_position)[0]
                : std::numeric_limits<double>::quiet_NaN();

        variables_prev.temperature = T_ref;
        variables.mechanical_strain
            .emplace<MathLib::KelvinVector::KelvinVectorType<DisplacementDim>>(
                eps);
        variables.temperature = T_ref;

        auto const rho =
            solid_phase[MPL::PropertyType::density].template value<double>(
                variables, x_position, t, dt);

        auto&& solution = this->solid_material_.integrateStress(
            variables_prev, variables, t, x_position, dt,
            *material_state.material_state_variables);

        if (!solution)
        {
            OGS_FATAL("Computation of local constitutive relation failed.");
        }

        MathLib::KelvinVector::KelvinMatrixType<DisplacementDim> C;
        std::tie(sigma, material_state.material_state_variables, C) =
            std::move(*solution);

        return {C, rho};
    }

    void assemble(double const /*t*/, double const /*dt*/,
                  std::vector<double> const& /*local_x*/,
                  std::vector<double> const& /*local_x_prev*/,
                  std::vector<double>& /*local_M_data*/,
                  std::vector<double>& /*local_K_data*/,
                  std::vector<double>& /*local_b_data*/) override
    {
        OGS_FATAL(
            "SmallDeformationLocalAssembler: assembly without jacobian is not "
            "implemented.");
    }

    void assembleWithJacobian(double const t, double const dt,
                              std::vector<double> const& local_x,
                              std::vector<double> const& local_x_prev,
                              std::vector<double>& /*local_M_data*/,
                              std::vector<double>& /*local_K_data*/,
                              std::vector<double>& local_b_data,
                              std::vector<double>& local_Jac_data) override
    {
        auto const local_matrix_size = local_x.size();

        auto local_Jac = MathLib::createZeroedMatrix<StiffnessMatrixType>(
            local_Jac_data, local_matrix_size, local_matrix_size);

        auto local_b = MathLib::createZeroedVector<NodalDisplacementVectorType>(
            local_b_data, local_matrix_size);

        auto [u] = localDOF(local_x);
        auto [u_prev] = localDOF(local_x_prev);

        unsigned const n_integration_points =
            this->integration_method_.getNumberOfPoints();

        ParameterLib::SpatialPosition x_position;
        x_position.setElementID(this->element_.getID());

        auto const& b = this->process_data_.specific_body_force;

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            x_position.setIntegrationPoint(ip);
            auto const& w = _ip_data[ip].integration_weight;
            auto const& N = _ip_data[ip].N;
            auto const& dNdx = _ip_data[ip].dNdx;

            auto const x_coord =
                NumLib::interpolateXCoordinate<ShapeFunction,
                                               ShapeMatricesType>(
                    this->element_, N);
            auto const B = LinearBMatrix::computeBMatrix<
                DisplacementDim, ShapeFunction::NPOINTS,
                typename BMatricesType::BMatrixType>(
                dNdx, N, x_coord, this->is_axially_symmetric_);

            auto const& sigma = _ip_data[ip].sigma;

            auto const [C, rho] = updateConstitutiveRelations(
                u, u_prev, x_position, t, dt, _ip_data[ip],
                this->material_states_[ip]);

            local_b.noalias() -=
                (B.transpose() * sigma - N_u_op(N).transpose() * rho * b) * w;
            local_Jac.noalias() += B.transpose() * C * B * w;
        }
    }

    void postTimestepConcrete(Eigen::VectorXd const& local_x,
                              Eigen::VectorXd const& local_x_prev,
                              double const t, double const dt,
                              bool const /*use_monolithic_scheme*/,
                              int const /*process_id*/) override
    {
        unsigned const n_integration_points =
            this->integration_method_.getNumberOfPoints();

        ParameterLib::SpatialPosition x_position;
        x_position.setElementID(this->element_.getID());

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            x_position.setIntegrationPoint(ip);

            updateConstitutiveRelations(local_x, local_x_prev, x_position, t,
                                        dt, _ip_data[ip],
                                        this->material_states_[ip]);

            auto& eps = _ip_data[ip].eps;
            auto& sigma = _ip_data[ip].sigma;
            auto& state = *this->material_states_[ip].material_state_variables;

            // Update free energy density needed for material forces.
            _ip_data[ip].free_energy_density =
                this->solid_material_.computeFreeEnergyDensity(
                    t, x_position, dt, eps, sigma, state);

            _ip_data[ip].pushBackState();
            this->material_states_[ip].pushBackState();
        }
    }

    std::vector<double> const& getMaterialForces(
        std::vector<double> const& local_x,
        std::vector<double>& nodal_values) override
    {
        return ProcessLib::SmallDeformation::getMaterialForces<
            DisplacementDim, ShapeFunction, ShapeMatricesType,
            typename BMatricesType::NodalForceVectorType,
            NodalDisplacementVectorType, GradientVectorType,
            GradientMatrixType>(local_x, nodal_values,
                                this->integration_method_, _ip_data,
                                this->element_, this->is_axially_symmetric_);
    }

    Eigen::Map<const Eigen::RowVectorXd> getShapeMatrix(
        const unsigned integration_point) const override
    {
        auto const& N = _secondary_data.N[integration_point];

        // assumes N is stored contiguously in memory
        return Eigen::Map<const Eigen::RowVectorXd>(N.data(), N.size());
    }

    std::vector<double> const& getIntPtFreeEnergyDensity(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        cache.clear();
        cache.reserve(_ip_data.size());

        transform(
            cbegin(_ip_data), cend(_ip_data), back_inserter(cache),
            [](auto const& ip_data) { return ip_data.free_energy_density; });

        return cache;
    }

    std::size_t setSigma(double const* values)
    {
        return ProcessLib::setIntegrationPointKelvinVectorData<DisplacementDim>(
            values, _ip_data, &IpData::sigma);
    }

    // TODO (naumov) This method is same as getIntPtSigma but for arguments and
    // the ordering of the cache_mat.
    // There should be only one.
    std::vector<double> getSigma() const override
    {
        return ProcessLib::getIntegrationPointKelvinVectorData<DisplacementDim>(
            _ip_data, &IpData::sigma);
    }

    std::vector<double> const& getIntPtSigma(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        return ProcessLib::getIntegrationPointKelvinVectorData<DisplacementDim>(
            _ip_data, &IpData::sigma, cache);
    }

    std::vector<double> const& getIntPtEpsilon(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        return ProcessLib::getIntegrationPointKelvinVectorData<DisplacementDim>(
            _ip_data, &IpData::eps, cache);
    }

    void computeSecondaryVariableConcrete(
        double const /*t*/, double const /*dt*/, Eigen::VectorXd const& /*x*/,
        Eigen::VectorXd const& /*x_prev*/) override
    {
        int const elem_id = this->element_.getID();
        ParameterLib::SpatialPosition x_position;
        x_position.setElementID(elem_id);
        unsigned const n_integration_points =
            this->integration_method_.getNumberOfPoints();

        auto sigma_sum = MathLib::KelvinVector::tensorToKelvin<DisplacementDim>(
            Eigen::Matrix<double, 3, 3>::Zero());

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            x_position.setIntegrationPoint(ip);
            auto const& sigma = _ip_data[ip].sigma;
            sigma_sum += sigma;
        }

        Eigen::Matrix<double, 3, 3, 0, 3, 3> const sigma_avg =
            MathLib::KelvinVector::kelvinVectorToTensor(sigma_sum) /
            n_integration_points;

        Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 3, 3>> e_s(
            sigma_avg);

        Eigen::Map<Eigen::Vector3d>(
            &(*this->process_data_.principal_stress_values)[elem_id * 3], 3) =
            e_s.eigenvalues();

        auto eigen_vectors = e_s.eigenvectors();

        for (auto i = 0; i < 3; i++)
        {
            Eigen::Map<Eigen::Vector3d>(
                &(*this->process_data_.principal_stress_vector[i])[elem_id * 3],
                3) = eigen_vectors.col(i);
        }
    }

private:
    static constexpr auto localDOF(std::vector<double> const& x)
    {
        return NumLib::localDOF<
            NumLib::Vectorial<ShapeFunction, DisplacementDim>>(x);
    }

private:
    std::vector<IpData, Eigen::aligned_allocator<IpData>> _ip_data;

    SecondaryData<typename ShapeMatrices::ShapeType> _secondary_data;
};

}  // namespace SmallDeformation
}  // namespace ProcessLib
