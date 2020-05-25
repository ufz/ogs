/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>
#include <vector>

#include "LocalAssemblerInterface.h"
#include "MaterialLib/PhysicalConstant.h"
#include "MaterialLib/SolidModels/SelectSolidConstitutiveRelation.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "NumLib/Extrapolation/ExtrapolatableElement.h"
#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "ParameterLib/Parameter.h"
#include "ProcessLib/Deformation/BMatrixPolicy.h"
#include "ProcessLib/Deformation/GMatrixPolicy.h"
#include "ProcessLib/Deformation/LinearBMatrix.h"
#include "ProcessLib/LocalAssemblerInterface.h"
#include "ProcessLib/LocalAssemblerTraits.h"
#include "ProcessLib/Utils/InitShapeMatrices.h"
#include "ProcessLib/Utils/SetOrGetIntegrationPointData.h"
#include "SmallDeformationProcessData.h"

namespace ProcessLib
{
namespace SmallDeformation
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

    typename BMatricesType::KelvinVectorType sigma, sigma_prev;
    typename BMatricesType::KelvinVectorType eps, eps_prev;
    double free_energy_density = 0;

    MaterialLib::Solids::MechanicsBase<DisplacementDim> const& solid_material;
    std::unique_ptr<typename MaterialLib::Solids::MechanicsBase<
        DisplacementDim>::MaterialStateVariables>
        material_state_variables;

    double integration_weight;
    typename ShapeMatricesType::NodalRowVectorType N;
    typename ShapeMatricesType::GlobalDimNodalMatrixType dNdx;

    void pushBackState()
    {
        eps_prev = eps;
        sigma_prev = sigma;
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

    SmallDeformationLocalAssembler(SmallDeformationLocalAssembler const&) =
        delete;
    SmallDeformationLocalAssembler(SmallDeformationLocalAssembler&&) = delete;

    SmallDeformationLocalAssembler(
        MeshLib::Element const& e,
        std::size_t const /*local_matrix_size*/,
        bool const is_axially_symmetric,
        unsigned const integration_order,
        SmallDeformationProcessData<DisplacementDim>& process_data)
        : process_data_(process_data),
          integration_method_(integration_order),
          element_(e),
          is_axially_symmetric_(is_axially_symmetric)
    {
        unsigned const n_integration_points =
            integration_method_.getNumberOfPoints();

        ip_data_.reserve(n_integration_points);
        secondary_data_.N.resize(n_integration_points);

        auto const shape_matrices =
            initShapeMatrices<ShapeFunction, ShapeMatricesType,
                              IntegrationMethod, DisplacementDim>(
                e, is_axially_symmetric, integration_method_);

        auto& solid_material =
            MaterialLib::Solids::selectSolidConstitutiveRelation(
                process_data_.solid_materials,
                process_data_.material_ids,
                e.getID());

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            ip_data_.emplace_back(solid_material);
            auto& ip_data = ip_data_[ip];
            auto const& sm = shape_matrices[ip];
            ip_data_[ip].integration_weight =
                integration_method_.getWeightedPoint(ip).getWeight() *
                sm.integralMeasure * sm.detJ;

            ip_data.N = sm.N;
            ip_data.dNdx = sm.dNdx;

            static const int kelvin_vector_size =
                MathLib::KelvinVector::KelvinVectorDimensions<
                    DisplacementDim>::value;
            // Initialize current time step values
            ip_data.sigma.setZero(kelvin_vector_size);
            ip_data.eps.setZero(kelvin_vector_size);

            // Previous time step values are not initialized and are set later.
            ip_data.sigma_prev.resize(kelvin_vector_size);
            ip_data.eps_prev.resize(kelvin_vector_size);

            secondary_data_.N[ip] = shape_matrices[ip].N;
        }
    }

    /// Returns number of read integration points.
    std::size_t setIPDataInitialConditions(std::string const& name,
                                           double const* values,
                                           int const integration_order) override
    {
        if (integration_order !=
            static_cast<int>(integration_method_.getIntegrationOrder()))
        {
            OGS_FATAL(
                "Setting integration point initial conditions; The integration "
                "order of the local assembler for element {:d} is different "
                "from the integration order in the initial condition.",
                element_.getID());
        }

        if (name == "sigma_ip")
        {
            if (process_data_.initial_stress != nullptr)
            {
                OGS_FATAL(
                    "Setting initial conditions for stress from integration "
                    "point data and from a parameter '{:s}' is not possible "
                    "simultaneously.",
                    process_data_.initial_stress->name);
            }
            return setSigma(values);
        }

        return 0;
    }

    void initializeConcrete() override
    {
        unsigned const n_integration_points =
            integration_method_.getNumberOfPoints();
        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            auto& ip_data = ip_data_[ip];

            /// Set initial stress from parameter.
            if (process_data_.initial_stress != nullptr)
            {
                ParameterLib::SpatialPosition const x_position{
                    boost::none, element_.getID(), ip,
                    MathLib::Point3d(interpolateCoordinates<ShapeFunction,
                                                            ShapeMatricesType>(
                        element_, ip_data.N))};

                ip_data.sigma =
                    MathLib::KelvinVector::symmetricTensorToKelvinVector<
                        DisplacementDim>((*process_data_.initial_stress)(
                        std::numeric_limits<
                            double>::quiet_NaN() /* time independent */,
                        x_position));
            }

            ip_data.pushBackState();
        }
    }

    void assemble(double const /*t*/, double const /*dt*/,
                  std::vector<double> const& /*local_x*/,
                  std::vector<double> const& /*local_xdot*/,
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
                              std::vector<double> const& /*local_xdot*/,
                              const double /*dxdot_dx*/, const double /*dx_dx*/,
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

        unsigned const n_integration_points =
            integration_method_.getNumberOfPoints();

        ParameterLib::SpatialPosition x_position;
        x_position.setElementID(element_.getID());

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            x_position.setIntegrationPoint(ip);
            auto const& w = ip_data_[ip].integration_weight;
            auto const& N = ip_data_[ip].N;
            auto const& dNdx = ip_data_[ip].dNdx;

            typename ShapeMatricesType::template MatrixType<DisplacementDim,
                                                            displacement_size>
                N_u_op = ShapeMatricesType::template MatrixType<
                    DisplacementDim,
                    displacement_size>::Zero(DisplacementDim,
                                             displacement_size);
            for (int i = 0; i < DisplacementDim; ++i)
            {
                N_u_op
                    .template block<1, displacement_size / DisplacementDim>(
                        i, i * displacement_size / DisplacementDim)
                    .noalias() = N;
            }

            auto const x_coord =
                interpolateXCoordinate<ShapeFunction, ShapeMatricesType>(
                    element_, N);
            auto const B = LinearBMatrix::computeBMatrix<
                DisplacementDim, ShapeFunction::NPOINTS,
                typename BMatricesType::BMatrixType>(dNdx, N, x_coord,
                                                     is_axially_symmetric_);

            auto const& eps_prev = ip_data_[ip].eps_prev;
            auto const& sigma_prev = ip_data_[ip].sigma_prev;

            auto& eps = ip_data_[ip].eps;
            auto& sigma = ip_data_[ip].sigma;
            auto& state = ip_data_[ip].material_state_variables;

            eps.noalias() =
                B *
                Eigen::Map<typename BMatricesType::NodalForceVectorType const>(
                    local_x.data(), ShapeFunction::NPOINTS * DisplacementDim);

            auto&& solution = ip_data_[ip].solid_material.integrateStress(
                t, x_position, dt, eps_prev, eps, sigma_prev, *state,
                process_data_.reference_temperature);

            if (!solution)
            {
                OGS_FATAL("Computation of local constitutive relation failed.");
            }

            MathLib::KelvinVector::KelvinMatrixType<DisplacementDim> C;
            std::tie(sigma, state, C) = std::move(*solution);

            auto const rho = process_data_.solid_density(t, x_position)[0];
            auto const& b = process_data_.specific_body_force;
            local_b.noalias() -=
                (B.transpose() * sigma - N_u_op.transpose() * rho * b) * w;
            local_Jac.noalias() += B.transpose() * C * B * w;
        }
    }

    void postTimestepConcrete(std::vector<double> const& /*local_x*/,
                              double const t, double const dt) override
    {
        unsigned const n_integration_points =
            integration_method_.getNumberOfPoints();

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            ip_data_[ip].pushBackState();
        }

        ParameterLib::SpatialPosition x_position;
        x_position.setElementID(element_.getID());

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            x_position.setIntegrationPoint(ip);

            auto& ip_data = ip_data_[ip];

            // Update free energy density needed for material forces.
            ip_data.free_energy_density =
                ip_data.solid_material.computeFreeEnergyDensity(
                    t, x_position, dt, ip_data.eps, ip_data.sigma,
                    *ip_data.material_state_variables);
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
            GradientMatrixType>(local_x, nodal_values, integration_method_,
                                ip_data_, element_, is_axially_symmetric_);
    }

    Eigen::Map<const Eigen::RowVectorXd> getShapeMatrix(
        const unsigned integration_point) const override
    {
        auto const& N = secondary_data_.N[integration_point];

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
        cache.reserve(ip_data_.size());

        for (auto const& ip_data : ip_data_)
        {
            cache.push_back(ip_data.free_energy_density);
        }

        return cache;
    }

    std::size_t setSigma(double const* values)
    {
        return ProcessLib::setIntegrationPointKelvinVectorData<DisplacementDim>(
            values, ip_data_, &IpData::sigma);
    }

    // TODO (naumov) This method is same as getIntPtSigma but for arguments and
    // the ordering of the cache_mat.
    // There should be only one.
    std::vector<double> getSigma() const override
    {
        return ProcessLib::getIntegrationPointKelvinVectorData<DisplacementDim>(
            ip_data_, &IpData::sigma);
    }

    std::vector<double> const& getIntPtSigma(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        return ProcessLib::getIntegrationPointKelvinVectorData<DisplacementDim>(
            ip_data_, &IpData::sigma, cache);
    }

    std::vector<double> const& getIntPtEpsilon(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        return ProcessLib::getIntegrationPointKelvinVectorData<DisplacementDim>(
            ip_data_, &IpData::eps, cache);
    }

    unsigned getNumberOfIntegrationPoints() const override
    {
        return integration_method_.getNumberOfPoints();
    }

    typename MaterialLib::Solids::MechanicsBase<
        DisplacementDim>::MaterialStateVariables const&
    getMaterialStateVariablesAt(unsigned integration_point) const override
    {
        return *ip_data_[integration_point].material_state_variables;
    }

    void computeSecondaryVariableConcrete(
        double const /*t*/, double const /*dt*/,
        std::vector<double> const& /*x*/,
        std::vector<double> const& /*x_dot*/) override
    {
        int const elem_id = element_.getID();
        ParameterLib::SpatialPosition x_position;
        x_position.setElementID(elem_id);
        unsigned const n_integration_points =
            integration_method_.getNumberOfPoints();

        auto sigma_sum = MathLib::KelvinVector::tensorToKelvin<DisplacementDim>(
            Eigen::Matrix<double, 3, 3>::Zero());

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            x_position.setIntegrationPoint(ip);
            auto const& sigma = ip_data_[ip].sigma;
            sigma_sum += sigma;
        }

        Eigen::Matrix<double, 3, 3, 0, 3, 3> const sigma_avg =
            MathLib::KelvinVector::kelvinVectorToTensor(sigma_sum) /
            n_integration_points;

        Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 3, 3>> e_s(
            sigma_avg);

        Eigen::Map<Eigen::Vector3d>(
            &(*process_data_.principal_stress_values)[elem_id * 3], 3) =
            e_s.eigenvalues();

        auto eigen_vectors = e_s.eigenvectors();

        for (auto i = 0; i < 3; i++)
        {
            Eigen::Map<Eigen::Vector3d>(
                &(*process_data_.principal_stress_vector[i])[elem_id * 3], 3) =
                eigen_vectors.col(i);
        }
    }

private:
    SmallDeformationProcessData<DisplacementDim>& process_data_;

    std::vector<IpData, Eigen::aligned_allocator<IpData>> ip_data_;

    IntegrationMethod integration_method_;
    MeshLib::Element const& element_;
    SecondaryData<typename ShapeMatrices::ShapeType> secondary_data_;
    bool const is_axially_symmetric_;

    static const int displacement_size =
        ShapeFunction::NPOINTS * DisplacementDim;
};

}  // namespace SmallDeformation
}  // namespace ProcessLib
