/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>
#include <vector>

#include "MaterialLib/PhysicalConstant.h"
#include "MaterialLib/SolidModels/LinearElasticIsotropic.h"
#include "MathLib/KelvinVector.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "NumLib/DOF/DOFTableUtil.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "ProcessLib/Deformation/BMatrixPolicy.h"
#include "ProcessLib/Deformation/LinearBMatrix.h"
#include "ProcessLib/LocalAssemblerTraits.h"
#include "ProcessLib/Parameter/Parameter.h"
#include "ProcessLib/Utils/InitShapeMatrices.h"

#include "ThermoHydroMechanicsProcessData.h"
#include "LocalAssemblerInterface.h"

namespace ProcessLib
{
namespace ThermoHydroMechanics
{
template <typename BMatricesType, typename ShapeMatrixTypeDisplacement,
          typename ShapeMatricesTypePressure, int DisplacementDim, int NPoints>
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

    typename ShapeMatrixTypeDisplacement::template MatrixType<
        DisplacementDim, NPoints * DisplacementDim>
        N_u_op;
    typename BMatricesType::KelvinVectorType sigma_eff, sigma_eff_prev;
    typename BMatricesType::KelvinVectorType eps, eps_prev;
    typename BMatricesType::KelvinVectorType eps_m, eps_m_prev;

    typename ShapeMatrixTypeDisplacement::NodalRowVectorType N_u;
    typename ShapeMatrixTypeDisplacement::GlobalDimNodalMatrixType dNdx_u;

    typename ShapeMatricesTypePressure::NodalRowVectorType N_p;
    typename ShapeMatricesTypePressure::GlobalDimNodalMatrixType dNdx_p;

    MaterialLib::Solids::MechanicsBase<DisplacementDim> const& solid_material;
    std::unique_ptr<typename MaterialLib::Solids::MechanicsBase<
        DisplacementDim>::MaterialStateVariables>
        material_state_variables;
    double integration_weight;

    void pushBackState()
    {
        eps_m_prev = eps_m;
        sigma_eff_prev = sigma_eff;
        material_state_variables->pushBackState();
    }

    template <typename DisplacementVectorType>
    typename BMatricesType::KelvinMatrixType updateConstitutiveRelation(
        double const t,
        SpatialPosition const& x_position,
        double const dt,
        DisplacementVectorType const& /*u*/,
        double const T)
    {
        auto&& solution = solid_material.integrateStress(
            t, x_position, dt, eps_m_prev, eps_m, sigma_eff_prev,
            *material_state_variables, T);

        if (!solution)
            OGS_FATAL("Computation of local constitutive relation failed.");

        MathLib::KelvinVector::KelvinMatrixType<DisplacementDim> C;
        std::tie(sigma_eff, material_state_variables, C) = std::move(*solution);

        return C;
    }

    template <typename DisplacementVectorType>
    typename BMatricesType::KelvinMatrixType updateConstitutiveRelationThermal(
        double const t,
        SpatialPosition const& x_position,
        double const dt,
        DisplacementVectorType const& /*u*/,
        double const T,
        double const thermal_strain)
    {
        auto const& identity2 = MathLib::KelvinVector::Invariants<
            MathLib::KelvinVector::KelvinVectorDimensions<
                DisplacementDim>::value>::identity2;

        // assume isotropic thermal expansion
        eps_m.noalias() = eps - thermal_strain * identity2;
        auto&& solution = solid_material.integrateStress(
            t, x_position, dt, eps_m_prev, eps_m, sigma_eff_prev,
            *material_state_variables, T);

        if (!solution)
            OGS_FATAL("Computation of local constitutive relation failed.");

        MathLib::KelvinVector::KelvinMatrixType<DisplacementDim> C;
        std::tie(sigma_eff, material_state_variables, C) = std::move(*solution);

        return C;
    }

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

/// Used by for extrapolation of the integration point values. It is ordered
/// (and stored) by integration points.
template <typename ShapeMatrixType>
struct SecondaryData
{
    std::vector<ShapeMatrixType, Eigen::aligned_allocator<ShapeMatrixType>> N_u;
};

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
class ThermoHydroMechanicsLocalAssembler : public LocalAssemblerInterface
{
public:
    using ShapeMatricesTypeDisplacement =
        ShapeMatrixPolicyType<ShapeFunctionDisplacement, DisplacementDim>;

    // Types for pressure.
    using ShapeMatricesTypePressure =
        ShapeMatrixPolicyType<ShapeFunctionPressure, DisplacementDim>;

    static int const KelvinVectorSize =
        MathLib::KelvinVector::KelvinVectorDimensions<DisplacementDim>::value;
    using Invariants = MathLib::KelvinVector::Invariants<KelvinVectorSize>;

    ThermoHydroMechanicsLocalAssembler(
        ThermoHydroMechanicsLocalAssembler const&) = delete;
    ThermoHydroMechanicsLocalAssembler(ThermoHydroMechanicsLocalAssembler&&) =
        delete;

    ThermoHydroMechanicsLocalAssembler(
        MeshLib::Element const& e,
        std::size_t const /*local_matrix_size*/,
        bool const is_axially_symmetric,
        unsigned const integration_order,
        ThermoHydroMechanicsProcessData<DisplacementDim>& process_data);

    void assemble(double const /*t*/, std::vector<double> const& /*local_x*/,
                  std::vector<double>& /*local_M_data*/,
                  std::vector<double>& /*local_K_data*/,
                  std::vector<double>& /*local_rhs_data*/) override
    {
        OGS_FATAL(
            "ThermoHydroMechanicsLocalAssembler: assembly without Jacobian is "
            "not implemented.");
    }

    void assembleWithJacobian(double const t,
                              std::vector<double> const& local_x,
                              std::vector<double> const& local_xdot,
                              const double /*dxdot_dx*/, const double /*dx_dx*/,
                              std::vector<double>& /*local_M_data*/,
                              std::vector<double>& /*local_K_data*/,
                              std::vector<double>& local_rhs_data,
                              std::vector<double>& local_Jac_data) override;

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

    void computeSecondaryVariableConcrete(
        double const t, std::vector<double> const& local_x) override;
    void postNonLinearSolverConcrete(std::vector<double> const& local_x,
                                     double const t,
                                     bool const use_monolithic_scheme) override;

    Eigen::Map<const Eigen::RowVectorXd> getShapeMatrix(
        const unsigned integration_point) const override
    {
        auto const& N_u = _secondary_data.N_u[integration_point];

        // assumes N is stored contiguously in memory
        return Eigen::Map<const Eigen::RowVectorXd>(N_u.data(), N_u.size());
    }

    std::vector<double> const& getIntPtDarcyVelocity(
        const double t,
        GlobalVector const& current_solution,
        NumLib::LocalToGlobalIndexMap const& dof_table,
        std::vector<double>& cache) const override;

private:
    std::size_t setSigma(double const* values)
        {
            auto const kelvin_vector_size =
                MathLib::KelvinVector::KelvinVectorDimensions<
                    DisplacementDim>::value;
            unsigned const n_integration_points =
                _integration_method.getNumberOfPoints();

            std::vector<double> ip_sigma_values;
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
                auto const& sigma = _ip_data[ip].sigma_eff;
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
                    auto const& sigma = _ip_data[ip].sigma_eff;
                    cache_mat.col(ip) =
                        MathLib::KelvinVector::kelvinVectorToSymmetricTensor(sigma);
                }

                return cache;
            }

            virtual std::vector<double> const& getIntPtEpsilon(
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

private:
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
    bool const _is_axially_symmetric;
    SecondaryData<
        typename ShapeMatricesTypeDisplacement::ShapeMatrices::ShapeType>
        _secondary_data;

    static const int temperature_index = 0;
    static const int temperature_size = ShapeFunctionPressure::NPOINTS;
    static const int pressure_index = ShapeFunctionPressure::NPOINTS;
    static const int pressure_size = ShapeFunctionPressure::NPOINTS;
    static const int displacement_index = ShapeFunctionPressure::NPOINTS * 2;
    static const int displacement_size =
        ShapeFunctionDisplacement::NPOINTS * DisplacementDim;
};

}  // namespace ThermoHydroMechanics
}  // namespace ProcessLib

#include "ThermoHydroMechanicsFEM-impl.h"
