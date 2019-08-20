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

#include "MaterialLib/SolidModels/LinearElasticIsotropic.h"
#include "MaterialLib/SolidModels/LinearElasticIsotropicPhaseField.h"
#include "MaterialLib/SolidModels/SelectSolidConstitutiveRelation.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "ParameterLib/SpatialPosition.h"
#include "ProcessLib/Deformation/BMatrixPolicy.h"
#include "ProcessLib/Deformation/LinearBMatrix.h"
#include "ProcessLib/Utils/InitShapeMatrices.h"

#include "LocalAssemblerInterface.h"
#include "PhaseFieldProcessData.h"

namespace ProcessLib
{
namespace PhaseField
{
template <typename BMatricesType, typename ShapeMatrixType, int DisplacementDim>
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

    typename ShapeMatrixType::NodalRowVectorType N;
    typename ShapeMatrixType::GlobalDimNodalMatrixType dNdx;

    typename BMatricesType::KelvinVectorType eps, eps_prev;

    typename BMatricesType::KelvinVectorType sigma_tensile, sigma_compressive,
        sigma;
    double strain_energy_tensile, elastic_energy;

    MaterialLib::Solids::MechanicsBase<DisplacementDim> const& solid_material;
    std::unique_ptr<typename MaterialLib::Solids::MechanicsBase<
        DisplacementDim>::MaterialStateVariables>
        material_state_variables;

    typename BMatricesType::KelvinMatrixType C_tensile, C_compressive;
    double integration_weight;
    double history_variable, history_variable_prev;

    void pushBackState()
    {
        history_variable_prev =
            std::max(history_variable_prev, history_variable);
        eps_prev = eps;
        material_state_variables->pushBackState();
    }

    template <typename DisplacementVectorType>
    void updateConstitutiveRelation(double const t,
                                    ParameterLib::SpatialPosition const& x,
                                    double const /*dt*/,
                                    DisplacementVectorType const& /*u*/,
                                    double const degradation)
    {
        auto linear_elastic_mp =
            static_cast<MaterialLib::Solids::LinearElasticIsotropic<
                DisplacementDim> const&>(solid_material)
                .getMaterialProperties();

        auto const bulk_modulus = linear_elastic_mp.bulk_modulus(t, x);
        auto const mu = linear_elastic_mp.mu(t, x);

        std::tie(sigma, sigma_tensile, C_tensile, C_compressive,
                 strain_energy_tensile, elastic_energy) =
            MaterialLib::Solids::Phasefield::calculateDegradedStress<
                DisplacementDim>(degradation, bulk_modulus, mu, eps);
    }
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

    static constexpr int kelvin_vector_size =
        MathLib::KelvinVector::size<DisplacementDim>();
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

    PhaseFieldLocalAssembler(PhaseFieldLocalAssembler const&) = delete;
    PhaseFieldLocalAssembler(PhaseFieldLocalAssembler&&) = delete;

    PhaseFieldLocalAssembler(
        MeshLib::Element const& e,
        std::size_t const /*local_matrix_size*/,
        bool const is_axially_symmetric,
        unsigned const integration_order,
        PhaseFieldProcessData<DisplacementDim>& process_data)
        : _process_data(process_data),
          _integration_method(integration_order),
          _element(e),
          _is_axially_symmetric(is_axially_symmetric)
    {
        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        _ip_data.reserve(n_integration_points);
        _secondary_data.N.resize(n_integration_points);

        auto& solid_material =
            MaterialLib::Solids::selectSolidConstitutiveRelation(
                _process_data.solid_materials,
                _process_data.material_ids,
                e.getID());
        if (!dynamic_cast<MaterialLib::Solids::LinearElasticIsotropic<
                DisplacementDim> const*>(&solid_material))
        {
            OGS_FATAL(
                "For phasefield process only linear elastic solid material "
                "support is implemented.");
        }

        auto const shape_matrices =
            initShapeMatrices<ShapeFunction, ShapeMatricesType,
                              IntegrationMethod, DisplacementDim>(
                e, is_axially_symmetric, _integration_method);

        ParameterLib::SpatialPosition x_position;
        x_position.setElementID(_element.getID());

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            _ip_data.emplace_back(solid_material);
            auto& ip_data = _ip_data[ip];
            ip_data.integration_weight =
                _integration_method.getWeightedPoint(ip).getWeight() *
                shape_matrices[ip].integralMeasure * shape_matrices[ip].detJ;

            static const int kelvin_vector_size =
                MathLib::KelvinVector::size<DisplacementDim>();
            ip_data.eps.setZero(kelvin_vector_size);
            ip_data.eps_prev.resize(kelvin_vector_size);
            ip_data.C_tensile.setZero(kelvin_vector_size, kelvin_vector_size);
            ip_data.C_compressive.setZero(kelvin_vector_size,
                                          kelvin_vector_size);
            ip_data.sigma_tensile.setZero(kelvin_vector_size);
            ip_data.sigma_compressive.setZero(kelvin_vector_size);
            ip_data.history_variable =
                _process_data.history_field(0, x_position)[0];
            ip_data.history_variable_prev =
                _process_data.history_field(0, x_position)[0];
            ip_data.sigma.setZero(kelvin_vector_size);
            ip_data.strain_energy_tensile = 0.0;
            ip_data.elastic_energy = 0.0;

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

    void assembleWithJacobianForStaggeredScheme(
        double const t, std::vector<double> const& local_xdot,
        const double dxdot_dx, const double dx_dx,
        std::vector<double>& local_M_data, std::vector<double>& local_K_data,
        std::vector<double>& local_b_data, std::vector<double>& local_Jac_data,
        LocalCoupledSolutions const& local_coupled_solutions) override;

    void initializeConcrete() override
    {
        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            _ip_data[ip].pushBackState();
        }
    }

    void postTimestepConcrete(std::vector<double> const& /*local_x*/
                              ) override
    {
        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            _ip_data[ip].pushBackState();
        }
    }

    void computeCrackIntegral(
        std::size_t mesh_item_id,
        std::vector<
            std::reference_wrapper<NumLib::LocalToGlobalIndexMap>> const&
            dof_tables,
        GlobalVector const& x, double const t, double& crack_volume,
        CoupledSolutionsForStaggeredScheme const* const cpl_xs) override;

    void computeEnergy(
        std::size_t mesh_item_id,
        std::vector<
            std::reference_wrapper<NumLib::LocalToGlobalIndexMap>> const&
            dof_tables,
        GlobalVector const& x, double const t, double& elastic_energy,
        double& surface_energy, double& pressure_work,
        CoupledSolutionsForStaggeredScheme const* const cpl_xs) override;

    Eigen::Map<const Eigen::RowVectorXd> getShapeMatrix(
        const unsigned integration_point) const override
    {
        auto const& N = _secondary_data.N[integration_point];

        // assumes N is stored contiguously in memory
        return Eigen::Map<const Eigen::RowVectorXd>(N.data(), N.size());
    }

private:
    std::vector<double> const& getIntPtSigma(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        static const int kelvin_vector_size =
            MathLib::KelvinVector::size<DisplacementDim>();
        auto const num_intpts = _ip_data.size();

        cache.clear();
        auto cache_mat = MathLib::createZeroedMatrix<Eigen::Matrix<
            double, kelvin_vector_size, Eigen::Dynamic, Eigen::RowMajor>>(
            cache, kelvin_vector_size, num_intpts);

        for (unsigned ip = 0; ip < num_intpts; ++ip)
        {
            auto const& sigma = _ip_data[ip].sigma;
            cache_mat.col(ip) =
                MathLib::KelvinVector::kelvinVectorToSymmetricTensor(sigma);
        }

        return cache;
    }

    std::vector<double> const& getIntPtEpsilon(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        auto const kelvin_vector_size =
            MathLib::KelvinVector::size<DisplacementDim>();
        auto const num_intpts = _ip_data.size();

        cache.clear();
        auto cache_mat = MathLib::createZeroedMatrix<Eigen::Matrix<
            double, kelvin_vector_size, Eigen::Dynamic, Eigen::RowMajor>>(
            cache, kelvin_vector_size, num_intpts);

        for (unsigned ip = 0; ip < num_intpts; ++ip)
        {
            auto const& eps = _ip_data[ip].eps;
            cache_mat.col(ip) =
                MathLib::KelvinVector::kelvinVectorToSymmetricTensor(eps);
        }

        return cache;
    }

    void assembleWithJacobianPhaseFiledEquations(
        double const t, std::vector<double> const& local_xdot,
        const double dxdot_dx, const double dx_dx,
        std::vector<double>& local_M_data, std::vector<double>& local_K_data,
        std::vector<double>& local_b_data, std::vector<double>& local_Jac_data,
        LocalCoupledSolutions const& local_coupled_solutions);

    void assembleWithJacobianForDeformationEquations(
        double const t, std::vector<double> const& local_xdot,
        const double dxdot_dx, const double dx_dx,
        std::vector<double>& local_M_data, std::vector<double>& local_K_data,
        std::vector<double>& local_b_data, std::vector<double>& local_Jac_data,
        LocalCoupledSolutions const& local_coupled_solutions);

    PhaseFieldProcessData<DisplacementDim>& _process_data;

    std::vector<
        IntegrationPointData<BMatricesType, ShapeMatricesType, DisplacementDim>,
        Eigen::aligned_allocator<IntegrationPointData<
            BMatricesType, ShapeMatricesType, DisplacementDim>>>
        _ip_data;

    IntegrationMethod _integration_method;
    MeshLib::Element const& _element;
    SecondaryData<typename ShapeMatrices::ShapeType> _secondary_data;
    bool const _is_axially_symmetric;

    static const int phasefield_index = 0;
    static const int phasefield_size = ShapeFunction::NPOINTS;
    static const int displacement_index = ShapeFunction::NPOINTS;
    static const int displacement_size =
        ShapeFunction::NPOINTS * DisplacementDim;
};

}  // namespace PhaseField
}  // namespace ProcessLib

#include "PhaseFieldFEM-impl.h"
