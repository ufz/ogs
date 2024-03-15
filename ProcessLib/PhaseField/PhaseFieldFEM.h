/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>
#include <vector>

#include "LocalAssemblerInterface.h"
#include "MaterialLib/SolidModels/LinearElasticIsotropic.h"
#include "MaterialLib/SolidModels/LinearElasticIsotropicPhaseField.h"
#include "MaterialLib/SolidModels/LinearElasticOrthotropic.h"
#include "MaterialLib/SolidModels/LinearElasticOrthotropicPhaseField.h"
#include "MaterialLib/SolidModels/SelectSolidConstitutiveRelation.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/InitShapeMatrices.h"
#include "NumLib/Fem/Integration/GenericIntegrationMethod.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "ParameterLib/CoordinateSystem.h"
#include "ParameterLib/SpatialPosition.h"
#include "PhaseFieldProcessData.h"
#include "ProcessLib/Deformation/BMatrixPolicy.h"
#include "ProcessLib/Deformation/LinearBMatrix.h"
#include "ProcessLib/Utils/SetOrGetIntegrationPointData.h"

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

    typename BMatricesType::KelvinVectorType eps, eps_prev, eps_tensile;

    typename BMatricesType::KelvinVectorType sigma_tensile, sigma_compressive,
        sigma;
    double strain_energy_tensile, elastic_energy;

    MaterialLib::Solids::MechanicsBase<DisplacementDim> const& solid_material;
    std::unique_ptr<typename MaterialLib::Solids::MechanicsBase<
        DisplacementDim>::MaterialStateVariables>
        material_state_variables;

    typename BMatricesType::KelvinMatrixType D, C_tensile, C_compressive;

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
                                    double const degradation,
                                    EnergySplitModel const energy_split_model)
    {
        auto linear_elastic_mp =
            static_cast<MaterialLib::Solids::LinearElasticIsotropic<
                DisplacementDim> const&>(solid_material)
                .getMaterialProperties();

        auto const bulk_modulus = linear_elastic_mp.bulk_modulus(t, x);
        auto const mu = linear_elastic_mp.mu(t, x);

        switch (energy_split_model)
        {
            case EnergySplitModel::Isotropic:
            {
                std::tie(sigma, sigma_tensile, D, strain_energy_tensile,
                         elastic_energy, C_tensile, C_compressive) =
                    MaterialLib::Solids::Phasefield::
                        calculateIsotropicDegradedStress<DisplacementDim>(
                            degradation, bulk_modulus, mu, eps);
                break;
            }
            case EnergySplitModel::VolDev:
            {
                std::tie(sigma, sigma_tensile, D, strain_energy_tensile,
                         elastic_energy, C_tensile, C_compressive) =
                    MaterialLib::Solids::Phasefield::
                        calculateVolDevDegradedStress<DisplacementDim>(
                            degradation, bulk_modulus, mu, eps);
                break;
            }
            case EnergySplitModel::EffectiveStress:
            {
                std::tie(sigma, sigma_tensile, D, strain_energy_tensile,
                         elastic_energy, C_tensile,
                         C_compressive) = MaterialLib::Solids::Phasefield::
                    calculateIsotropicDegradedStressWithRankineEnergy<
                        DisplacementDim>(degradation, bulk_modulus, mu, eps);
                break;
            }
            case EnergySplitModel::OrthoVolDev:
            {
                double temperature = 0.;
                MathLib::KelvinVector::KelvinMatrixType<DisplacementDim>
                    C_ortho = static_cast<
                                  MaterialLib::Solids::LinearElasticOrthotropic<
                                      DisplacementDim> const&>(solid_material)
                                  .getElasticTensor(t, x, temperature);

                std::tie(eps_tensile, sigma, sigma_tensile, sigma_compressive,
                         D, strain_energy_tensile, elastic_energy, C_tensile,
                         C_compressive) = MaterialLib::Solids::Phasefield::
                    calculateOrthoVolDevDegradedStress<DisplacementDim>(
                        degradation, eps, C_ortho);
                break;
            }
            case EnergySplitModel::OrthoMasonry:
            {
                double temperature = 0.;
                MathLib::KelvinVector::KelvinMatrixType<DisplacementDim>
                    C_ortho = static_cast<
                                  MaterialLib::Solids::LinearElasticOrthotropic<
                                      DisplacementDim> const&>(solid_material)
                                  .getElasticTensor(t, x, temperature);

                std::tie(eps_tensile, sigma, sigma_tensile, D,
                         strain_energy_tensile, elastic_energy, C_tensile,
                         C_compressive) = MaterialLib::Solids::Phasefield::
                    calculateOrthoMasonryDegradedStress<DisplacementDim>(
                        degradation, eps, C_ortho);
                break;
            }
        }

        history_variable =
            std::max(history_variable_prev, strain_energy_tensile);
    }
};

/// Used for the extrapolation of the integration point values. It is
/// ordered (and stored) by integration points.
template <typename ShapeMatrixType>
struct SecondaryData
{
    std::vector<ShapeMatrixType, Eigen::aligned_allocator<ShapeMatrixType>> N;
};

template <typename ShapeFunction, int DisplacementDim>
class PhaseFieldLocalAssembler : public PhaseFieldLocalAssemblerInterface
{
private:
    static constexpr int displacement_index = 0;
    static constexpr int displacement_size =
        ShapeFunction::NPOINTS * DisplacementDim;
    static constexpr int phasefield_index =
        displacement_index + displacement_size;
    static constexpr int phasefield_size = ShapeFunction::NPOINTS;

public:
    using ShapeMatricesType =
        ShapeMatrixPolicyType<ShapeFunction, DisplacementDim>;
    // Types for displacement.
    // (Higher order elements = ShapeFunction).
    using ShapeMatrices = typename ShapeMatricesType::ShapeMatrices;
    using BMatricesType = BMatrixPolicyType<ShapeFunction, DisplacementDim>;

    using NodalForceVectorType = typename BMatricesType::NodalForceVectorType;

    using DeformationVector =
        typename ShapeMatricesType::template VectorType<displacement_size>;
    using PhaseFieldVector =
        typename ShapeMatricesType::template VectorType<phasefield_size>;
    using PhaseFieldMatrix =
        typename ShapeMatricesType::template MatrixType<phasefield_size,
                                                        phasefield_size>;
    using NodalMatrixType = typename ShapeMatricesType::NodalMatrixType;
    using NodalVectorType = typename ShapeMatricesType::NodalVectorType;
    using IpData =
        IntegrationPointData<BMatricesType, ShapeMatricesType, DisplacementDim>;

    PhaseFieldLocalAssembler(PhaseFieldLocalAssembler const&) = delete;
    PhaseFieldLocalAssembler(PhaseFieldLocalAssembler&&) = delete;

    PhaseFieldLocalAssembler(
        MeshLib::Element const& e,
        std::size_t const /*local_matrix_size*/,
        NumLib::GenericIntegrationMethod const& integration_method,
        bool const is_axially_symmetric,
        PhaseFieldProcessData<DisplacementDim>& process_data)
        : _process_data(process_data),
          _integration_method(integration_method),
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

        auto const shape_matrices =
            NumLib::initShapeMatrices<ShapeFunction, ShapeMatricesType,
                                      DisplacementDim>(e, is_axially_symmetric,
                                                       _integration_method);

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
                MathLib::KelvinVector::kelvin_vector_dimensions(
                    DisplacementDim);
            ip_data.eps_tensile.setZero(kelvin_vector_size);
            ip_data.eps.setZero(kelvin_vector_size);
            ip_data.eps_prev.resize(kelvin_vector_size);
            ip_data.D.setZero(kelvin_vector_size, kelvin_vector_size);

            ip_data.sigma_tensile.setZero(kelvin_vector_size);
            ip_data.sigma_compressive.setZero(kelvin_vector_size);
            ip_data.sigma.setZero(kelvin_vector_size);
            ip_data.strain_energy_tensile = 0.0;
            ip_data.elastic_energy = 0.0;

            ip_data.N = shape_matrices[ip].N;
            ip_data.dNdx = shape_matrices[ip].dNdx;

            _secondary_data.N[ip] = shape_matrices[ip].N;
        }
    }

    void assemble(double const /*t*/, double const /*dt*/,
                  std::vector<double> const& /*local_x*/,
                  std::vector<double> const& /*local_x_prev*/,
                  std::vector<double>& /*local_M_data*/,
                  std::vector<double>& /*local_K_data*/,
                  std::vector<double>& /*local_rhs_data*/) override
    {
        OGS_FATAL(
            "PhaseFieldLocalAssembler: assembly without jacobian is not "
            "implemented.");
    }

    void assembleWithJacobianForStaggeredScheme(
        double const t, double const dt, Eigen::VectorXd const& local_x,
        Eigen::VectorXd const& local_x_prev, int const process_id,
        std::vector<double>& local_M_data, std::vector<double>& local_K_data,
        std::vector<double>& local_b_data,
        std::vector<double>& local_Jac_data) override;

    void initializeConcrete() override
    {
        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            _ip_data[ip].pushBackState();
        }
    }

    void postTimestepConcrete(Eigen::VectorXd const& /*local_x*/,
                              Eigen::VectorXd const& /*local_x_prev*/,
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

    void computeCrackIntegral(
        std::size_t mesh_item_id,
        std::vector<
            std::reference_wrapper<NumLib::LocalToGlobalIndexMap>> const&
            dof_tables,
        std::vector<GlobalVector*> const& x, double const t,
        double& crack_volume) override;

    void computeEnergy(std::size_t mesh_item_id,
                       std::vector<std::reference_wrapper<
                           NumLib::LocalToGlobalIndexMap>> const& dof_tables,
                       std::vector<GlobalVector*> const& x, double const t,
                       double& elastic_energy, double& surface_energy,
                       double& pressure_work) override;

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
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        return ProcessLib::getIntegrationPointKelvinVectorData<DisplacementDim>(
            _ip_data, &IpData::sigma, cache);
    }

    std::vector<double> const& getIntPtSigmaTensile(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        return ProcessLib::getIntegrationPointKelvinVectorData<DisplacementDim>(
            _ip_data, &IpData::sigma_tensile, cache);
    }

    std::vector<double> const& getIntPtSigmaCompressive(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        return ProcessLib::getIntegrationPointKelvinVectorData<DisplacementDim>(
            _ip_data, &IpData::sigma_compressive, cache);
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

    std::vector<double> const& getIntPtEpsilonTensile(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        return ProcessLib::getIntegrationPointKelvinVectorData<DisplacementDim>(
            _ip_data, &IpData::eps_tensile, cache);
    }

    void assembleWithJacobianPhaseFieldEquations(
        double const t, double const dt, Eigen::VectorXd const& local_x,
        std::vector<double>& local_b_data, std::vector<double>& local_Jac_data);

    void assembleWithJacobianForDeformationEquations(
        double const t, double const dt, Eigen::VectorXd const& local_x,
        std::vector<double>& local_b_data, std::vector<double>& local_Jac_data);

    PhaseFieldProcessData<DisplacementDim>& _process_data;

    std::vector<IpData, Eigen::aligned_allocator<IpData>> _ip_data;

    NumLib::GenericIntegrationMethod const& _integration_method;
    MeshLib::Element const& _element;
    SecondaryData<typename ShapeMatrices::ShapeType> _secondary_data;
    bool const _is_axially_symmetric;

    static const int phase_process_id = 1;
    static const int mechanics_process_id = 0;
};

}  // namespace PhaseField
}  // namespace ProcessLib

#include "PhaseFieldFEM-impl.h"
