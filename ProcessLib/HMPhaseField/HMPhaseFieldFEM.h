// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <memory>
#include <vector>

#include "HMPhaseFieldProcessData.h"
#include "LocalAssemblerInterface.h"
#include "MaterialLib/SolidModels/LinearElasticIsotropic.h"
#include "MaterialLib/SolidModels/LinearElasticIsotropicPhaseField.h"
#include "MaterialLib/SolidModels/LinearElasticOrthotropic.h"
#include "MaterialLib/SolidModels/LinearElasticOrthotropicPhaseField.h"
#include "MaterialLib/SolidModels/PhaseFieldBase.h"
#include "MaterialLib/SolidModels/SelectSolidConstitutiveRelation.h"
#include "MathLib/EigenBlockMatrixView.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "MeshLib/Elements/Element.h"
#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/InitShapeMatrices.h"
#include "NumLib/Fem/Integration/GenericIntegrationMethod.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "ParameterLib/CoordinateSystem.h"
#include "ParameterLib/SpatialPosition.h"
#include "ProcessLib/Deformation/BMatrixPolicy.h"
#include "ProcessLib/Deformation/LinearBMatrix.h"
#include "ProcessLib/Utils/SetOrGetIntegrationPointData.h"

namespace ProcessLib
{
namespace HMPhaseField
{
namespace MPL = MaterialPropertyLib;
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
    double width_ip;
    Eigen::Vector<double, DisplacementDim> normal_ip;

    MaterialLib::Solids::MechanicsBase<DisplacementDim> const& solid_material;
    std::unique_ptr<typename MaterialLib::Solids::MechanicsBase<
        DisplacementDim>::MaterialStateVariables>
        material_state_variables;

    typename BMatricesType::KelvinMatrixType D, C_tensile, C_compressive;

    double integration_weight;
    double coupling_pressure = 0.0;
    double fracture_enhanced_porosity = 0.0;
    double biot_coefficient, biot_coefficient_prev, biot_modulus_inv,
        biot_modulus_inv_prev;

    void pushBackState()
    {
        eps_prev = eps;
        biot_coefficient_prev = biot_coefficient;
        biot_modulus_inv_prev = biot_modulus_inv;
        material_state_variables->pushBackState();
    }

    template <typename DisplacementVectorType>
    void updateConstitutiveRelation(
        double const t,
        ParameterLib::SpatialPosition const& x,
        double const /*dt*/,
        DisplacementVectorType const& /*u*/,
        double const degradation,
        MaterialLib::Solids::Phasefield::EnergySplitModel const
            energy_split_model)
    {
        MaterialLib::Solids::Phasefield::calculateStress<
            decltype(sigma), decltype(D), DisplacementDim>(
            sigma, sigma_tensile, sigma_compressive, eps_tensile, D, C_tensile,
            C_compressive, strain_energy_tensile, elastic_energy, degradation,
            eps, energy_split_model, t, x, solid_material);
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
class HMPhaseFieldLocalAssembler : public HMPhaseFieldLocalAssemblerInterface
{
private:
    static constexpr int phasefield_index = 0;
    static constexpr int phasefield_size = ShapeFunction::NPOINTS;
    static constexpr int pressure_index = phasefield_index + phasefield_size;
    static constexpr int pressure_size = ShapeFunction::NPOINTS;
    static constexpr int displacement_index = pressure_index + pressure_size;
    static constexpr int displacement_size =
        ShapeFunction::NPOINTS * DisplacementDim;

public:
    using ShapeMatricesType =
        ShapeMatrixPolicyType<ShapeFunction, DisplacementDim>;

    using ShapeMatrices = typename ShapeMatricesType::ShapeMatrices;

    using BMatricesType = BMatrixPolicyType<ShapeFunction, DisplacementDim>;

    using NodalForceVectorType = typename BMatricesType::NodalForceVectorType;

    using GlobalDimVectorType = typename ShapeMatricesType::GlobalDimVectorType;

    using DeformationVector =
        typename ShapeMatricesType::template VectorType<displacement_size>;
    using DeformationMatrix =
        typename ShapeMatricesType::template MatrixType<displacement_size,
                                                        displacement_size>;
    using PhaseFieldVector =
        typename ShapeMatricesType::template VectorType<phasefield_size>;
    using PhaseFieldMatrix =
        typename ShapeMatricesType::template MatrixType<phasefield_size,
                                                        phasefield_size>;
    using PressureVector =
        typename ShapeMatricesType::template VectorType<pressure_size>;
    using PressureMatrix =
        typename ShapeMatricesType::template MatrixType<pressure_size,
                                                        pressure_size>;

    using NodalMatrixType = typename ShapeMatricesType::NodalMatrixType;
    using NodalVectorType = typename ShapeMatricesType::NodalVectorType;
    using IpData =
        IntegrationPointData<BMatricesType, ShapeMatricesType, DisplacementDim>;

    static int const KelvinVectorSize =
        MathLib::KelvinVector::kelvin_vector_dimensions(DisplacementDim);

    static constexpr auto& N_u_op = MathLib::eigenBlockMatrixView<
        DisplacementDim, typename ShapeMatricesType::NodalRowVectorType>;

    using Invariants = MathLib::KelvinVector::Invariants<KelvinVectorSize>;

    HMPhaseFieldLocalAssembler(HMPhaseFieldLocalAssembler const&) = delete;
    HMPhaseFieldLocalAssembler(HMPhaseFieldLocalAssembler&&) = delete;

    HMPhaseFieldLocalAssembler(
        MeshLib::Element const& e,
        std::size_t const /*local_matrix_size*/,
        NumLib::GenericIntegrationMethod const& integration_method,
        bool const is_axially_symmetric,
        HMPhaseFieldProcessData<DisplacementDim>& process_data)
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
            ip_data.C_tensile.setZero(kelvin_vector_size, kelvin_vector_size);
            ip_data.C_compressive.setZero(kelvin_vector_size,
                                          kelvin_vector_size);

            ip_data.sigma_tensile.setZero(kelvin_vector_size);
            ip_data.sigma_compressive.setZero(kelvin_vector_size);
            ip_data.sigma.setZero(kelvin_vector_size);
            ip_data.strain_energy_tensile = 0.0;
            ip_data.elastic_energy = 0.0;
            ip_data.width_ip = 0.0;

            ip_data.N = shape_matrices[ip].N;
            ip_data.dNdx = shape_matrices[ip].dNdx;

            _secondary_data.N[ip] = shape_matrices[ip].N;
        }
    }

    void assembleWithJacobianForStaggeredScheme(
        double const t, double const dt, Eigen::VectorXd const& local_x,
        Eigen::VectorXd const& local_x_prev, int const process_id,
        std::vector<double>& local_b_data,
        std::vector<double>& local_Jac_data) override;

    void initializeConcrete() override
    {
        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            _ip_data[ip].pushBackState();

            // Specify the direction of preexisting fracture if it exists. The
            // default value is zero.
            auto& normal_ip = _ip_data[ip].normal_ip;
            auto const fracture_normal =
                _process_data.specific_fracture_direction;
            normal_ip = fracture_normal;
        }
        // Specify the aperture of preexisting fractures if they exist. The
        // default value is zero.
        ParameterLib::SpatialPosition x_position;
        x_position.setElementID(_element.getID());
        auto const width_init = _process_data.width_init(0, x_position)[0];
        (*_process_data.width)[_element.getID()] = width_init;
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

    void postNonLinearSolverConcrete(Eigen::VectorXd const& local_x,
                                     Eigen::VectorXd const& local_x_prev,
                                     double const t, double const dt,
                                     int const process_id) override;

    void approximateFractureWidth(
        std::size_t mesh_item_id,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_tables,
        std::vector<GlobalVector*> const& x, double const t,
        double const dt) override;

    void computeEnergy(
        std::size_t mesh_item_id,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_tables,
        std::vector<GlobalVector*> const& x, double const t,
        double& elastic_energy, double& surface_energy,
        double& pressure_work) override;

    inline double heaviside(double const v) { return (v < 0) ? 0.0 : 1.0; }

    Eigen::Map<const Eigen::RowVectorXd> getShapeMatrix(
        const unsigned integration_point) const override
    {
        auto const& N = _secondary_data.N[integration_point];

        // assumes N is stored contiguously in memory
        return Eigen::Map<const Eigen::RowVectorXd>(N.data(), N.size());
    }

    std::vector<double> const& getIntPtWidth(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const override;

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

    std::vector<double> const& getIntPtEpsilon(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        return ProcessLib::getIntegrationPointKelvinVectorData<DisplacementDim>(
            _ip_data, &IpData::eps, cache);
    }

    void assembleWithJacobianHydroEquations(
        const double t, double const dt, Eigen::VectorXd const& local_x,
        Eigen::VectorXd const& local_x_prev, std::vector<double>& local_b_data,
        std::vector<double>& local_Jac_data);

    void assembleWithJacobianPhaseFieldEquations(
        double const t, double const dt, Eigen::VectorXd const& local_x,
        std::vector<double>& local_b_data, std::vector<double>& local_Jac_data);

    void assembleWithJacobianForDeformationEquations(
        double const t, double const dt, Eigen::VectorXd const& local_x,
        std::vector<double>& local_b_data, std::vector<double>& local_Jac_data);

    HMPhaseFieldProcessData<DisplacementDim>& _process_data;

    std::vector<IpData, Eigen::aligned_allocator<IpData>> _ip_data;

    NumLib::GenericIntegrationMethod const& _integration_method;
    MeshLib::Element const& _element;
    SecondaryData<typename ShapeMatrices::ShapeType> _secondary_data;
    bool const _is_axially_symmetric;
};

}  // namespace HMPhaseField
}  // namespace ProcessLib

#include "HMPhaseFieldFEM-impl.h"
