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

#include "IntegrationPointData.h"
#include "LocalAssemblerInterface.h"
#include "MaterialLib/PhysicalConstant.h"
#include "MaterialLib/SolidModels/LinearElasticIsotropic.h"
#include "MathLib/EigenBlockMatrixView.h"
#include "MathLib/KelvinVector.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "NumLib/DOF/DOFTableUtil.h"
#include "NumLib/DOF/LocalDOF.h"
#include "NumLib/Fem/InitShapeMatrices.h"
#include "NumLib/Fem/Integration/GenericIntegrationMethod.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "ParameterLib/Parameter.h"
#include "ProcessLib/Deformation/BMatrixPolicy.h"
#include "ProcessLib/Deformation/LinearBMatrix.h"
#include "ProcessLib/LocalAssemblerTraits.h"
#include "ProcessLib/Utils/SetOrGetIntegrationPointData.h"
#include "ProcessLib/Utils/TransposeInPlace.h"
#include "ThermoHydroMechanicsProcessData.h"

namespace ProcessLib
{
namespace ThermoHydroMechanics
{
namespace MPL = MaterialPropertyLib;

/// Used by for extrapolation of the integration point values. It is ordered
/// (and stored) by integration points.
template <typename ShapeMatrixType>
struct SecondaryData
{
    std::vector<ShapeMatrixType, Eigen::aligned_allocator<ShapeMatrixType>> N_u;
};

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
class ThermoHydroMechanicsLocalAssembler
    : public LocalAssemblerInterface<DisplacementDim>
{
public:
    using ShapeMatricesTypeDisplacement =
        ShapeMatrixPolicyType<ShapeFunctionDisplacement, DisplacementDim>;

    // Types for pressure.
    using ShapeMatricesTypePressure =
        ShapeMatrixPolicyType<ShapeFunctionPressure, DisplacementDim>;

    using GlobalDimMatrixType =
        typename ShapeMatricesTypePressure::GlobalDimMatrixType;

    using GlobalDimVectorType =
        typename ShapeMatricesTypePressure::GlobalDimVectorType;

    static int const KelvinVectorSize =
        MathLib::KelvinVector::kelvin_vector_dimensions(DisplacementDim);
    using Invariants = MathLib::KelvinVector::Invariants<KelvinVectorSize>;

    using SymmetricTensor = Eigen::Matrix<double, KelvinVectorSize, 1>;

    static constexpr auto& N_u_op = MathLib::eigenBlockMatrixView<
        DisplacementDim,
        typename ShapeMatricesTypeDisplacement::NodalRowVectorType>;

    ThermoHydroMechanicsLocalAssembler(
        ThermoHydroMechanicsLocalAssembler const&) = delete;
    ThermoHydroMechanicsLocalAssembler(ThermoHydroMechanicsLocalAssembler&&) =
        delete;

    ThermoHydroMechanicsLocalAssembler(
        MeshLib::Element const& e,
        std::size_t const /*local_matrix_size*/,
        NumLib::GenericIntegrationMethod const& integration_method,
        bool const is_axially_symmetric,
        ThermoHydroMechanicsProcessData<DisplacementDim>& process_data);

    /// Returns number of read integration points.
    std::size_t setIPDataInitialConditions(
        std::string_view const name,
        double const* values,
        int const integration_order) override;

    void assemble(double const /*t*/, double const /*dt*/,
                  std::vector<double> const& /*local_x*/,
                  std::vector<double> const& /*local_x_prev*/,
                  std::vector<double>& /*local_M_data*/,
                  std::vector<double>& /*local_K_data*/,
                  std::vector<double>& /*local_rhs_data*/) override
    {
        OGS_FATAL(
            "ThermoHydroMechanicsLocalAssembler: assembly without Jacobian is "
            "not implemented.");
    }

    void assembleWithJacobian(double const t, double const dt,
                              std::vector<double> const& local_x,
                              std::vector<double> const& local_x_prev,
                              std::vector<double>& /*local_M_data*/,
                              std::vector<double>& /*local_K_data*/,
                              std::vector<double>& local_rhs_data,
                              std::vector<double>& local_Jac_data) override;

    void initializeConcrete() override
    {
        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            auto& ip_data = _ip_data[ip];

            ParameterLib::SpatialPosition const x_position{
                std::nullopt, _element.getID(), ip,
                MathLib::Point3d(
                    NumLib::interpolateCoordinates<
                        ShapeFunctionDisplacement,
                        ShapeMatricesTypeDisplacement>(_element, ip_data.N_u))};

            /// Set initial stress from parameter.
            if (_process_data.initial_stress.value)
            {
                ip_data.sigma_eff =
                    MathLib::KelvinVector::symmetricTensorToKelvinVector<
                        DisplacementDim>((*_process_data.initial_stress.value)(
                        std::numeric_limits<double>::quiet_NaN() /* time
                                                                    independent
                                                                  */
                        ,
                        x_position));
            }

            double const t = 0;  // TODO (naumov) pass t from top
            ip_data.solid_material.initializeInternalStateVariables(
                t, x_position, *ip_data.material_state_variables);

            ip_data.pushBackState();
        }
    }

    void setInitialConditionsConcrete(Eigen::VectorXd const local_x,
                                      double const t,
                                      int const process_id) override;

    void preTimestepConcrete(std::vector<double> const& /*local_x*/,
                             double const /*t*/, double const /*dt*/) override
    {
        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            _ip_data_output[ip].velocity.setConstant(
                DisplacementDim, std::numeric_limits<double>::quiet_NaN());
        }
    }

    void postTimestepConcrete(Eigen::VectorXd const& local_x,
                              Eigen::VectorXd const& local_x_prev,
                              double const t, double const dt,
                              int const /*process_id*/) override
    {
        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        auto const [T_prev, p_prev, u_prev] = localDOF(local_x_prev);

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            auto& ip_data = _ip_data[ip];
            auto const& N_u = ip_data.N_u;
            auto const& dNdx_u = ip_data.dNdx_u;

            ParameterLib::SpatialPosition const x_position{
                std::nullopt, _element.getID(), ip,
                MathLib::Point3d(
                    NumLib::interpolateCoordinates<
                        ShapeFunctionDisplacement,
                        ShapeMatricesTypeDisplacement>(_element, N_u))};

            updateConstitutiveRelations(local_x, local_x_prev, x_position, t,
                                        dt, _ip_data[ip], _ip_data_output[ip]);

            auto const x_coord =
                NumLib::interpolateXCoordinate<ShapeFunctionDisplacement,
                                               ShapeMatricesTypeDisplacement>(
                    _element, N_u);
            auto const B = LinearBMatrix::computeBMatrix<
                DisplacementDim, ShapeFunctionDisplacement::NPOINTS,
                typename BMatricesType::BMatrixType>(dNdx_u, N_u, x_coord,
                                                     _is_axially_symmetric);

            ConstitutiveRelationsValues<DisplacementDim> crv;

            MathLib::KelvinVector::KelvinVectorType<DisplacementDim> const
                eps_prev = B * u_prev;

            _ip_data[ip].eps0 =
                _ip_data[ip].eps0_prev +
                (1 - _ip_data[ip].phi_fr_prev / _ip_data[ip].porosity) *
                    (eps_prev - _ip_data[ip].eps0_prev);
            _ip_data[ip].pushBackState();
        }
    }

    void computeSecondaryVariableConcrete(
        double const t, double const dt, Eigen::VectorXd const& local_x,
        Eigen::VectorXd const& local_x_prev) override;

    Eigen::Map<const Eigen::RowVectorXd> getShapeMatrix(
        const unsigned integration_point) const override
    {
        auto const& N_u = _secondary_data.N_u[integration_point];

        // assumes N is stored contiguously in memory
        return Eigen::Map<const Eigen::RowVectorXd>(N_u.data(), N_u.size());
    }

    std::vector<double> const& getIntPtDarcyVelocity(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const override;

    // TODO (naumov) This method is same as getIntPtSigma but for arguments and
    // the ordering of the cache_mat.
    // There should be only one.
    std::vector<double> getSigma() const override
    {
        constexpr int kelvin_vector_size =
            MathLib::KelvinVector::kelvin_vector_dimensions(DisplacementDim);

        return transposeInPlace<kelvin_vector_size>(
            [this](std::vector<double>& values)
            { return getIntPtSigma(0, {}, {}, values); });
    }

    std::vector<double> const& getIntPtFluidDensity(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const override;

    std::vector<double> const& getIntPtViscosity(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const override;

private:
    using BMatricesType =
        BMatrixPolicyType<ShapeFunctionDisplacement, DisplacementDim>;
    using IpData =
        IntegrationPointData<BMatricesType, ShapeMatricesTypeDisplacement,
                             ShapeMatricesTypePressure, DisplacementDim,
                             ShapeFunctionDisplacement::NPOINTS>;

private:
    ConstitutiveRelationsValues<DisplacementDim> updateConstitutiveRelations(
        Eigen::Ref<Eigen::VectorXd const> const local_x,
        Eigen::Ref<Eigen::VectorXd const> const local_x_prev,
        ParameterLib::SpatialPosition const& x_position, double const t,
        double const dt, IpData& ip_data,
        IntegrationPointDataForOutput<DisplacementDim>& ip_data_output) const;

    std::size_t setSigma(double const* values)
    {
        return ProcessLib::setIntegrationPointKelvinVectorData<DisplacementDim>(
            values, _ip_data, &IpData::sigma_eff);
    }

    std::vector<double> const& getIntPtSigma(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        return ProcessLib::getIntegrationPointKelvinVectorData<DisplacementDim>(
            _ip_data, &IpData::sigma_eff, cache);
    }

    std::vector<double> const& getIntPtSigmaIce(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        return ProcessLib::getIntegrationPointKelvinVectorData<DisplacementDim>(
            _ip_data, &IpData::sigma_eff_ice, cache);
    }

    std::vector<double> getEpsilonM() const override
    {
        constexpr int kelvin_vector_size =
            MathLib::KelvinVector::kelvin_vector_dimensions(DisplacementDim);

        return transposeInPlace<kelvin_vector_size>(
            [this](std::vector<double>& values)
            { return getIntPtEpsilonM(0, {}, {}, values); });
    }

    virtual std::vector<double> const& getIntPtEpsilonM(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        return ProcessLib::getIntegrationPointKelvinVectorData<DisplacementDim>(
            _ip_data, &IpData::eps_m, cache);
    }

    std::vector<double> getEpsilon() const override
    {
        constexpr int kelvin_vector_size =
            MathLib::KelvinVector::kelvin_vector_dimensions(DisplacementDim);

        return transposeInPlace<kelvin_vector_size>(
            [this](std::vector<double>& values)
            { return getIntPtEpsilon(0, {}, {}, values); });
    }

    virtual std::vector<double> const& getIntPtEpsilon(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        return ProcessLib::getIntegrationPointKelvinVectorData<DisplacementDim>(
            _ip_data, &IpData::eps, cache);
    }

    std::vector<double> const& getIntPtIceVolume(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        return ProcessLib::getIntegrationPointScalarData(
            _ip_data, &IpData::phi_fr, cache);
    }

    unsigned getNumberOfIntegrationPoints() const override
    {
        return _integration_method.getNumberOfPoints();
    }

    int getMaterialID() const override
    {
        return _process_data.material_ids == nullptr
                   ? 0
                   : (*_process_data.material_ids)[_element.getID()];
    }

    std::vector<double> getMaterialStateVariableInternalState(
        std::function<std::span<double>(
            typename MaterialLib::Solids::MechanicsBase<DisplacementDim>::
                MaterialStateVariables&)> const& get_values_span,
        int const& n_components) const override
    {
        return ProcessLib::getIntegrationPointDataMaterialStateVariables(
            _ip_data, &IpData::material_state_variables, get_values_span,
            n_components);
    }

    typename MaterialLib::Solids::MechanicsBase<
        DisplacementDim>::MaterialStateVariables const&
    getMaterialStateVariablesAt(unsigned integration_point) const override
    {
        return *_ip_data[integration_point].material_state_variables;
    }

private:
    template <typename SolutionVector>
    static constexpr auto localDOF(SolutionVector const& x)
    {
        return NumLib::localDOF<
            ShapeFunctionPressure, ShapeFunctionPressure,
            NumLib::Vectorial<ShapeFunctionDisplacement, DisplacementDim>>(x);
    }

    ThermoHydroMechanicsProcessData<DisplacementDim>& _process_data;

    std::vector<IpData, Eigen::aligned_allocator<IpData>> _ip_data;
    std::vector<IntegrationPointDataForOutput<DisplacementDim>,
                Eigen::aligned_allocator<
                    IntegrationPointDataForOutput<DisplacementDim>>>
        _ip_data_output;

    NumLib::GenericIntegrationMethod const& _integration_method;
    MeshLib::Element const& _element;
    bool const _is_axially_symmetric;
    SecondaryData<
        typename ShapeMatricesTypeDisplacement::ShapeMatrices::ShapeType>
        _secondary_data;

    // The shape function of pressure has the same form with the shape function
    // of temperature
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
