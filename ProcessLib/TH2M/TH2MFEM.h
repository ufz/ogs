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

#include "ConstitutiveRelations/ConstitutiveModels.h"
#include "IntegrationPointData.h"
#include "LocalAssemblerInterface.h"
#include "MaterialLib/PhysicalConstant.h"
#include "MaterialLib/SolidModels/LinearElasticIsotropic.h"
#include "MathLib/EigenBlockMatrixView.h"
#include "MathLib/KelvinVector.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "NumLib/DOF/DOFTableUtil.h"
#include "NumLib/Fem/InitShapeMatrices.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "ParameterLib/Parameter.h"
#include "ProcessLib/Deformation/BMatrixPolicy.h"
#include "ProcessLib/Deformation/LinearBMatrix.h"
#include "ProcessLib/LocalAssemblerTraits.h"
#include "ProcessLib/Utils/SetOrGetIntegrationPointData.h"
#include "ProcessLib/Utils/TransposeInPlace.h"

namespace ProcessLib
{
namespace TH2M
{
/// Used by for extrapolation of the integration point values. It is ordered
/// (and stored) by integration points.
template <typename ShapeMatrixType>
struct SecondaryData
{
    std::vector<ShapeMatrixType, Eigen::aligned_allocator<ShapeMatrixType>> N_u;
};

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
class TH2MLocalAssembler : public LocalAssemblerInterface<DisplacementDim>
{
public:
    using ShapeMatricesTypeDisplacement =
        ShapeMatrixPolicyType<ShapeFunctionDisplacement, DisplacementDim>;

    using ShapeMatricesTypePressure =
        ShapeMatrixPolicyType<ShapeFunctionPressure, DisplacementDim>;

    template <int N>
    using VectorType =
        typename ShapeMatricesTypePressure::template VectorType<N>;

    template <int M, int N>
    using MatrixType =
        typename ShapeMatricesTypePressure::template MatrixType<M, N>;

    using GlobalDimMatrixType =
        typename ShapeMatricesTypePressure::GlobalDimMatrixType;

    using GlobalDimVectorType =
        typename ShapeMatricesTypePressure::GlobalDimVectorType;

    static int const KelvinVectorSize =
        MathLib::KelvinVector::kelvin_vector_dimensions(DisplacementDim);
    using SymmetricTensor = Eigen::Matrix<double, KelvinVectorSize, 1>;

    static constexpr auto& N_u_op = MathLib::eigenBlockMatrixView<
        DisplacementDim,
        typename ShapeMatricesTypeDisplacement::NodalRowVectorType>;

    using Invariants = MathLib::KelvinVector::Invariants<KelvinVectorSize>;

    TH2MLocalAssembler(TH2MLocalAssembler const&) = delete;
    TH2MLocalAssembler(TH2MLocalAssembler&&) = delete;

    TH2MLocalAssembler(
        MeshLib::Element const& e,
        std::size_t const /*local_matrix_size*/,
        NumLib::GenericIntegrationMethod const& integration_method,
        bool const is_axially_symmetric,
        TH2MProcessData<DisplacementDim>& process_data);

private:
    /// \return the number of read integration points.
    std::size_t setIPDataInitialConditions(
        std::string_view const name,
        double const* values,
        int const integration_order) override;

    void setInitialConditionsConcrete(Eigen::VectorXd const local_x,
                                      double const t,
                                      int const process_id) override;

    void assemble(double const /*t*/, double const /*dt*/,
                  std::vector<double> const& /*local_x*/,
                  std::vector<double> const& /*local_x_prev*/,
                  std::vector<double>& /*local_M_data*/,
                  std::vector<double>& /*local_K_data*/,
                  std::vector<double>& /*local_rhs_data*/) override;

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
            this->integration_method_.getNumberOfPoints();

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            auto& ip_data = _ip_data[ip];

            ParameterLib::SpatialPosition const x_position{
                std::nullopt, this->element_.getID(), ip,
                MathLib::Point3d(NumLib::interpolateCoordinates<
                                 ShapeFunctionDisplacement,
                                 ShapeMatricesTypeDisplacement>(this->element_,
                                                                ip_data.N_u))};

            /// Set initial stress from parameter.
            if (this->process_data_.initial_stress.value)
            {
                this->current_states_[ip].eff_stress_data.sigma.noalias() =
                    MathLib::KelvinVector::symmetricTensorToKelvinVector<
                        DisplacementDim>(
                        (*this->process_data_.initial_stress.value)(
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

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            this->prev_states_[ip] = this->current_states_[ip];
        }
    }

    void postTimestepConcrete(Eigen::VectorXd const& /*local_x*/,
                              Eigen::VectorXd const& /*local_x_prev*/,
                              double const /*t*/, double const /*dt*/,
                              int const /*process_id*/) override
    {
        unsigned const n_integration_points =
            this->integration_method_.getNumberOfPoints();

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            _ip_data[ip].pushBackState();
            this->material_states_[ip].pushBackState();
        }

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            this->prev_states_[ip] = this->current_states_[ip];
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

    std::vector<double> const& getIntPtDarcyVelocityGas(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const override;
    std::vector<double> const& getIntPtDarcyVelocityLiquid(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const override;
    std::vector<double> const& getIntPtDiffusionVelocityVapourGas(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const override;
    std::vector<double> const& getIntPtDiffusionVelocityGasGas(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const override;
    std::vector<double> const& getIntPtDiffusionVelocitySoluteLiquid(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const override;
    std::vector<double> const& getIntPtDiffusionVelocityLiquidLiquid(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const override;

    std::tuple<
        std::vector<ConstitutiveRelations::ConstitutiveData<DisplacementDim>>,
        std::vector<
            ConstitutiveRelations::ConstitutiveTempData<DisplacementDim>>>
    updateConstitutiveVariables(
        Eigen::VectorXd const& local_x, Eigen::VectorXd const& local_x_prev,
        double const t, double const dt,
        ConstitutiveRelations::ConstitutiveModels<DisplacementDim> const&
            models);

    // TODO: Here is some refactoring potential. All secondary variables could
    // be stored in some container to avoid defining one method for each
    // variable.

    virtual std::vector<double> const& getIntPtLiquidDensity(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        return ProcessLib::getIntegrationPointScalarData(_ip_data,
                                                         &IpData::rhoLR, cache);
    }

    virtual std::vector<double> const& getIntPtGasDensity(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        return ProcessLib::getIntegrationPointScalarData(_ip_data,
                                                         &IpData::rhoGR, cache);
    }

    virtual std::vector<double> const& getIntPtSolidDensity(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        return ProcessLib::getIntegrationPointScalarData(_ip_data,
                                                         &IpData::rhoSR, cache);
    }

    virtual std::vector<double> const& getIntPtVapourPressure(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        return ProcessLib::getIntegrationPointScalarData(_ip_data,
                                                         &IpData::pWGR, cache);
    }

    virtual std::vector<double> const& getIntPtPorosity(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        return ProcessLib::getIntegrationPointScalarData(_ip_data, &IpData::phi,
                                                         cache);
    }

    virtual std::vector<double> const& getIntPtMoleFractionGas(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        return ProcessLib::getIntegrationPointScalarData(_ip_data,
                                                         &IpData::xnCG, cache);
    }
    virtual std::vector<double> const& getIntPtMassFractionGas(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        return ProcessLib::getIntegrationPointScalarData(_ip_data,
                                                         &IpData::xmCG, cache);
    }
    virtual std::vector<double> const& getIntPtMassFractionLiquid(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        return ProcessLib::getIntegrationPointScalarData(_ip_data,
                                                         &IpData::xmWL, cache);
    }

    virtual std::vector<double> const& getIntPtEnthalpyGas(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        return ProcessLib::getIntegrationPointScalarData(_ip_data, &IpData::h_G,
                                                         cache);
    }
    virtual std::vector<double> const& getIntPtEnthalpyLiquid(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        return ProcessLib::getIntegrationPointScalarData(_ip_data, &IpData::h_L,
                                                         cache);
    }
    virtual std::vector<double> const& getIntPtEnthalpySolid(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        return ProcessLib::getIntegrationPointScalarData(_ip_data, &IpData::h_S,
                                                         cache);
    }

private:
    using BMatricesType =
        BMatrixPolicyType<ShapeFunctionDisplacement, DisplacementDim>;
    using IpData =
        IntegrationPointData<BMatricesType, ShapeMatricesTypeDisplacement,
                             ShapeMatricesTypePressure, DisplacementDim,
                             ShapeFunctionDisplacement::NPOINTS>;
    std::vector<IpData, Eigen::aligned_allocator<IpData>> _ip_data;

    SecondaryData<
        typename ShapeMatricesTypeDisplacement::ShapeMatrices::ShapeType>
        _secondary_data;

    // The shape function of pressure has the same form with the shape function
    // of temperature
    static const int gas_pressure_index = 0;
    static const int gas_pressure_size = ShapeFunctionPressure::NPOINTS;
    static const int capillary_pressure_index = ShapeFunctionPressure::NPOINTS;
    static const int capillary_pressure_size = ShapeFunctionPressure::NPOINTS;
    static const int temperature_index = 2 * ShapeFunctionPressure::NPOINTS;
    static const int temperature_size = ShapeFunctionPressure::NPOINTS;
    static const int displacement_index = ShapeFunctionPressure::NPOINTS * 3;
    static const int displacement_size =
        ShapeFunctionDisplacement::NPOINTS * DisplacementDim;

    static const int C_index = 0;
    static const int C_size = ShapeFunctionPressure::NPOINTS;
    static const int W_index = ShapeFunctionPressure::NPOINTS;
    static const int W_size = ShapeFunctionPressure::NPOINTS;
};

}  // namespace TH2M
}  // namespace ProcessLib

#include "TH2MFEM-impl.h"
