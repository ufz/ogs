/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
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
#include "ProcessLib/TH2M/ConstitutiveRelations/Porosity.h"
#include "ProcessLib/TH2M/ConstitutiveRelations/TransportPorosity.h"
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
                              std::vector<double>& local_rhs_data,
                              std::vector<double>& local_Jac_data) override;

    void initializeConcrete() override
    {
        unsigned const n_integration_points =
            this->integration_method_.getNumberOfPoints();
        auto const time_independent = std::numeric_limits<double>::quiet_NaN();
        auto const& medium =
            *this->process_data_.media_map.getMedium(this->element_.getID());

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            auto& ip_data = _ip_data[ip];

            ParameterLib::SpatialPosition const x_position{
                std::nullopt, this->element_.getID(),
                MathLib::Point3d(NumLib::interpolateCoordinates<
                                 ShapeFunctionDisplacement,
                                 ShapeMatricesTypeDisplacement>(this->element_,
                                                                ip_data.N_u))};
            auto& current_state = this->current_states_[ip];

            /// Set initial stress from parameter.
            if (this->process_data_.initial_stress.value)
            {
                current_state.eff_stress_data.sigma_eff.noalias() =
                    MathLib::KelvinVector::symmetricTensorToKelvinVector<
                        DisplacementDim>(
                        (*this->process_data_.initial_stress.value)(
                            std::numeric_limits<
                                double>::quiet_NaN() /* time independent */,
                            x_position));
            }

            if (*this->process_data_.initialize_porosity_from_medium_property)
            {
                // Initial porosity. Could be read from integration point data
                // or mesh.
                current_state.porosity_data.phi =
                    // std::get<iConstitutiveRelations::PorosityData>(current_state).phi
                    // =
                    medium.property(MaterialPropertyLib::porosity)
                        .template initialValue<double>(x_position,
                                                       time_independent);

                if (medium.hasProperty(MaterialPropertyLib::PropertyType::transport_porosity))
                {
                    current_state.transport_porosity_data.phi =
                        medium.property(MaterialPropertyLib::transport_porosity)
                            .template initialValue<double>(x_position,
                                                           time_independent);
                }
                else
                {
                    current_state.transport_porosity_data.phi =
                        current_state.porosity_data.phi;
                }
            }

            double const t = 0;  // TODO (naumov) pass t from top
            auto& material_state = this->material_states_[ip];
            this->solid_material_.initializeInternalStateVariables(
                t, x_position, *material_state.material_state_variables);

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

    std::tuple<
        std::vector<ConstitutiveRelations::ConstitutiveData<DisplacementDim>>,
        std::vector<
            ConstitutiveRelations::ConstitutiveTempData<DisplacementDim>>>
    updateConstitutiveVariables(
        Eigen::VectorXd const& local_x, Eigen::VectorXd const& local_x_prev,
        double const t, double const dt,
        ConstitutiveRelations::ConstitutiveModels<DisplacementDim> const&
            models);

    std::vector<ConstitutiveRelations::DerivativesData<DisplacementDim>>
    updateConstitutiveVariablesDerivatives(
        Eigen::VectorXd const& local_x, Eigen::VectorXd const& local_x_prev,
        double const t, double const dt,
        std::vector<
            ConstitutiveRelations::ConstitutiveData<DisplacementDim>> const&
            ip_constitutive_data,
        std::vector<
            ConstitutiveRelations::ConstitutiveTempData<DisplacementDim>> const&
            ip_constitutive_variables,
        ConstitutiveRelations::ConstitutiveModels<DisplacementDim> const&
            models);

    virtual std::optional<VectorSegment> getVectorDeformationSegment()
        const override
    {
        return std::optional<VectorSegment>{
            {displacement_index, displacement_size}};
    }

private:
    using BMatricesType =
        BMatrixPolicyType<ShapeFunctionDisplacement, DisplacementDim>;
    using IpData = IntegrationPointData<ShapeMatricesTypeDisplacement,
                                        ShapeMatricesTypePressure>;
    std::vector<IpData> _ip_data;

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
