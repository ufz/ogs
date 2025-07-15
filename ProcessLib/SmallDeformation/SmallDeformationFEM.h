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

#include <limits>
#include <memory>
#include <optional>
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
    double integration_weight;
    typename ShapeMatricesType::NodalRowVectorType N_u;
    typename ShapeMatricesType::GlobalDimNodalMatrixType dNdx_u;

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
    using BBarMatrixType = typename BMatricesType::BBarMatrixType;
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

        ip_data_.resize(n_integration_points);
        secondary_data_.N.resize(n_integration_points);

        auto const shape_matrices =
            NumLib::initShapeMatrices<ShapeFunction, ShapeMatricesType,
                                      DisplacementDim>(
                e, is_axially_symmetric, this->integration_method_);

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            auto& ip_data = ip_data_[ip];
            auto const& sm = shape_matrices[ip];
            ip_data_[ip].integration_weight =
                this->integration_method_.getWeightedPoint(ip).getWeight() *
                sm.integralMeasure * sm.detJ;

            ip_data.N_u = sm.N;
            ip_data.dNdx_u = sm.dNdx;

            secondary_data_.N[ip] = shape_matrices[ip].N;
        }
    }

    void initializeConcrete() override
    {
        unsigned const n_integration_points =
            this->integration_method_.getNumberOfPoints();
        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            auto const& ip_data = ip_data_[ip];

            ParameterLib::SpatialPosition const x_position{
                std::nullopt, this->element_.getID(),
                MathLib::Point3d(
                    NumLib::interpolateCoordinates<ShapeFunction,
                                                   ShapeMatricesType>(
                        this->element_, ip_data.N_u))};

            /// Set initial stress from parameter.
            if (this->process_data_.initial_stress != nullptr)
            {
                this->current_states_[ip].stress_data.sigma.noalias() =
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

            material_state.pushBackState();
        }

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            this->prev_states_[ip] = this->current_states_[ip];
        }
    }

    void setInitialConditionsConcrete(Eigen::VectorXd const local_x,
                                      double const /*t*/,
                                      int const /*process_id*/) override
    {
        unsigned const n_integration_points =
            this->integration_method_.getNumberOfPoints();
        auto const B_dil_bar = getDilatationalBBarMatrix();

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            auto const& N = ip_data_[ip].N_u;
            auto const& dNdx = ip_data_[ip].dNdx_u;
            ParameterLib::SpatialPosition const x_position{
                std::nullopt, this->element_.getID(),
                MathLib::Point3d(
                    NumLib::interpolateCoordinates<ShapeFunction,
                                                   ShapeMatricesType>(
                        this->element_, N))};

            auto const x_coord =
                NumLib::interpolateXCoordinate<ShapeFunction,
                                               ShapeMatricesType>(
                    this->element_, N);
            auto const B = LinearBMatrix::computeBMatrixPossiblyWithBbar<
                DisplacementDim, ShapeFunction::NPOINTS, BBarMatrixType,
                typename BMatricesType::BMatrixType>(
                dNdx, N, B_dil_bar, x_coord, this->is_axially_symmetric_);

            this->output_data_[ip].eps_data.eps.noalias() = B * local_x;
        }
    }

    typename ConstitutiveRelations::ConstitutiveData<DisplacementDim>
    updateConstitutiveRelations(
        BMatrixType const& B, Eigen::Ref<Eigen::VectorXd const> const& u,
        Eigen::Ref<Eigen::VectorXd const> const& u_prev,
        ParameterLib::SpatialPosition const& x_position, double const t,
        double const dt,
        typename ConstitutiveRelations::ConstitutiveSetting<DisplacementDim>&
            CS,
        MaterialPropertyLib::Medium const& medium,
        typename ConstitutiveRelations::StatefulData<DisplacementDim>&
            current_state,
        typename ConstitutiveRelations::StatefulDataPrev<DisplacementDim> const&
            prev_state,
        MaterialStateData<DisplacementDim>& material_state,
        typename ConstitutiveRelations::OutputData<DisplacementDim>&
            output_data) const
    {
        double const T_ref =
            this->process_data_.reference_temperature
                ? (*this->process_data_.reference_temperature)(t, x_position)[0]
                : std::numeric_limits<double>::quiet_NaN();

        typename ConstitutiveRelations::ConstitutiveModels<DisplacementDim>
            models{this->process_data_, this->solid_material_};
        typename ConstitutiveRelations::ConstitutiveTempData<DisplacementDim>
            tmp;
        typename ConstitutiveRelations::ConstitutiveData<DisplacementDim> CD;

        CS.eval(models, t, dt, x_position,  //
                medium,                     //
                T_ref, B * u, B * u_prev,   //
                current_state, prev_state, material_state, tmp, output_data,
                CD);

        return CD;
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

    std::optional<BBarMatrixType> getDilatationalBBarMatrix() const
    {
        if (!(this->process_data_.use_b_bar))
        {
            return std::nullopt;
        }

        return LinearBMatrix::computeDilatationalBbar<
            DisplacementDim, ShapeFunction::NPOINTS, ShapeFunction,
            BBarMatrixType, ShapeMatricesType, IpData>(
            ip_data_, this->element_, this->integration_method_,
            this->is_axially_symmetric_);
    }

    void assembleWithJacobian(double const t, double const dt,
                              std::vector<double> const& local_x,
                              std::vector<double> const& local_x_prev,
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

        typename ConstitutiveRelations::ConstitutiveSetting<DisplacementDim>
            constitutive_setting;
        auto const& medium =
            *this->process_data_.media_map.getMedium(this->element_.getID());

        auto const B_dil_bar = getDilatationalBBarMatrix();

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            auto const& w = ip_data_[ip].integration_weight;
            auto const& N = ip_data_[ip].N_u;
            auto const& dNdx = ip_data_[ip].dNdx_u;

            ParameterLib::SpatialPosition const x_position{
                std::nullopt, this->element_.getID(),
                MathLib::Point3d(
                    NumLib::interpolateCoordinates<ShapeFunction,
                                                   ShapeMatricesType>(
                        this->element_, N))};

            auto const x_coord =
                NumLib::interpolateXCoordinate<ShapeFunction,
                                               ShapeMatricesType>(
                    this->element_, N);
            auto const B = LinearBMatrix::computeBMatrixPossiblyWithBbar<
                DisplacementDim, ShapeFunction::NPOINTS, BBarMatrixType,
                typename BMatricesType::BMatrixType>(
                dNdx, N, B_dil_bar, x_coord, this->is_axially_symmetric_);

            auto const CD = updateConstitutiveRelations(
                B, u, u_prev, x_position, t, dt, constitutive_setting, medium,
                this->current_states_[ip], this->prev_states_[ip],
                this->material_states_[ip], this->output_data_[ip]);

            auto const& sigma = this->current_states_[ip].stress_data.sigma;
            auto const& b = *CD.volumetric_body_force;
            auto const& C = CD.s_mech_data.stiffness_tensor;

            local_b.noalias() -=
                (B.transpose() * sigma - N_u_op(N).transpose() * b) * w;
            local_Jac.noalias() += B.transpose() * C * B * w;
        }
    }

    void postTimestepConcrete(Eigen::VectorXd const& local_x,
                              Eigen::VectorXd const& local_x_prev,
                              double const t, double const dt,
                              int const /*process_id*/) override
    {
        unsigned const n_integration_points =
            this->integration_method_.getNumberOfPoints();

        typename ConstitutiveRelations::ConstitutiveSetting<DisplacementDim>
            constitutive_setting;

        auto const& medium =
            *this->process_data_.media_map.getMedium(this->element_.getID());

        auto const B_dil_bar = getDilatationalBBarMatrix();

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            auto const& N = ip_data_[ip].N_u;
            auto const& dNdx = ip_data_[ip].dNdx_u;

            auto const x_coord =
                NumLib::interpolateXCoordinate<ShapeFunction,
                                               ShapeMatricesType>(
                    this->element_, N);

            ParameterLib::SpatialPosition const x_position{
                std::nullopt, this->element_.getID(),
                MathLib::Point3d(
                    NumLib::interpolateCoordinates<ShapeFunction,
                                                   ShapeMatricesType>(
                        this->element_, N))};

            auto const B = LinearBMatrix::computeBMatrixPossiblyWithBbar<
                DisplacementDim, ShapeFunction::NPOINTS, BBarMatrixType,
                typename BMatricesType::BMatrixType>(
                dNdx, N, B_dil_bar, x_coord, this->is_axially_symmetric_);

            updateConstitutiveRelations(
                B, local_x, local_x_prev, x_position, t, dt,
                constitutive_setting, medium, this->current_states_[ip],
                this->prev_states_[ip], this->material_states_[ip],
                this->output_data_[ip]);

            this->material_states_[ip].pushBackState();
        }

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            this->prev_states_[ip] = this->current_states_[ip];
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
                                this->integration_method_, ip_data_,
                                this->current_states_, this->output_data_,
                                this->element_, this->is_axially_symmetric_);
    }

    Eigen::Map<const Eigen::RowVectorXd> getShapeMatrix(
        const unsigned integration_point) const override
    {
        auto const& N = secondary_data_.N[integration_point];

        // assumes N is stored contiguously in memory
        return Eigen::Map<const Eigen::RowVectorXd>(N.data(), N.size());
    }

private:
    static constexpr auto localDOF(std::vector<double> const& x)
    {
        return NumLib::localDOF<
            NumLib::Vectorial<ShapeFunction, DisplacementDim>>(x);
    }

private:
    std::vector<IpData, Eigen::aligned_allocator<IpData>> ip_data_;

    SecondaryData<typename ShapeMatrices::ShapeType> secondary_data_;
};

}  // namespace SmallDeformation
}  // namespace ProcessLib
