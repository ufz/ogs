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

#include <limits>
#include <memory>
#include <vector>

#include "IntegrationPointData.h"
#include "LocalAssemblerInterface.h"
#include "MaterialLib/PhysicalConstant.h"
#include "MathLib/EigenBlockMatrixView.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "MathLib/VectorizedTensor.h"
#include "NumLib/DOF/LocalDOF.h"
#include "NumLib/Extrapolation/ExtrapolatableElement.h"
#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/InitShapeMatrices.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "ParameterLib/Parameter.h"
#include "ProcessLib/Deformation/BMatrixPolicy.h"
#include "ProcessLib/Deformation/GMatrix.h"
#include "ProcessLib/Deformation/GMatrixPolicy.h"
#include "ProcessLib/Deformation/LinearBMatrix.h"
#include "ProcessLib/Deformation/NonLinearBMatrix.h"
#include "ProcessLib/LocalAssemblerInterface.h"
#include "ProcessLib/LocalAssemblerTraits.h"

namespace ProcessLib
{
namespace LargeDeformation
{
namespace MPL = MaterialPropertyLib;

/// Used for the extrapolation of the integration point values. It is ordered
/// (and stored) by integration points.
template <typename ShapeMatrixType>
struct SecondaryData
{
    std::vector<ShapeMatrixType, Eigen::aligned_allocator<ShapeMatrixType>> N;
};

template <int DisplacementDim, typename ShapeMatricesType>
Eigen::Matrix<double, MPL::tensorSize(DisplacementDim),
              MPL::tensorSize(DisplacementDim)>
computeSigmaGeom(Eigen::Matrix3d const& sigma_tensor)
{
    static constexpr auto& sigma_geom_op = MathLib::eigenBlockMatrixView<
        DisplacementDim,
        Eigen::Matrix<double, DisplacementDim, DisplacementDim>>;

    using SigmaGeom = Eigen::Matrix<double, MPL::tensorSize(DisplacementDim),
                                    MPL::tensorSize(DisplacementDim)>;
    if constexpr (DisplacementDim == 2)
    {
        SigmaGeom sigma_geom = SigmaGeom::Zero(5, 5);
        sigma_geom.template block<4, 4>(0, 0) =
            sigma_geom_op(sigma_tensor.template block<2, 2>(0, 0).eval());
        sigma_geom(4, 4) = sigma_tensor(2, 2);

        return sigma_geom;
    }

    if constexpr (DisplacementDim == 3)
    {
        return sigma_geom_op(sigma_tensor);
    }
}

template <typename ShapeFunction, int DisplacementDim>
class LargeDeformationLocalAssembler
    : public LargeDeformationLocalAssemblerInterface<DisplacementDim>
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

    LargeDeformationLocalAssembler(LargeDeformationLocalAssembler const&) =
        delete;
    LargeDeformationLocalAssembler(LargeDeformationLocalAssembler&&) = delete;

    LargeDeformationLocalAssembler(
        MeshLib::Element const& e,
        std::size_t const /*local_matrix_size*/,
        NumLib::GenericIntegrationMethod const& integration_method,
        bool const is_axially_symmetric,
        LargeDeformationProcessData<DisplacementDim>& process_data)
        : LargeDeformationLocalAssemblerInterface<DisplacementDim>(
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

            _secondary_data.N[ip] = shape_matrices[ip].N;
        }
    }

    void initializeConcrete() override
    {
        unsigned const n_integration_points =
            this->integration_method_.getNumberOfPoints();
        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            auto const& ip_data = _ip_data[ip];

            ParameterLib::SpatialPosition const x_position{
                std::nullopt, this->element_.getID(), ip,
                MathLib::Point3d(
                    NumLib::interpolateCoordinates<ShapeFunction,
                                                   ShapeMatricesType>(
                        this->element_, ip_data.N))};

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
            this->prev_states_[ip] = this->current_states_[ip];
        }
    }

    typename ConstitutiveRelations::ConstitutiveData<DisplacementDim>
    updateConstitutiveRelations(
        Eigen::Ref<Eigen::VectorXd const> const& u,
        Eigen::Ref<Eigen::VectorXd const> const& u_prev,
        ParameterLib::SpatialPosition const& x_position, double const t,
        double const dt,
        IntegrationPointData<BMatricesType, ShapeMatricesType, DisplacementDim>&
            ip_data,
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
        auto const& N = ip_data.N;
        auto const& dNdx = ip_data.dNdx;
        auto const x_coord =
            NumLib::interpolateXCoordinate<ShapeFunction, ShapeMatricesType>(
                this->element_, N);

        // For the 2D case the 33-component is needed (and the four entries
        // of the non-symmetric matrix); In 3d there are nine entries.
        GradientMatrixType G(
            DisplacementDim * DisplacementDim + (DisplacementDim == 2 ? 1 : 0),
            ShapeFunction::NPOINTS * DisplacementDim);
        Deformation::computeGMatrix<DisplacementDim, ShapeFunction::NPOINTS>(
            dNdx, G, this->is_axially_symmetric_, N, x_coord);

        GradientVectorType const grad_u = G * u;

        auto const B_0 =
            LinearBMatrix::computeBMatrix<DisplacementDim,
                                          ShapeFunction::NPOINTS,
                                          typename BMatricesType::BMatrixType>(
                dNdx, N, x_coord, this->is_axially_symmetric_);

        auto const B_N = NonLinearBMatrix::computeBMatrix<
            DisplacementDim, ShapeFunction::NPOINTS,
            typename BMatricesType::BMatrixType>(dNdx, N, grad_u, x_coord,
                                                 this->is_axially_symmetric_);

        auto const B = B_0 + 0.5 * B_N;  // For Green-Lagrange strain.

        double const T_ref =
            this->process_data_.reference_temperature
                ? (*this->process_data_.reference_temperature)(t, x_position)[0]
                : std::numeric_limits<double>::quiet_NaN();

        typename ConstitutiveRelations::ConstitutiveModels<DisplacementDim>
            models(this->process_data_, this->solid_material_);
        typename ConstitutiveRelations::ConstitutiveTempData<DisplacementDim>
            tmp;
        typename ConstitutiveRelations::ConstitutiveData<DisplacementDim> CD;

        output_data.eps_data.eps = B * u;
        output_data.deformation_gradient_data.deformation_gradient =
            grad_u + MathLib::VectorizedTensor::identity<DisplacementDim>();
        output_data.deformation_gradient_data.volume_ratio =
            MathLib::VectorizedTensor::determinant(
                output_data.deformation_gradient_data.deformation_gradient);

        CS.eval(
            models, t, dt, x_position,              //
            medium,                                 //
            T_ref,                                  //
            output_data.deformation_gradient_data,  //
            G * u_prev + MathLib::VectorizedTensor::identity<DisplacementDim>(),
            current_state, prev_state, material_state, tmp, CD);

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
            "LargeDeformationLocalAssembler: assembly without jacobian is not "
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

        typename ConstitutiveRelations::ConstitutiveSetting<DisplacementDim>
            constitutive_setting;
        auto const& medium =
            *this->process_data_.media_map.getMedium(this->element_.getID());

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

            // For the 2D case the 33-component is needed (and the four entries
            // of the non-symmetric matrix); In 3d there are nine entries.
            GradientMatrixType G(DisplacementDim * DisplacementDim +
                                     (DisplacementDim == 2 ? 1 : 0),
                                 ShapeFunction::NPOINTS * DisplacementDim);
            Deformation::computeGMatrix<DisplacementDim,
                                        ShapeFunction::NPOINTS>(
                dNdx, G, this->is_axially_symmetric_, N, x_coord);

            GradientVectorType const grad_u = G * u;

            auto const B_0 = LinearBMatrix::computeBMatrix<
                DisplacementDim, ShapeFunction::NPOINTS,
                typename BMatricesType::BMatrixType>(
                dNdx, N, x_coord, this->is_axially_symmetric_);

            auto const B_N = NonLinearBMatrix::computeBMatrix<
                DisplacementDim, ShapeFunction::NPOINTS,
                typename BMatricesType::BMatrixType>(
                dNdx, N, grad_u, x_coord, this->is_axially_symmetric_);

            auto const B = B_0 + B_N;

            auto const CD = updateConstitutiveRelations(
                u, u_prev, x_position, t, dt, _ip_data[ip],
                constitutive_setting, medium, this->current_states_[ip],
                this->prev_states_[ip], this->material_states_[ip],
                this->output_data_[ip]);

            auto const& sigma = this->current_states_[ip].stress_data.sigma;
            auto const& b = *CD.volumetric_body_force;
            auto const& C = CD.s_mech_data.stiffness_tensor;

            local_b.noalias() -=
                (B.transpose() * sigma - N_u_op(N).transpose() * b) * w;

            auto const sigma_geom =
                computeSigmaGeom<DisplacementDim, ShapeMatricesType>(
                    MathLib::KelvinVector::kelvinVectorToTensor(sigma));

            local_Jac.noalias() += G.transpose() * sigma_geom * G * w;

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

        ParameterLib::SpatialPosition x_position;
        x_position.setElementID(this->element_.getID());

        typename ConstitutiveRelations::ConstitutiveSetting<DisplacementDim>
            constitutive_setting;

        auto& medium =
            *this->process_data_.media_map.getMedium(this->element_.getID());

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            x_position.setIntegrationPoint(ip);

            updateConstitutiveRelations(
                local_x, local_x_prev, x_position, t, dt, _ip_data[ip],
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

    Eigen::Map<const Eigen::RowVectorXd> getShapeMatrix(
        const unsigned integration_point) const override
    {
        auto const& N = _secondary_data.N[integration_point];

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
    std::vector<IpData, Eigen::aligned_allocator<IpData>> _ip_data;

    SecondaryData<typename ShapeMatrices::ShapeType> _secondary_data;

    static const int displacement_size =
        ShapeFunction::NPOINTS * DisplacementDim;
};

}  // namespace LargeDeformation
}  // namespace ProcessLib
