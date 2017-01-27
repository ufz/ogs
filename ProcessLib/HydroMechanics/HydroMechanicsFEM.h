/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <iostream>
#include <memory>
#include <vector>

#include "MaterialLib/SolidModels/KelvinVector.h"
#include "MaterialLib/SolidModels/LinearElasticIsotropic.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "ProcessLib/Deformation/BMatrixPolicy.h"
#include "ProcessLib/Deformation/LinearBMatrix.h"
#include "ProcessLib/LocalAssemblerInterface.h"
#include "ProcessLib/LocalAssemblerTraits.h"
#include "ProcessLib/Parameter/Parameter.h"
#include "ProcessLib/Utils/InitShapeMatrices.h"

#include "HydroMechanicsProcessData.h"

namespace ProcessLib
{
namespace HydroMechanics
{
template <typename BMatricesType, typename ShapeMatrixTypeDisplacement,
          typename ShapeMatricesTypePressure, int DisplacementDim, int NPoints>
struct IntegrationPointData final
{
    explicit IntegrationPointData(
        MaterialLib::Solids::MechanicsBase<DisplacementDim>& solid_material)
        : solid_material(solid_material),
          material_state_variables(
              solid_material.createMaterialStateVariables())
    {
    }

#if defined(_MSC_VER) && _MSC_VER < 1900
    // The default generated move-ctor is correctly generated for other
    // compilers.
    explicit IntegrationPointData(IntegrationPointData&& other)
        : b_matrices(std::move(other.b_matrices)),
          sigma_eff(std::move(other.sigma_eff)),
          sigma_eff_prev(std::move(other.sigma_eff_prev)),
          eps(std::move(other.eps)),
          eps_prev(std::move(other.eps_prev)),
          solid_material(other.solid_material),
          material_state_variables(std::move(other.material_state_variables)),
          C(std::move(other.C)),
          integration_weight(std::move(other.integration_weight))
    {
    }
#endif  // _MSC_VER

    typename ShapeMatrixTypeDisplacement::template MatrixType<
        DisplacementDim, NPoints * DisplacementDim>
        N_u;
    //typename ShapeMatrixTypeDisplacement::NodalRowVectorType N_u;
    typename BMatricesType::BMatrixType b_matrices;
    typename BMatricesType::KelvinVectorType sigma_eff, sigma_eff_prev;
    typename BMatricesType::KelvinVectorType eps, eps_prev;

    typename ShapeMatricesTypePressure::NodalRowVectorType N_p;
    typename ShapeMatricesTypePressure::GlobalDimNodalMatrixType dNdx_p;

    MaterialLib::Solids::MechanicsBase<DisplacementDim>& solid_material;
    std::unique_ptr<typename MaterialLib::Solids::MechanicsBase<
        DisplacementDim>::MaterialStateVariables>
        material_state_variables;

    typename BMatricesType::KelvinMatrixType C;
    double integration_weight;

    void pushBackState()
    {
        eps_prev = eps;
        sigma_eff_prev = sigma_eff;
        material_state_variables->pushBackState();
    }

    template <typename DisplacementVectorType>
    void updateConstitutiveRelation(double const t,
                                    SpatialPosition const& x_position,
                                    double const dt,
                                    DisplacementVectorType const& u)
    {
        eps.noalias() = b_matrices * u;
        solid_material.computeConstitutiveRelation(
            t, x_position, dt, eps_prev, eps, sigma_eff_prev, sigma_eff, C,
            *material_state_variables);
    }

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
class HydroMechanicsLocalAssembler : public ProcessLib::LocalAssemblerInterface
{
public:
    using ShapeMatricesTypeDisplacement =
        ShapeMatrixPolicyType<ShapeFunctionDisplacement, DisplacementDim>;

    // Types for pressure.
    using ShapeMatricesTypePressure =
        ShapeMatrixPolicyType<ShapeFunctionPressure, DisplacementDim>;

    HydroMechanicsLocalAssembler(HydroMechanicsLocalAssembler const&) = delete;
    HydroMechanicsLocalAssembler(HydroMechanicsLocalAssembler&&) = delete;

    HydroMechanicsLocalAssembler(
        MeshLib::Element const& e,
        std::size_t const /*local_matrix_size*/,
        bool is_axially_symmetric,
        unsigned const integration_order,
        HydroMechanicsProcessData<DisplacementDim>& process_data)
        : _process_data(process_data),
          _integration_method(integration_order),
          _element(e)
    {
        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        _ip_data.reserve(n_integration_points);

        auto const shape_matrices_u =
            initShapeMatrices<ShapeFunctionDisplacement,
                              ShapeMatricesTypeDisplacement, IntegrationMethod,
                              DisplacementDim>(e, is_axially_symmetric,
                                               _integration_method);

        auto const shape_matrices_p =
            initShapeMatrices<ShapeFunctionPressure, ShapeMatricesTypePressure,
                              IntegrationMethod, DisplacementDim>(
                e, is_axially_symmetric, _integration_method);

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            // displacement (subscript u)
            _ip_data.emplace_back(*_process_data.material);
            auto& ip_data = _ip_data[ip];
            auto const& sm = shape_matrices_u[ip];
            _ip_data[ip].integration_weight =
                _integration_method.getWeightedPoint(ip).getWeight() * sm.integralMeasure *
                shape_matrices_u[ip].detJ;
            ip_data.b_matrices.resize(
                kelvin_vector_size,
                ShapeFunctionDisplacement::NPOINTS * DisplacementDim);

            auto const x_coord =
                interpolateXCoordinate<ShapeFunctionDisplacement,
                                       ShapeMatricesTypeDisplacement>(
                    e, shape_matrices_u[ip].N);
            LinearBMatrix::computeBMatrix<DisplacementDim,
                                          ShapeFunctionDisplacement::NPOINTS>(
                shape_matrices_u[ip].dNdx, ip_data.b_matrices,
                is_axially_symmetric, shape_matrices_u[ip].N, x_coord);

            ip_data.sigma_eff.resize(kelvin_vector_size);
            ip_data.sigma_eff_prev.resize(kelvin_vector_size);
            ip_data.eps.resize(kelvin_vector_size);
            ip_data.eps_prev.resize(kelvin_vector_size);
            ip_data.C.resize(kelvin_vector_size, kelvin_vector_size);

            ip_data.N_u = ShapeMatricesTypeDisplacement::template MatrixType<
                DisplacementDim, displacement_size>::Zero(DisplacementDim,
                                                          displacement_size);
            for (int i = 0; i < DisplacementDim; ++i)
                ip_data.N_u
                    .template block<1, displacement_size / DisplacementDim>(
                        i, i * displacement_size / DisplacementDim)
                    .noalias() = shape_matrices_u[ip].N;

            ip_data.N_p = shape_matrices_p[ip].N;
            ip_data.dNdx_p = shape_matrices_p[ip].dNdx;
        }
    }

    void assemble(double const /*t*/, std::vector<double> const& /*local_x*/,
                  std::vector<double>& /*local_M_data*/,
                  std::vector<double>& /*local_K_data*/,
                  std::vector<double>& /*local_rhs_data*/) override
    {
        OGS_FATAL(
            "HydroMechanicsLocalAssembler: assembly without jacobian is not "
            "implemented.");
    }

    void assembleWithJacobian(double const t,
                              std::vector<double> const& local_x,
                              std::vector<double> const& local_xdot,
                              const double /*dxdot_dx*/, const double /*dx_dx*/,
                              std::vector<double>& /*local_M_data*/,
                              std::vector<double>& /*local_K_data*/,
                              std::vector<double>& local_rhs_data,
                              std::vector<double>& local_Jac_data) override
    {
        assert(local_x.size() == pressure_size + displacement_size);

        auto p =
            Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
                pressure_size> const>(local_x.data() + pressure_index,
                                      pressure_size);

        auto u = Eigen::Map<typename ShapeMatricesTypeDisplacement::
                                template VectorType<displacement_size> const>(
            local_x.data() + displacement_index, displacement_size);

        auto p_dot =
            Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
                pressure_size> const>(local_xdot.data() + pressure_index,
                                      pressure_size);
        auto u_dot =
            Eigen::Map<typename ShapeMatricesTypeDisplacement::
                           template VectorType<displacement_size> const>(
                local_xdot.data() + displacement_index, displacement_size);

        auto local_Jac = MathLib::createZeroedMatrix<
            typename ShapeMatricesTypeDisplacement::template MatrixType<
                displacement_size + pressure_size,
                displacement_size + pressure_size>>(
            local_Jac_data, displacement_size + pressure_size,
            displacement_size + pressure_size);

        auto local_rhs = MathLib::createZeroedVector<
            typename ShapeMatricesTypeDisplacement::template VectorType<
                displacement_size + pressure_size>>(
            local_rhs_data, displacement_size + pressure_size);

        typename ShapeMatricesTypePressure::NodalMatrixType laplace_p;
        laplace_p.setZero(pressure_size, pressure_size);

        typename ShapeMatricesTypePressure::NodalMatrixType storage_p;
        storage_p.setZero(pressure_size, pressure_size);

        typename ShapeMatricesTypeDisplacement::template MatrixType<
            displacement_size, pressure_size>
            Kup;
        Kup.setZero(displacement_size, pressure_size);

        double const& dt = _process_data.dt;

        SpatialPosition x_position;
        x_position.setElementID(_element.getID());

        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();
        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            x_position.setIntegrationPoint(ip);
            auto const& w = _ip_data[ip].integration_weight;

            auto const& N_p = _ip_data[ip].N_p;
            auto const& N_u = _ip_data[ip].N_u;
            auto const& dNdx_p = _ip_data[ip].dNdx_p;

            auto const& B = _ip_data[ip].b_matrices;
            auto const& sigma_eff = _ip_data[ip].sigma_eff;

            auto const& C = _ip_data[ip].C;

            double const S =
                _process_data.specific_storage(t, x_position)[0];
            double const K_over_mu =
                _process_data.intrinsic_permeability(t, x_position)[0] /
                _process_data.fluid_viscosity(t, x_position)[0];
            auto const alpha = _process_data.biot_coefficient(t, x_position)[0];
            auto const rho_sr = _process_data.solid_density(t, x_position)[0];
            auto const rho_fr = _process_data.fluid_density(t, x_position)[0];
            auto const porosity = _process_data.porosity(t, x_position)[0];
            auto const& b = _process_data.specific_body_force;
            auto const& identity2 = MaterialLib::SolidModels::Invariants<
                kelvin_vector_size>::identity2;

            //
            // displacement equation, displacement part
            //
            _ip_data[ip].updateConstitutiveRelation(t, x_position, dt, u);

            local_Jac
                .template block<displacement_size, displacement_size>(
                    displacement_index, displacement_index)
                .noalias() += B.transpose() * C * B * w;

            double const rho = rho_sr * (1 - porosity) + porosity * rho_fr;
            local_rhs.template segment<displacement_size>(displacement_index)
                .noalias() -=
                (B.transpose() * sigma_eff - N_u.transpose() * rho * b) * w;

            //
            // displacement equation, pressure part
            //
            Kup.noalias() += B.transpose() * alpha * identity2 * N_p * w;

            //
            // pressure equation, pressure part.
            //
            laplace_p.noalias() += dNdx_p.transpose() * K_over_mu * dNdx_p * w;

            storage_p.noalias() += N_p.transpose() * S * N_p * w;

            local_rhs.template segment<pressure_size>(pressure_index)
                .noalias() += dNdx_p.transpose() * rho_fr * K_over_mu * b * w;

            //
            // pressure equation, displacement part.
            //
            // Reusing Kup.transpose().
        }
        // displacement equation, pressure part
        local_Jac
            .template block<displacement_size, pressure_size>(
                displacement_index, pressure_index)
            .noalias() -= Kup;

        // pressure equation, pressure part.
        local_Jac
            .template block<pressure_size, pressure_size>(pressure_index,
                                                          pressure_index)
            .noalias() += laplace_p + storage_p / dt;

        // pressure equation, displacement part.
        local_Jac
            .template block<pressure_size, displacement_size>(
                pressure_index, displacement_index)
            .noalias() += Kup.transpose() / dt;

        // pressure equation
        local_rhs.template segment<pressure_size>(pressure_index)
            .noalias() -=
            laplace_p * p + storage_p * p_dot + Kup.transpose() * u_dot;

        // displacement equation
        local_rhs.template segment<displacement_size>(displacement_index)
            .noalias() += Kup * p;
    }

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

private:
    HydroMechanicsProcessData<DisplacementDim>& _process_data;

    using BMatricesType =
        BMatrixPolicyType<ShapeFunctionDisplacement, DisplacementDim>;
    using IpData =
        IntegrationPointData<BMatricesType, ShapeMatricesTypeDisplacement,
                             ShapeMatricesTypePressure, DisplacementDim,
                             ShapeFunctionDisplacement::NPOINTS>;
    std::vector<IpData, Eigen::aligned_allocator<IpData>> _ip_data;

    IntegrationMethod _integration_method;
    MeshLib::Element const& _element;

    static const int pressure_index = 0;
    static const int pressure_size = ShapeFunctionPressure::NPOINTS;
    static const int displacement_index = ShapeFunctionPressure::NPOINTS;
    static const int displacement_size =
        ShapeFunctionDisplacement::NPOINTS * DisplacementDim;
    static const int kelvin_vector_size =
        KelvinVectorDimensions<DisplacementDim>::value;
};

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, unsigned GlobalDim, int DisplacementDim>
class LocalAssemblerData final
    : public HydroMechanicsLocalAssembler<ShapeFunctionDisplacement,
                                          ShapeFunctionPressure,
                                          IntegrationMethod, DisplacementDim>
{
public:
    LocalAssemblerData(LocalAssemblerData const&) = delete;
    LocalAssemblerData(LocalAssemblerData&&) = delete;

    LocalAssemblerData(MeshLib::Element const& e,
                       std::size_t const local_matrix_size,
                       bool is_axially_symmetric,
                       unsigned const integration_order,
                       HydroMechanicsProcessData<DisplacementDim>& process_data)
        : HydroMechanicsLocalAssembler<ShapeFunctionDisplacement,
                                       ShapeFunctionPressure, IntegrationMethod,
                                       DisplacementDim>(
              e, local_matrix_size, is_axially_symmetric, integration_order,
              process_data)
    {
    }
};

}  // namespace HydroMechanics
}  // namespace ProcessLib
