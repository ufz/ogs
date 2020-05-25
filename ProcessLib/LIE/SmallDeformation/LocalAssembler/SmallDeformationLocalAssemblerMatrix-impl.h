/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "SmallDeformationLocalAssemblerMatrix.h"

#include <valarray>
#include <vector>

#include <Eigen/Eigen>

#include "MaterialLib/PhysicalConstant.h"
#include "MaterialLib/SolidModels/SelectSolidConstitutiveRelation.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"

#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "ProcessLib/Deformation/BMatrixPolicy.h"
#include "ProcessLib/Deformation/LinearBMatrix.h"
#include "ProcessLib/Utils/InitShapeMatrices.h"

#include "ProcessLib/LIE/SmallDeformation/SmallDeformationProcessData.h"

#include "IntegrationPointDataMatrix.h"
#include "SecondaryData.h"
#include "SmallDeformationLocalAssemblerInterface.h"

namespace ProcessLib
{
namespace LIE
{
namespace SmallDeformation
{
template <typename ShapeFunction, typename IntegrationMethod,
          int DisplacementDim>
SmallDeformationLocalAssemblerMatrix<ShapeFunction, IntegrationMethod,
                                     DisplacementDim>::
    SmallDeformationLocalAssemblerMatrix(
        MeshLib::Element const& e,
        std::size_t const /*local_matrix_size*/,
        bool const is_axially_symmetric,
        unsigned const integration_order,
        SmallDeformationProcessData<DisplacementDim>& process_data)
    : process_data_(process_data),
      integration_method_(integration_order),
      element_(e),
      is_axially_symmetric_(is_axially_symmetric)
{
    unsigned const n_integration_points =
        integration_method_.getNumberOfPoints();

    ip_data_.reserve(n_integration_points);
    secondary_data_.N.resize(n_integration_points);

    auto const shape_matrices =
        initShapeMatrices<ShapeFunction, ShapeMatricesType, IntegrationMethod,
                          DisplacementDim>(e, is_axially_symmetric,
                                           integration_method_);

    auto& solid_material = MaterialLib::Solids::selectSolidConstitutiveRelation(
        process_data_.solid_materials, process_data_.material_ids, e.getID());

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        ip_data_.emplace_back(solid_material);
        auto& ip_data = ip_data_[ip];
        auto const& sm = shape_matrices[ip];
        ip_data.N = sm.N;
        ip_data.dNdx = sm.dNdx;
        ip_data.integration_weight =
            integration_method_.getWeightedPoint(ip).getWeight() *
            sm.integralMeasure * sm.detJ;

        // Initialize current time step values
        static const int kelvin_vector_size =
            MathLib::KelvinVector::KelvinVectorDimensions<
                DisplacementDim>::value;
        ip_data.sigma_.setZero(kelvin_vector_size);
        ip_data.eps_.setZero(kelvin_vector_size);

        // Previous time step values are not initialized and are set later.
        ip_data.sigma_prev_.resize(kelvin_vector_size);
        ip_data.eps_prev_.resize(kelvin_vector_size);

        ip_data.C_.resize(kelvin_vector_size, kelvin_vector_size);

        secondary_data_.N[ip] = sm.N;
    }
}

template <typename ShapeFunction, typename IntegrationMethod,
          int DisplacementDim>
void SmallDeformationLocalAssemblerMatrix<ShapeFunction, IntegrationMethod,
                                          DisplacementDim>::
    assembleWithJacobian(double const t, double const dt,
                         std::vector<double> const& local_x,
                         std::vector<double> const& /*local_xdot*/,
                         const double /*dxdot_dx*/, const double /*dx_dx*/,
                         std::vector<double>& /*local_M_data*/,
                         std::vector<double>& /*local_K_data*/,
                         std::vector<double>& local_b_data,
                         std::vector<double>& local_Jac_data)
{
    assert(element_.getDimension() == DisplacementDim);

    auto const local_matrix_size = local_x.size();

    auto local_Jac = MathLib::createZeroedMatrix<StiffnessMatrixType>(
        local_Jac_data, local_matrix_size, local_matrix_size);

    auto local_b = MathLib::createZeroedVector<NodalDisplacementVectorType>(
        local_b_data, local_matrix_size);

    unsigned const n_integration_points =
        integration_method_.getNumberOfPoints();

    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(element_.getID());

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);
        auto const& w = ip_data_[ip].integration_weight;

        auto const& N = ip_data_[ip].N;
        auto const& dNdx = ip_data_[ip].dNdx;
        auto const x_coord =
            interpolateXCoordinate<ShapeFunction, ShapeMatricesType>(element_,
                                                                     N);
        auto const B =
            LinearBMatrix::computeBMatrix<DisplacementDim,
                                          ShapeFunction::NPOINTS,
                                          typename BMatricesType::BMatrixType>(
                dNdx, N, x_coord, is_axially_symmetric_);

        auto const& eps_prev = ip_data_[ip].eps_prev_;
        auto const& sigma_prev = ip_data_[ip].sigma_prev_;

        auto& eps = ip_data_[ip].eps_;
        auto& sigma = ip_data_[ip].sigma_;
        auto& state = ip_data_[ip].material_state_variables_;

        eps.noalias() =
            B * Eigen::Map<typename BMatricesType::NodalForceVectorType const>(
                    local_x.data(), ShapeFunction::NPOINTS * DisplacementDim);

        auto&& solution = ip_data_[ip].solid_material_.integrateStress(
            t, x_position, dt, eps_prev, eps, sigma_prev, *state,
            process_data_.reference_temperature_);

        if (!solution)
        {
            OGS_FATAL("Computation of local constitutive relation failed.");
        }

        MathLib::KelvinVector::KelvinMatrixType<DisplacementDim> C;
        std::tie(sigma, state, C) = std::move(*solution);

        local_b.noalias() -= B.transpose() * sigma * w;
        local_Jac.noalias() += B.transpose() * C * B * w;
    }
}

template <typename ShapeFunction, typename IntegrationMethod,
          int DisplacementDim>
void SmallDeformationLocalAssemblerMatrix<ShapeFunction, IntegrationMethod,
                                          DisplacementDim>::
    computeSecondaryVariableConcreteWithVector(
        double const /*t*/, Eigen::VectorXd const& /*local_x*/)
{
    // Compute average value per element
    const int n = DisplacementDim == 2 ? 4 : 6;
    Eigen::VectorXd ele_stress = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd ele_strain = Eigen::VectorXd::Zero(n);

    unsigned const n_integration_points =
        integration_method_.getNumberOfPoints();
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        auto& ip_data = ip_data_[ip];

        ele_stress += ip_data.sigma_;
        ele_strain += ip_data.eps_;
    }
    ele_stress /= n_integration_points;
    ele_strain /= n_integration_points;

    (*process_data_.mesh_prop_stress_xx_)[element_.getID()] = ele_stress[0];
    (*process_data_.mesh_prop_stress_yy_)[element_.getID()] = ele_stress[1];
    (*process_data_.mesh_prop_stress_zz_)[element_.getID()] = ele_stress[2];
    (*process_data_.mesh_prop_stress_xy_)[element_.getID()] = ele_stress[3];
    if (DisplacementDim == 3)
    {
        (*process_data_.mesh_prop_stress_yz_)[element_.getID()] = ele_stress[4];
        (*process_data_.mesh_prop_stress_xz_)[element_.getID()] = ele_stress[5];
    }

    (*process_data_.mesh_prop_strain_xx_)[element_.getID()] = ele_strain[0];
    (*process_data_.mesh_prop_strain_yy_)[element_.getID()] = ele_strain[1];
    (*process_data_.mesh_prop_strain_zz_)[element_.getID()] = ele_strain[2];
    (*process_data_.mesh_prop_strain_xy_)[element_.getID()] = ele_strain[3];
    if (DisplacementDim == 3)
    {
        (*process_data_.mesh_prop_strain_yz_)[element_.getID()] = ele_strain[4];
        (*process_data_.mesh_prop_strain_xz_)[element_.getID()] = ele_strain[5];
    }
}

}  // namespace SmallDeformation
}  // namespace LIE
}  // namespace ProcessLib
