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

#include "SmallDeformationLocalAssemblerMatrixNearFracture.h"

#include <valarray>
#include <vector>

#include <Eigen/Eigen>

#include "MaterialLib/PhysicalConstant.h"

#include "MathLib/KelvinVector.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "MathLib/Point3d.h"

#include "MeshLib/Node.h"

#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"

#include "ProcessLib/Deformation/BMatrixPolicy.h"
#include "ProcessLib/Deformation/LinearBMatrix.h"
#include "ProcessLib/Utils/InitShapeMatrices.h"

#include "ProcessLib/LIE/Common/LevelSetFunction.h"
#include "ProcessLib/LIE/Common/Utils.h"

#include "IntegrationPointDataMatrix.h"
#include "SecondaryData.h"
#include "SmallDeformationLocalAssemblerInterface.h"

namespace ProcessLib
{
namespace LIE
{
namespace SmallDeformation
{
template <typename ShapeFunction,
          typename IntegrationMethod,
          int DisplacementDim>
SmallDeformationLocalAssemblerMatrixNearFracture<ShapeFunction,
                                                 IntegrationMethod,
                                                 DisplacementDim>::
    SmallDeformationLocalAssemblerMatrixNearFracture(
        MeshLib::Element const& e,
        std::size_t const n_variables,
        std::size_t const /*local_matrix_size*/,
        std::vector<unsigned> const& dofIndex_to_localIndex,
        bool const is_axially_symmetric,
        unsigned const integration_order,
        SmallDeformationProcessData<DisplacementDim>& process_data)
    : SmallDeformationLocalAssemblerInterface(
          n_variables * ShapeFunction::NPOINTS * DisplacementDim,
          dofIndex_to_localIndex),
      process_data_(process_data),
      integration_method_(integration_order),
      element_(e),
      is_axially_symmetric_(is_axially_symmetric)
{
    std::vector<
        ShapeMatrices,
        Eigen::aligned_allocator<typename ShapeMatricesType::ShapeMatrices>>
        shape_matrices = initShapeMatrices<ShapeFunction,
                                           ShapeMatricesType,
                                           IntegrationMethod,
                                           DisplacementDim>(
            e, is_axially_symmetric, integration_method_);

    unsigned const n_integration_points =
        integration_method_.getNumberOfPoints();

    ip_data_.reserve(n_integration_points);
    secondary_data_.N.resize(n_integration_points);

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

    for (auto fid : process_data.vec_ele_connected_fractureIDs_[e.getID()])
    {
        fracID_to_local_.insert({fid, fracture_props_.size()});
        fracture_props_.push_back(&process_data_.fracture_properties[fid]);
    }

    for (auto jid : process_data.vec_ele_connected_junctionIDs_[e.getID()])
    {
        junction_props_.push_back(&process_data_.junction_properties[jid]);
    }
}

template <typename ShapeFunction, typename IntegrationMethod,
          int DisplacementDim>
void SmallDeformationLocalAssemblerMatrixNearFracture<
    ShapeFunction, IntegrationMethod,
    DisplacementDim>::assembleWithJacobian(double const t, double const dt,
                                           Eigen::VectorXd const& local_u,
                                           Eigen::VectorXd& local_b,
                                           Eigen::MatrixXd& local_J)
{
    assert(element_.getDimension() == DisplacementDim);

    auto const N_DOF_PER_VAR = ShapeFunction::NPOINTS * DisplacementDim;
    auto const n_fractures = fracture_props_.size();
    auto const n_junctions = junction_props_.size();
    auto const n_enrich_var = n_fractures + n_junctions;

    using BlockVectorType =
        typename Eigen::VectorXd::FixedSegmentReturnType<N_DOF_PER_VAR>::Type;
    using BlockMatrixType =
        Eigen::Block<Eigen::MatrixXd, N_DOF_PER_VAR, N_DOF_PER_VAR>;

    //--------------------------------------------------------------------------------------
    // prepare sub vectors, matrices for regular displacement (u) and
    // displacement jumps (g)
    //
    // example with two fractures with one intersection:
    //     |b(u)|
    // b = |b(g1)|
    //     |b(g2)|
    //     |b(j1)|
    //
    //     |J(u,u)  J(u,g1)  J(u,g2)  J(u,j1) |
    // J = |J(g1,u) J(g1,g1) J(g1,g2) J(g1,j1)|
    //     |J(g2,u) J(g2,g1) J(g2,g2) J(g2,j1)|
    //     |J(j1,u) J(j1,g1) J(j1,g2) J(j1,j1)|
    //--------------------------------------------------------------------------------------
    auto local_b_u = local_b.segment<N_DOF_PER_VAR>(0);
    std::vector<BlockVectorType> vec_local_b_g;
    for (unsigned i = 0; i < n_enrich_var; i++)
    {
        vec_local_b_g.push_back(
            local_b.segment<N_DOF_PER_VAR>(N_DOF_PER_VAR * (i + 1)));
    }

    auto local_J_uu = local_J.block<N_DOF_PER_VAR, N_DOF_PER_VAR>(0, 0);
    std::vector<BlockMatrixType> vec_local_J_ug;
    std::vector<BlockMatrixType> vec_local_J_gu;
    std::vector<std::vector<BlockMatrixType>> vec_local_J_gg(n_enrich_var);
    for (unsigned i = 0; i < n_enrich_var; i++)
    {
        auto sub_ug = local_J.block<N_DOF_PER_VAR, N_DOF_PER_VAR>(
            0, N_DOF_PER_VAR * (i + 1));
        vec_local_J_ug.push_back(sub_ug);

        auto sub_gu = local_J.block<N_DOF_PER_VAR, N_DOF_PER_VAR>(
            N_DOF_PER_VAR * (i + 1), 0);
        vec_local_J_gu.push_back(sub_gu);

        for (unsigned j = 0; j < n_enrich_var; j++)
        {
            auto sub_gg = local_J.block<N_DOF_PER_VAR, N_DOF_PER_VAR>(
                N_DOF_PER_VAR * (i + 1), N_DOF_PER_VAR * (j + 1));
            vec_local_J_gg[i].push_back(sub_gg);
        }
    }

    auto const nodal_u = local_u.segment<N_DOF_PER_VAR>(0);
    std::vector<BlockVectorType> vec_nodal_g;
    for (unsigned i = 0; i < n_enrich_var; i++)
    {
        auto sub = const_cast<Eigen::VectorXd&>(local_u).segment<N_DOF_PER_VAR>(
            N_DOF_PER_VAR * (i + 1));
        vec_nodal_g.push_back(sub);
    }

    //------------------------------------------------
    // integration
    //------------------------------------------------
    unsigned const n_integration_points =
        integration_method_.getNumberOfPoints();

    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(element_.getID());

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);

        auto& ip_data = ip_data_[ip];
        auto const& w = ip_data_[ip].integration_weight;

        auto const& N = ip_data.N;
        auto const& dNdx = ip_data.dNdx;

        // levelset functions
        Eigen::Vector3d const ip_physical_coords(
            computePhysicalCoordinates(element_, N).getCoords());
        std::vector<double> const levelsets(
            uGlobalEnrichments(fracture_props_, junction_props_,
                               fracID_to_local_, ip_physical_coords));

        // u = u^hat + sum_i(enrich^br_i(x) * [u]i_) + sum_i(enrich^junc_i(x) *
        // [u]i_)
        NodalDisplacementVectorType nodal_total_u = nodal_u;
        for (unsigned i = 0; i < n_enrich_var; i++)
        {
            nodal_total_u += levelsets[i] * vec_nodal_g[i];
        }

        auto const x_coord =
            interpolateXCoordinate<ShapeFunction, ShapeMatricesType>(element_,
                                                                     N);
        auto const B =
            LinearBMatrix::computeBMatrix<DisplacementDim,
                                          ShapeFunction::NPOINTS,
                                          typename BMatricesType::BMatrixType>(
                dNdx, N, x_coord, is_axially_symmetric_);

        // strain, stress
        auto const& eps_prev = ip_data.eps_prev_;
        auto const& sigma_prev = ip_data.sigma_prev_;

        auto& eps = ip_data.eps_;
        auto& sigma = ip_data.sigma_;
        auto& state = ip_data.material_state_variables_;

        eps.noalias() = B * nodal_total_u;

        auto&& solution = ip_data_[ip].solid_material_.integrateStress(
            t, x_position, dt, eps_prev, eps, sigma_prev, *state,
            process_data_.reference_temperature_);

        if (!solution)
        {
            OGS_FATAL("Computation of local constitutive relation failed.");
        }

        MathLib::KelvinVector::KelvinMatrixType<DisplacementDim> C;
        std::tie(sigma, state, C) = std::move(*solution);

        // r_u = B^T * Sigma = B^T * C * B * (u+phi*[u])
        // r_[u] = (phi*B)^T * Sigma = (phi*B)^T * C * B * (u+phi*[u])
        local_b_u.noalias() -= B.transpose() * sigma * w;
        for (unsigned i = 0; i < n_enrich_var; i++)
        {
            vec_local_b_g[i].noalias() -=
                levelsets[i] * B.transpose() * sigma * w;
        }

        // J_uu += B^T * C * B
        local_J_uu.noalias() += B.transpose() * C * B * w;

        for (unsigned i = 0; i < n_enrich_var; i++)
        {
            // J_u[u] += B^T * C * (levelset * B)
            vec_local_J_ug[i].noalias() +=
                B.transpose() * C * (levelsets[i] * B) * w;

            // J_[u]u += (levelset * B)^T * C * B
            vec_local_J_gu[i].noalias() +=
                (levelsets[i] * B.transpose()) * C * B * w;

            for (unsigned j = 0; j < n_enrich_var; j++)
            {
                // J_[u][u] += (levelset * B)^T * C * (levelset * B)
                vec_local_J_gg[i][j].noalias() +=
                    (levelsets[i] * B.transpose()) * C * (levelsets[j] * B) * w;
            }
        }
    }
}

template <typename ShapeFunction,
          typename IntegrationMethod,
          int DisplacementDim>
void SmallDeformationLocalAssemblerMatrixNearFracture<ShapeFunction,
                                                      IntegrationMethod,
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
