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

#include "SmallDeformationLocalAssemblerFracture.h"

#include <Eigen/Eigen>

#include "MathLib/LinAlg/Eigen/EigenMapTools.h"

#include "ProcessLib/Utils/InitShapeMatrices.h"

#include "ProcessLib/LIE/Common/LevelSetFunction.h"

namespace ProcessLib
{
namespace LIE
{
namespace SmallDeformation
{
template <typename ShapeFunction, typename IntegrationMethod,
          int DisplacementDim>
SmallDeformationLocalAssemblerFracture<ShapeFunction, IntegrationMethod,
                                       DisplacementDim>::
    SmallDeformationLocalAssemblerFracture(
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
      shape_matrices_(initShapeMatrices<ShapeFunction, ShapeMatricesType,
                                        IntegrationMethod, DisplacementDim>(
          e, is_axially_symmetric, integration_method_)),
      element_(e)
{
    assert(element_.getDimension() == DisplacementDim - 1);

    unsigned const n_integration_points =
        integration_method_.getNumberOfPoints();

    ip_data_.reserve(n_integration_points);
    secondary_data_.N.resize(n_integration_points);

    auto mat_id = (*process_data_.mesh_prop_materialIDs_)[e.getID()];
    auto frac_id = process_data_.map_materialID_to_fractureID_[mat_id];
    fracture_property_ = &process_data_.fracture_properties[frac_id];
    for (auto fid : process_data.vec_ele_connected_fractureIDs_[e.getID()])
    {
        fracID_to_local_.insert({fid, fracture_props_.size()});
        fracture_props_.push_back(&process_data_.fracture_properties[fid]);
    }

    for (auto jid : process_data.vec_ele_connected_junctionIDs_[e.getID()])
    {
        junction_props_.push_back(&process_data_.junction_properties[jid]);
    }

    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(element_.getID());
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);

        ip_data_.emplace_back(*process_data_.fracture_model_);
        auto const& sm = shape_matrices_[ip];
        auto& ip_data = ip_data_[ip];
        ip_data.integration_weight =
            integration_method_.getWeightedPoint(ip).getWeight() *
            sm.integralMeasure * sm.detJ;
        ip_data.h_matrices.setZero(DisplacementDim,
                                   ShapeFunction::NPOINTS * DisplacementDim);

        computeHMatrix<DisplacementDim, ShapeFunction::NPOINTS,
                       typename ShapeMatricesType::NodalRowVectorType,
                       HMatrixType>(sm.N, ip_data.h_matrices);

        // Initialize current time step values
        ip_data.w.setZero(DisplacementDim);
        ip_data.sigma.setZero(DisplacementDim);

        // Previous time step values are not initialized and are set later.
        ip_data.sigma_prev.resize(DisplacementDim);
        ip_data.w_prev.resize(DisplacementDim);

        ip_data.C.resize(DisplacementDim, DisplacementDim);

        ip_data.aperture0 = fracture_property_->aperture0(0, x_position)[0];
        ip_data.aperture_prev = ip_data.aperture0;

        secondary_data_.N[ip] = sm.N;
    }
}

template <typename ShapeFunction, typename IntegrationMethod,
          int DisplacementDim>
void SmallDeformationLocalAssemblerFracture<
    ShapeFunction, IntegrationMethod,
    DisplacementDim>::assembleWithJacobian(double const t, double const /*dt*/,
                                           Eigen::VectorXd const& local_u,
                                           Eigen::VectorXd& local_b,
                                           Eigen::MatrixXd& local_J)
{
    auto const N_DOF_PER_VAR = ShapeFunction::NPOINTS * DisplacementDim;
    auto const n_fractures = fracture_props_.size();
    auto const n_junctions = junction_props_.size();
    auto const n_enrich_var = n_fractures + n_junctions;

    //--------------------------------------------------------------------------------------
    // prepare sub vectors, matrices for regular displacement (u) and
    // displacement jumps (g)
    //
    // example with two fractures with one intersection:
    // b = |b(g1)|
    //     |b(g2)|
    //
    // J = |J(g1,g1) J(g1,g2)|
    //     |J(g2,g1) J(g2,g2)|
    //--------------------------------------------------------------------------------------

    using BlockVectorType =
        typename Eigen::VectorXd::FixedSegmentReturnType<N_DOF_PER_VAR>::Type;
    using BlockMatrixType =
        Eigen::Block<Eigen::MatrixXd, N_DOF_PER_VAR, N_DOF_PER_VAR>;

    std::vector<BlockVectorType> vec_local_b_g;
    for (unsigned i = 0; i < n_enrich_var; i++)
    {
        vec_local_b_g.push_back(
            local_b.segment<N_DOF_PER_VAR>(N_DOF_PER_VAR * i));
    }
    std::vector<std::vector<BlockMatrixType>> vec_local_J_gg(n_enrich_var);
    for (unsigned i = 0; i < n_enrich_var; i++)
    {
        for (unsigned j = 0; j < n_enrich_var; j++)
        {
            auto sub_gg = local_J.block<N_DOF_PER_VAR, N_DOF_PER_VAR>(
                N_DOF_PER_VAR * i, N_DOF_PER_VAR * j);
            vec_local_J_gg[i].push_back(sub_gg);
        }
    }

    std::vector<Eigen::VectorXd> vec_nodal_g;
    for (unsigned i = 0; i < n_enrich_var; i++)
    {
        auto sub = const_cast<Eigen::VectorXd&>(local_u).segment<N_DOF_PER_VAR>(
            N_DOF_PER_VAR * i);
        vec_nodal_g.push_back(sub);
    }

    //------------------------------------------------
    // integration
    //------------------------------------------------
    // the index of a normal (normal to a fracture plane) component
    // in a displacement vector
    int const index_normal = DisplacementDim - 1;
    auto const& R = fracture_property_->R;

    unsigned const n_integration_points =
        integration_method_.getNumberOfPoints();

    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(element_.getID());

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);

        auto& ip_data = ip_data_[ip];
        auto const& integration_weight = ip_data.integration_weight;
        auto const& H = ip_data.h_matrices;
        auto& mat = ip_data.fracture_material;
        auto& sigma = ip_data.sigma;
        auto const& sigma_prev = ip_data.sigma_prev;
        auto& w = ip_data.w;
        auto const& w_prev = ip_data.w_prev;
        auto& C = ip_data.C;
        auto& state = *ip_data.material_state_variables;
        auto& N = secondary_data_.N[ip];

        Eigen::Vector3d const ip_physical_coords(
            computePhysicalCoordinates(element_, N).getCoords());
        std::vector<double> const levelsets(duGlobalEnrichments(
            fracture_property_->fracture_id, fracture_props_, junction_props_,
            fracID_to_local_, ip_physical_coords));

        // du = du^hat + sum_i(enrich^br_i(x) * [u]i_) + sum_i(enrich^junc_i(x)
        // * [u]i_)
        Eigen::VectorXd nodal_gap(N_DOF_PER_VAR);
        nodal_gap.setZero();
        for (unsigned i = 0; i < n_enrich_var; i++)
        {
            nodal_gap += levelsets[i] * vec_nodal_g[i];
        }

        // displacement jumps
        w.noalias() = R * H * nodal_gap;

        // total aperture
        ip_data.aperture = ip_data.aperture0 + w[index_normal];

        // local C, local stress
        mat.computeConstitutiveRelation(
            t, x_position, ip_data.aperture0,
            Eigen::Matrix<double, DisplacementDim, 1>::Zero(),  // TODO (naumov)
                                                                // Replace with
                                                                // initial
                                                                // stress values
            w_prev, w, sigma_prev, sigma, C, state);

        // r_[u] += H^T*Stress
        for (unsigned i = 0; i < n_enrich_var; i++)
        {
            vec_local_b_g[i].noalias() -= levelsets[i] * H.transpose() *
                                          R.transpose() * sigma *
                                          integration_weight;
        }

        // J_[u][u] += H^T*C*H
        for (unsigned i = 0; i < n_enrich_var; i++)
        {
            for (unsigned j = 0; j < n_enrich_var; j++)
            {
                // J_[u][u] += (levelset * B)^T * C * (levelset * B)
                vec_local_J_gg[i][j].noalias() +=
                    (levelsets[i] * H.transpose() * R.transpose()) * C *
                    (levelsets[j] * R * H) * integration_weight;
            }
        }
    }
}

template <typename ShapeFunction, typename IntegrationMethod,
          int DisplacementDim>
void SmallDeformationLocalAssemblerFracture<ShapeFunction, IntegrationMethod,
                                            DisplacementDim>::
    computeSecondaryVariableConcreteWithVector(const double t,
                                               Eigen::VectorXd const& local_u)
{
    auto const N_DOF_PER_VAR = ShapeFunction::NPOINTS * DisplacementDim;
    auto const n_fractures = fracture_props_.size();
    auto const n_junctions = junction_props_.size();
    auto const n_enrich_var = n_fractures + n_junctions;

    auto const& R = fracture_property_->R;

    // the index of a normal (normal to a fracture plane) component
    // in a displacement vector
    int const index_normal = DisplacementDim - 1;

    unsigned const n_integration_points =
        integration_method_.getNumberOfPoints();

    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(element_.getID());

    std::vector<Eigen::VectorXd> vec_nodal_g;
    for (unsigned i = 0; i < n_enrich_var; i++)
    {
        auto sub = const_cast<Eigen::VectorXd&>(local_u).segment<N_DOF_PER_VAR>(
            N_DOF_PER_VAR * i);
        vec_nodal_g.push_back(sub);
    }

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);

        auto& ip_data = ip_data_[ip];
        auto const& H = ip_data.h_matrices;
        auto& mat = ip_data.fracture_material;
        auto& sigma = ip_data.sigma;
        auto const& sigma_prev = ip_data.sigma_prev;
        auto& w = ip_data.w;
        auto const& w_prev = ip_data.w_prev;
        auto& C = ip_data.C;
        auto& state = *ip_data.material_state_variables;
        auto& bm_ = ip_data.aperture;
        auto& N = secondary_data_.N[ip];

        Eigen::Vector3d const ip_physical_coords(
            computePhysicalCoordinates(element_, N).getCoords());
        std::vector<double> const levelsets(duGlobalEnrichments(
            fracture_property_->fracture_id, fracture_props_, junction_props_,
            fracID_to_local_, ip_physical_coords));

        // du = du^hat + sum_i(enrich^br_i(x) * [u]i_) + sum_i(enrich^junc_i(x)
        // * [u]i_)
        Eigen::VectorXd nodal_gap(N_DOF_PER_VAR);
        nodal_gap.setZero();
        for (unsigned i = 0; i < n_enrich_var; i++)
        {
            nodal_gap += levelsets[i] * vec_nodal_g[i];
        }

        // displacement jumps in local coordinates
        w.noalias() = R * H * nodal_gap;

        // aperture
        bm_ = ip_data.aperture0 + w[index_normal];
        if (bm_ < 0.0)
        {
            OGS_FATAL(
                "Element {:d}, gp {:d}: Fracture aperture is {:g}, but it must "
                "be non-negative.",
                element_.getID(), ip, bm_);
        }

        // local C, local stress
        mat.computeConstitutiveRelation(
            t, x_position, ip_data.aperture0,
            Eigen::Matrix<double, DisplacementDim, 1>::Zero(),  // TODO (naumov)
                                                                // Replace with
                                                                // initial
                                                                // stress values
            w_prev, w, sigma_prev, sigma, C, state);
    }

    double ele_b = 0;
    typename HMatricesType::ForceVectorType ele_sigma =
        HMatricesType::ForceVectorType::Zero(DisplacementDim);
    typename HMatricesType::ForceVectorType ele_w =
        HMatricesType::ForceVectorType::Zero(DisplacementDim);

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        auto& ip_data = ip_data_[ip];
        ele_b += ip_data.aperture;
        ele_w += ip_data.w;
        ele_sigma += ip_data.sigma;
    }
    ele_b /= n_integration_points;
    ele_w /= n_integration_points;
    ele_sigma /= n_integration_points;
    (*process_data_.mesh_prop_b_)[element_.getID()] = ele_b;
    (*process_data_.mesh_prop_w_n_)[element_.getID()] = ele_w[index_normal];
    (*process_data_.mesh_prop_w_s_)[element_.getID()] = ele_w[0];
    (*process_data_.mesh_prop_fracture_stress_normal_)[element_.getID()] =
        ele_sigma[index_normal];
    (*process_data_.mesh_prop_fracture_stress_shear_)[element_.getID()] =
        ele_sigma[0];
}

}  // namespace SmallDeformation
}  // namespace LIE
}  // namespace ProcessLib
