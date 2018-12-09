/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
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
      _process_data(process_data),
      _integration_method(integration_order),
      _shape_matrices(initShapeMatrices<ShapeFunction, ShapeMatricesType,
                                        IntegrationMethod, DisplacementDim>(
          e, is_axially_symmetric, _integration_method)),
      _element(e)
{
    assert(_element.getDimension() == DisplacementDim - 1);

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    _ip_data.reserve(n_integration_points);
    _secondary_data.N.resize(n_integration_points);

    auto mat_id = (*_process_data._mesh_prop_materialIDs)[e.getID()];
    auto frac_id = _process_data._map_materialID_to_fractureID[mat_id];
    _fracture_property = &_process_data._vec_fracture_property[frac_id];
    for (auto fid : process_data._vec_ele_connected_fractureIDs[e.getID()])
    {
        _fracID_to_local.insert({fid, _fracture_props.size()});
        _fracture_props.push_back(&_process_data._vec_fracture_property[fid]);
    }

    for (auto jid : process_data._vec_ele_connected_junctionIDs[e.getID()])
    {
        _junction_props.push_back(&_process_data._vec_junction_property[jid]);
    }

    SpatialPosition x_position;
    x_position.setElementID(_element.getID());
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);

        _ip_data.emplace_back(*_process_data._fracture_model);
        auto const& sm = _shape_matrices[ip];
        auto& ip_data = _ip_data[ip];
        ip_data.integration_weight =
            _integration_method.getWeightedPoint(ip).getWeight() *
            sm.integralMeasure * sm.detJ;
        ip_data._h_matrices.setZero(DisplacementDim,
                                    ShapeFunction::NPOINTS * DisplacementDim);

        computeHMatrix<DisplacementDim, ShapeFunction::NPOINTS,
                       typename ShapeMatricesType::NodalRowVectorType,
                       HMatrixType>(sm.N, ip_data._h_matrices);

        // Initialize current time step values
        ip_data._w.setZero(DisplacementDim);
        ip_data._sigma.setZero(DisplacementDim);

        // Previous time step values are not initialized and are set later.
        ip_data._sigma_prev.resize(DisplacementDim);
        ip_data._w_prev.resize(DisplacementDim);

        ip_data._C.resize(DisplacementDim, DisplacementDim);

        ip_data._aperture0 = _fracture_property->aperture0(0, x_position)[0];
        ip_data._aperture_prev = ip_data._aperture0;

        _secondary_data.N[ip] = sm.N;
    }
}

template <typename ShapeFunction, typename IntegrationMethod,
          int DisplacementDim>
void SmallDeformationLocalAssemblerFracture<
    ShapeFunction, IntegrationMethod,
    DisplacementDim>::assembleWithJacobian(double const t,
                                           Eigen::VectorXd const& local_u,
                                           Eigen::VectorXd& local_b,
                                           Eigen::MatrixXd& local_J)
{
    auto const N_DOF_PER_VAR = ShapeFunction::NPOINTS * DisplacementDim;
    auto const n_fractures = _fracture_props.size();
    auto const n_junctions = _junction_props.size();
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
    auto const& R = _fracture_property->R;

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    SpatialPosition x_position;
    x_position.setElementID(_element.getID());

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);

        auto& ip_data = _ip_data[ip];
        auto const& integration_weight = ip_data.integration_weight;
        auto const& H = ip_data._h_matrices;
        auto& mat = ip_data._fracture_material;
        auto& sigma = ip_data._sigma;
        auto const& sigma_prev = ip_data._sigma_prev;
        auto& w = ip_data._w;
        auto const& w_prev = ip_data._w_prev;
        auto& C = ip_data._C;
        auto& state = *ip_data._material_state_variables;
        auto& N = _secondary_data.N[ip];

        Eigen::Vector3d const ip_physical_coords(
            computePhysicalCoordinates(_element, N).getCoords());
        std::vector<double> const levelsets(duGlobalEnrichments(
            _fracture_property->fracture_id, _fracture_props, _junction_props,
            _fracID_to_local, ip_physical_coords));

        // du = du^hat + sum_i(enrich^br_i(x) * [u]_i) + sum_i(enrich^junc_i(x)
        // * [u]_i)
        Eigen::VectorXd nodal_gap(N_DOF_PER_VAR);
        nodal_gap.setZero();
        for (unsigned i = 0; i < n_enrich_var; i++)
        {
            nodal_gap += levelsets[i] * vec_nodal_g[i];
        }

        // displacement jumps
        w.noalias() = R * H * nodal_gap;

        // total aperture
        ip_data._aperture = ip_data._aperture0 + w[index_normal];

        // local C, local stress
        mat.computeConstitutiveRelation(
            t, x_position, ip_data._aperture0,
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
    auto const n_fractures = _fracture_props.size();
    auto const n_junctions = _junction_props.size();
    auto const n_enrich_var = n_fractures + n_junctions;

    auto const& R = _fracture_property->R;

    // the index of a normal (normal to a fracture plane) component
    // in a displacement vector
    int const index_normal = DisplacementDim - 1;

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    SpatialPosition x_position;
    x_position.setElementID(_element.getID());

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

        auto& ip_data = _ip_data[ip];
        auto const& H = ip_data._h_matrices;
        auto& mat = ip_data._fracture_material;
        auto& sigma = ip_data._sigma;
        auto const& sigma_prev = ip_data._sigma_prev;
        auto& w = ip_data._w;
        auto const& w_prev = ip_data._w_prev;
        auto& C = ip_data._C;
        auto& state = *ip_data._material_state_variables;
        auto& b_m = ip_data._aperture;
        auto& N = _secondary_data.N[ip];

        Eigen::Vector3d const ip_physical_coords(
            computePhysicalCoordinates(_element, N).getCoords());
        std::vector<double> const levelsets(duGlobalEnrichments(
            _fracture_property->fracture_id, _fracture_props, _junction_props,
            _fracID_to_local, ip_physical_coords));

        // du = du^hat + sum_i(enrich^br_i(x) * [u]_i) + sum_i(enrich^junc_i(x)
        // * [u]_i)
        Eigen::VectorXd nodal_gap(N_DOF_PER_VAR);
        nodal_gap.setZero();
        for (unsigned i = 0; i < n_enrich_var; i++)
        {
            nodal_gap += levelsets[i] * vec_nodal_g[i];
        }

        // displacement jumps in local coordinates
        w.noalias() = R * H * nodal_gap;

        // aperture
        b_m = ip_data._aperture0 + w[index_normal];
        if (b_m < 0.0)
            OGS_FATAL(
                "Element %d, gp %d: Fracture aperture is %g, but it must be "
                "non-negative.",
                _element.getID(), ip, b_m);

        // local C, local stress
        mat.computeConstitutiveRelation(
            t, x_position, ip_data._aperture0,
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
        ele_b += _ip_data[ip]._aperture;
        ele_w += _ip_data[ip]._w;
        ele_sigma += _ip_data[ip]._sigma;
    }
    ele_b /= n_integration_points;
    ele_w /= n_integration_points;
    ele_sigma /= n_integration_points;
    (*_process_data._mesh_prop_b)[_element.getID()] = ele_b;
    (*_process_data._mesh_prop_w_n)[_element.getID()] = ele_w[index_normal];
    (*_process_data._mesh_prop_w_s)[_element.getID()] = ele_w[0];
    (*_process_data._mesh_prop_fracture_stress_normal)[_element.getID()] =
        ele_sigma[index_normal];
    (*_process_data._mesh_prop_fracture_stress_shear)[_element.getID()] =
        ele_sigma[0];
}

}  // namespace SmallDeformation
}  // namespace LIE
}  // namespace ProcessLib
