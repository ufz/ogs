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
        std::size_t const /*local_matrix_size*/,
        std::vector<unsigned> const& dofIndex_to_localIndex,
        bool const is_axially_symmetric,
        unsigned const integration_order,
        SmallDeformationProcessData<DisplacementDim>& process_data)
    : SmallDeformationLocalAssemblerInterface(
          ShapeFunction::NPOINTS * DisplacementDim,  // no intersection
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
    _fracture_property = _process_data._vec_fracture_property[frac_id].get();

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

        ip_data._aperture0 = (*_fracture_property->aperture0)(0, x_position)[0];
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
    auto const& nodal_jump = local_u;

    auto const& R = _fracture_property->R;

    // the index of a normal (normal to a fracture plane) component
    // in a displacement vector
    int const index_normal = DisplacementDim - 1;

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

        // displacement jumps
        w.noalias() = R * H * nodal_jump;

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
        local_b.noalias() -=
            H.transpose() * R.transpose() * sigma * integration_weight;

        // J_[u][u] += H^T*C*H
        local_J.noalias() +=
            H.transpose() * R.transpose() * C * R * H * integration_weight;
    }
}

template <typename ShapeFunction, typename IntegrationMethod,
          int DisplacementDim>
void SmallDeformationLocalAssemblerFracture<ShapeFunction, IntegrationMethod,
                                            DisplacementDim>::
    computeSecondaryVariableConcreteWithVector(const double t,
                                               Eigen::VectorXd const& local_x)
{
    auto const& nodal_jump = local_x;

    auto const& R = _fracture_property->R;

    // the index of a normal (normal to a fracture plane) component
    // in a displacement vector
    int const index_normal = DisplacementDim - 1;

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    SpatialPosition x_position;
    x_position.setElementID(_element.getID());

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

        // displacement jumps in local coordinates
        w.noalias() = R * H * nodal_jump;

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
