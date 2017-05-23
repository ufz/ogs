/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
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
        bool is_axially_symmetric,
        unsigned const integration_order,
        SmallDeformationProcessData<DisplacementDim>& process_data)
    : _process_data(process_data),
      _integration_method(integration_order),
      _element(e)
{
    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    _ip_data.reserve(n_integration_points);
    _secondary_data.N.resize(n_integration_points);

    auto const shape_matrices =
        initShapeMatrices<ShapeFunction, ShapeMatricesType,
                          IntegrationMethod, DisplacementDim>(
                e, is_axially_symmetric, _integration_method);

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        _ip_data.emplace_back(*_process_data._material);
        auto& ip_data = _ip_data[ip];
        auto const& sm = shape_matrices[ip];
        ip_data._detJ = sm.detJ;
        ip_data._integralMeasure = sm.integralMeasure;
        ip_data._b_matrices.resize(
            KelvinVectorDimensions<DisplacementDim>::value,
            ShapeFunction::NPOINTS * DisplacementDim);

        auto const x_coord =
            interpolateXCoordinate<ShapeFunction, ShapeMatricesType>(e,
                                                                     sm.N);
        LinearBMatrix::computeBMatrix<DisplacementDim,
                                      ShapeFunction::NPOINTS>(
            sm.dNdx, ip_data._b_matrices, is_axially_symmetric, sm.N,
            x_coord);

        ip_data._sigma.resize(KelvinVectorDimensions<DisplacementDim>::value);
        ip_data._sigma_prev.resize(
            KelvinVectorDimensions<DisplacementDim>::value);
        ip_data._eps.resize(KelvinVectorDimensions<DisplacementDim>::value);
        ip_data._eps_prev.resize(
            KelvinVectorDimensions<DisplacementDim>::value);
        ip_data._C.resize(KelvinVectorDimensions<DisplacementDim>::value,
                          KelvinVectorDimensions<DisplacementDim>::value);

        _secondary_data.N[ip] = sm.N;
    }
}


template <typename ShapeFunction, typename IntegrationMethod,
          int DisplacementDim>
void SmallDeformationLocalAssemblerMatrix<ShapeFunction, IntegrationMethod,
                                    DisplacementDim>::
assembleWithJacobian(
    double const t,
    std::vector<double> const& local_x,
    std::vector<double> const& /*local_xdot*/,
    const double /*dxdot_dx*/, const double /*dx_dx*/,
    std::vector<double>& /*local_M_data*/,
    std::vector<double>& /*local_K_data*/,
    std::vector<double>& local_b_data,
    std::vector<double>& local_Jac_data)
{
    assert (_element.getDimension() == DisplacementDim);

    auto const local_matrix_size = local_x.size();

    auto local_Jac = MathLib::createZeroedMatrix<StiffnessMatrixType>(
        local_Jac_data, local_matrix_size, local_matrix_size);

    auto local_b = MathLib::createZeroedVector<NodalDisplacementVectorType>(
        local_b_data, local_matrix_size);

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    SpatialPosition x_position;
    x_position.setElementID(_element.getID());

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);
        auto const& wp = _integration_method.getWeightedPoint(ip);
        auto const& detJ = _ip_data[ip]._detJ;
        auto const& integralMeasure = _ip_data[ip]._integralMeasure;

        auto const& B = _ip_data[ip]._b_matrices;
        auto const& eps_prev = _ip_data[ip]._eps_prev;
        auto const& sigma_prev = _ip_data[ip]._sigma_prev;

        auto& eps = _ip_data[ip]._eps;
        auto& sigma = _ip_data[ip]._sigma;
        auto& state = _ip_data[ip]._material_state_variables;

        eps.noalias() =
            B *
            Eigen::Map<typename BMatricesType::NodalForceVectorType const>(
                local_x.data(), ShapeFunction::NPOINTS * DisplacementDim);

        auto&& solution = _ip_data[ip]._solid_material.integrateStress(
            t, x_position, _process_data.dt, eps_prev, eps, sigma_prev, *state);

        if (!solution)
            OGS_FATAL("Computation of local constitutive relation failed.");

        KelvinMatrixType<DisplacementDim> C;
        std::tie(sigma, state, C) = std::move(*solution);

        local_b.noalias() -=
            B.transpose() * sigma * detJ * wp.getWeight() * integralMeasure;
        local_Jac.noalias() +=
            B.transpose() * C * B * detJ * wp.getWeight() * integralMeasure;
    }
}


template <typename ShapeFunction, typename IntegrationMethod,
          int DisplacementDim>
void SmallDeformationLocalAssemblerMatrix<ShapeFunction, IntegrationMethod,
                                    DisplacementDim>::
postTimestepConcrete(std::vector<double> const& /*local_x*/)
{
    const int n = 3;
    std::valarray<double> ele_stress(0.0, n);
    std::valarray<double> ele_strain(0.0, n);
    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        auto& ip_data = _ip_data[ip];

        ele_stress[0] += ip_data._sigma[0];
        ele_stress[1] += ip_data._sigma[1];
        ele_stress[2] += ip_data._sigma[3];

        ele_strain[0] += ip_data._eps[0];
        ele_strain[1] += ip_data._eps[1];
        ele_strain[2] += ip_data._eps[3];
    }
    ele_stress /= n_integration_points;
    ele_strain /= n_integration_points;
    (*_process_data._mesh_prop_stress_xx)[_element.getID()] = ele_stress[0];
    (*_process_data._mesh_prop_stress_yy)[_element.getID()] = ele_stress[1];
    (*_process_data._mesh_prop_stress_xy)[_element.getID()] = ele_stress[2];
    (*_process_data._mesh_prop_strain_xx)[_element.getID()] = ele_strain[0];
    (*_process_data._mesh_prop_strain_yy)[_element.getID()] = ele_strain[1];
    (*_process_data._mesh_prop_strain_xy)[_element.getID()] = ele_strain[2];
}

}  // namespace SmallDeformation
}  // namespace LIE
}  // namespace ProcessLib
