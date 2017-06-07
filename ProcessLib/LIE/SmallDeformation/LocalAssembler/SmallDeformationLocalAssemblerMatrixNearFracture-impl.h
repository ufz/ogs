/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
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
      _process_data(process_data),
      _integration_method(integration_order),
      _element(e),
      _is_axially_symmetric(is_axially_symmetric)
{
    std::vector<
        ShapeMatrices,
        Eigen::aligned_allocator<typename ShapeMatricesType::ShapeMatrices>>
        shape_matrices = initShapeMatrices<ShapeFunction,
                                           ShapeMatricesType,
                                           IntegrationMethod,
                                           DisplacementDim>(
            e, is_axially_symmetric, _integration_method);

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    _ip_data.reserve(n_integration_points);
    _secondary_data.N.resize(n_integration_points);

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        _ip_data.emplace_back(*_process_data._material);
        auto& ip_data = _ip_data[ip];
        auto const& sm = shape_matrices[ip];
        ip_data.N = sm.N;
        ip_data.dNdx = sm.dNdx;
        ip_data._detJ = sm.detJ;
        ip_data._integralMeasure = sm.integralMeasure;

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

    for (auto fid : process_data._vec_ele_connected_fractureIDs[e.getID()])
        _fracture_props.push_back(_process_data._vec_fracture_property[fid].get());
}


template <typename ShapeFunction, typename IntegrationMethod,
          int DisplacementDim>
void SmallDeformationLocalAssemblerMatrixNearFracture<ShapeFunction, IntegrationMethod,
                                    DisplacementDim>::
assembleWithJacobian(
    double const t,
    Eigen::VectorXd const& local_u,
    Eigen::VectorXd& local_b,
    Eigen::MatrixXd& local_J)
{
    assert (_element.getDimension() == DisplacementDim);

    auto const N_DOF_PER_VAR = ShapeFunction::NPOINTS * DisplacementDim;
    auto const n_fractures = _fracture_props.size();

    using BlockVectorType = typename Eigen::VectorXd::FixedSegmentReturnType<N_DOF_PER_VAR>::Type;
    using BlockMatrixType = Eigen::Block<Eigen::MatrixXd,N_DOF_PER_VAR,N_DOF_PER_VAR>;

    //--------------------------------------------------------------------------------------
    // prepare sub vectors, matrices for regular displacement (u) and displacement jumps (g)
    //
    // example with two fractures:
    //     |b(u)|
    // b = |b(g1)|
    //     |b(g2)|
    //
    //     |J(u,u)  J(u,g1)  J(u,g2) |
    // J = |J(g1,u) J(g1,g1) J(g1,g2)|
    //     |J(g2,u) J(g2,g1) J(g2,g2)|
    //--------------------------------------------------------------------------------------
    auto local_b_u = local_b.segment<N_DOF_PER_VAR>(0);
    std::vector<BlockVectorType> vec_local_b_g;
    for (unsigned i=0; i<n_fractures; i++)
        vec_local_b_g.push_back(local_b.segment<N_DOF_PER_VAR>(N_DOF_PER_VAR*(i+1)));

    auto local_J_uu = local_J.block<N_DOF_PER_VAR, N_DOF_PER_VAR>(0, 0);
    std::vector<BlockMatrixType> vec_local_J_ug;
    std::vector<BlockMatrixType> vec_local_J_gu;
    std::vector<std::vector<BlockMatrixType>> vec_local_J_gg(n_fractures);
    for (unsigned i=0; i<n_fractures; i++)
    {
        auto sub_ug = local_J.block<N_DOF_PER_VAR, N_DOF_PER_VAR>(0, N_DOF_PER_VAR*(i+1));
        vec_local_J_ug.push_back(sub_ug);

        auto sub_gu = local_J.block<N_DOF_PER_VAR, N_DOF_PER_VAR>(N_DOF_PER_VAR*(i+1), 0);
        vec_local_J_gu.push_back(sub_gu);

        for (unsigned j=0; j<n_fractures; j++)
        {
            auto sub_gg =
                local_J.block<N_DOF_PER_VAR, N_DOF_PER_VAR>(N_DOF_PER_VAR * (i + 1), N_DOF_PER_VAR * (j + 1));
            vec_local_J_gg[i].push_back(sub_gg);
        }
    }

    auto const nodal_u = local_u.segment<N_DOF_PER_VAR>(0);
    std::vector<BlockVectorType> vec_nodal_g;
    for (unsigned i=0; i<n_fractures; i++)
    {
        auto sub = const_cast<Eigen::VectorXd&>(local_u).segment<N_DOF_PER_VAR>(N_DOF_PER_VAR*(i+1));
        vec_nodal_g.push_back(sub);
    }

    //------------------------------------------------
    // integration
    //------------------------------------------------
    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    SpatialPosition x_position;
    x_position.setElementID(_element.getID());

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);

        auto &ip_data = _ip_data[ip];
        auto const& wp = _integration_method.getWeightedPoint(ip);
        auto const ip_factor = ip_data._detJ * wp.getWeight() * ip_data._integralMeasure;

        auto const& N = ip_data.N;
        auto const& dNdx = ip_data.dNdx;

        // levelset functions
        auto const ip_physical_coords = computePhysicalCoordinates(_element, N);
        std::vector<double> levelsets(n_fractures);
        for (unsigned i = 0; i < n_fractures; i++)
            levelsets[i] = calculateLevelSetFunction(*_fracture_props[i],
                                                     ip_physical_coords.getCoords());

        // nodal displacement = u^hat + sum_i(levelset_i(x) * [u]_i)
        NodalDisplacementVectorType nodal_total_u = nodal_u;
        for (unsigned i=0; i<n_fractures; i++)
            nodal_total_u += levelsets[i] * vec_nodal_g[i];

        auto const x_coord =
            interpolateXCoordinate<ShapeFunction, ShapeMatricesType>(_element,
                                                                     N);
        auto const B =
            LinearBMatrix::computeBMatrix<DisplacementDim,
                                          ShapeFunction::NPOINTS,
                                          typename BMatricesType::BMatrixType>(
                dNdx, N, x_coord, _is_axially_symmetric);

        // strain, stress
        auto const& eps_prev = ip_data._eps_prev;
        auto const& sigma_prev = ip_data._sigma_prev;

        auto& eps = ip_data._eps;
        auto& sigma = ip_data._sigma;
        auto& state = ip_data._material_state_variables;

        eps.noalias() = B * nodal_total_u;

        auto&& solution = _ip_data[ip]._solid_material.integrateStress(
            t, x_position, _process_data.dt, eps_prev, eps, sigma_prev, *state);

        if (!solution)
            OGS_FATAL("Computation of local constitutive relation failed.");

        KelvinMatrixType<DisplacementDim> C;
        std::tie(sigma, state, C) = std::move(*solution);

        // r_u = B^T * Sigma = B^T * C * B * (u+phi*[u])
        // r_[u] = (phi*B)^T * Sigma = (phi*B)^T * C * B * (u+phi*[u])
        local_b_u.noalias() -= B.transpose() * sigma * ip_factor;
        for (unsigned i=0; i<n_fractures; i++)
            vec_local_b_g[i].noalias() -= levelsets[i] * B.transpose() * sigma * ip_factor;

        // J_uu += B^T * C * B
        local_J_uu.noalias() += B.transpose() * C * B * ip_factor;

        for (unsigned i=0; i<n_fractures; i++)
        {
            // J_u[u] += B^T * C * (levelset * B)
            vec_local_J_ug[i].noalias() +=
                B.transpose() * C * (levelsets[i] * B) * ip_factor;

            // J_[u]u += (levelset * B)^T * C * B
            vec_local_J_gu[i].noalias() +=
                (levelsets[i] * B.transpose()) * C * B * ip_factor;

            for (unsigned j=0; j<n_fractures; j++)
            {
                // J_[u][u] += (levelset * B)^T * C * (levelset * B)
                vec_local_J_gg[i][j].noalias() +=
                    (levelsets[i] * B.transpose()) * C * (levelsets[j] * B) *
                    ip_factor;
            }
        }
    }
}


template <typename ShapeFunction, typename IntegrationMethod,
          int DisplacementDim>
void SmallDeformationLocalAssemblerMatrixNearFracture<ShapeFunction, IntegrationMethod,
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
