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

#include <Eigen/Core>
#include <range/v3/range/conversion.hpp>
#include <range/v3/view/transform.hpp>
#include <vector>

#include "IntegrationPointDataMatrix.h"
#include "MaterialLib/PhysicalConstant.h"
#include "MathLib/KelvinVector.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "MathLib/Point3d.h"
#include "MeshLib/Node.h"
#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "NumLib/Fem/InitShapeMatrices.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "ProcessLib/Deformation/BMatrixPolicy.h"
#include "ProcessLib/Deformation/LinearBMatrix.h"
#include "ProcessLib/LIE/Common/LevelSetFunction.h"
#include "ProcessLib/LIE/Common/Utils.h"
#include "ProcessLib/Utils/SetOrGetIntegrationPointData.h"
#include "SecondaryData.h"
#include "SmallDeformationLocalAssemblerInterface.h"
#include "SmallDeformationLocalAssemblerMatrixNearFracture.h"

namespace ProcessLib
{
namespace LIE
{
namespace SmallDeformation
{
template <typename ShapeFunction,

          int DisplacementDim>
SmallDeformationLocalAssemblerMatrixNearFracture<ShapeFunction,

                                                 DisplacementDim>::
    SmallDeformationLocalAssemblerMatrixNearFracture(
        MeshLib::Element const& e,
        std::size_t const n_variables,
        std::size_t const /*local_matrix_size*/,
        std::vector<unsigned> const& dofIndex_to_localIndex,
        NumLib::GenericIntegrationMethod const& integration_method,
        bool const is_axially_symmetric,
        SmallDeformationProcessData<DisplacementDim>& process_data)
    : SmallDeformationLocalAssemblerInterface(
          n_variables * ShapeFunction::NPOINTS * DisplacementDim,
          dofIndex_to_localIndex),
      _process_data(process_data),
      _integration_method(integration_method),
      _element(e),
      _is_axially_symmetric(is_axially_symmetric)
{
    std::vector<ShapeMatrices, Eigen::aligned_allocator<
                                   typename ShapeMatricesType::ShapeMatrices>>
        shape_matrices =
            NumLib::initShapeMatrices<ShapeFunction, ShapeMatricesType,
                                      DisplacementDim>(e, is_axially_symmetric,
                                                       _integration_method);

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    _ip_data.reserve(n_integration_points);
    _secondary_data.N.resize(n_integration_points);

    auto& solid_material = MaterialLib::Solids::selectSolidConstitutiveRelation(
        _process_data.solid_materials, _process_data.material_ids, e.getID());

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        _ip_data.emplace_back(solid_material);
        auto& ip_data = _ip_data[ip];
        auto const& sm = shape_matrices[ip];
        ip_data.N = sm.N;
        ip_data.dNdx = sm.dNdx;
        ip_data.integration_weight =
            _integration_method.getWeightedPoint(ip).getWeight() *
            sm.integralMeasure * sm.detJ;

        // Initialize current time step values
        static const int kelvin_vector_size =
            MathLib::KelvinVector::kelvin_vector_dimensions(DisplacementDim);
        ip_data._sigma.setZero(kelvin_vector_size);
        ip_data._eps.setZero(kelvin_vector_size);

        // Previous time step values are not initialized and are set later.
        ip_data._sigma_prev.resize(kelvin_vector_size);
        ip_data._eps_prev.resize(kelvin_vector_size);

        ip_data._C.resize(kelvin_vector_size, kelvin_vector_size);

        _secondary_data.N[ip] = sm.N;
    }

    for (auto fid : process_data._vec_ele_connected_fractureIDs[e.getID()])
    {
        _fracID_to_local.insert({fid, _fracture_props.size()});
        _fracture_props.push_back(&_process_data.fracture_properties[fid]);
    }

    _junction_props = process_data._vec_ele_connected_junctionIDs[e.getID()] |
                      ranges::views::transform(
                          [&](auto const jid)
                          { return &_process_data.junction_properties[jid]; }) |
                      ranges::to<std::vector>;
}

template <typename ShapeFunction, int DisplacementDim>
void SmallDeformationLocalAssemblerMatrixNearFracture<
    ShapeFunction,
    DisplacementDim>::assembleWithJacobian(double const t, double const dt,
                                           Eigen::VectorXd const& local_u,
                                           Eigen::VectorXd& local_b,
                                           Eigen::MatrixXd& local_J)
{
    assert(_element.getDimension() == DisplacementDim);

    auto const N_DOF_PER_VAR = ShapeFunction::NPOINTS * DisplacementDim;
    auto const n_fractures = _fracture_props.size();
    auto const n_junctions = _junction_props.size();
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
        _integration_method.getNumberOfPoints();

    MPL::VariableArray variables;
    MPL::VariableArray variables_prev;
    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(_element.getID());

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);

        auto& ip_data = _ip_data[ip];
        auto const& w = _ip_data[ip].integration_weight;

        auto const& N = ip_data.N;
        auto const& dNdx = ip_data.dNdx;

        // levelset functions
        Eigen::Vector3d const ip_physical_coords(
            computePhysicalCoordinates(_element, N).data());
        std::vector<double> const levelsets(
            uGlobalEnrichments(_fracture_props, _junction_props,
                               _fracID_to_local, ip_physical_coords));

        // u = u^hat + sum_i(enrich^br_i(x) * [u]_i) + sum_i(enrich^junc_i(x) *
        // [u]_i)
        NodalDisplacementVectorType nodal_total_u = nodal_u;
        for (unsigned i = 0; i < n_enrich_var; i++)
        {
            nodal_total_u += levelsets[i] * vec_nodal_g[i];
        }

        auto const x_coord =
            NumLib::interpolateXCoordinate<ShapeFunction, ShapeMatricesType>(
                _element, N);
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

        variables.mechanical_strain
            .emplace<MathLib::KelvinVector::KelvinVectorType<DisplacementDim>>(
                eps);

        variables_prev.stress
            .emplace<MathLib::KelvinVector::KelvinVectorType<DisplacementDim>>(
                sigma_prev);
        variables_prev.mechanical_strain
            .emplace<MathLib::KelvinVector::KelvinVectorType<DisplacementDim>>(
                eps_prev);
        variables_prev.temperature = _process_data._reference_temperature;

        auto&& solution = _ip_data[ip]._solid_material.integrateStress(
            variables_prev, variables, t, x_position, dt, *state);

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

          int DisplacementDim>
void SmallDeformationLocalAssemblerMatrixNearFracture<ShapeFunction,

                                                      DisplacementDim>::
    computeSecondaryVariableConcreteWithVector(
        double const /*t*/, Eigen::VectorXd const& /*local_x*/)
{
    // Compute average value per element
    using KV = MathLib::KelvinVector::KelvinVectorType<DisplacementDim>;
    KV sigma_avg = KV::Zero();
    auto const e_id = _element.getID();

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        sigma_avg += _ip_data[ip]._sigma;
    }
    sigma_avg /= n_integration_points;

    Eigen::Map<KV>(
        &(*_process_data.element_stresses)[e_id * KV::RowsAtCompileTime]) =
        MathLib::KelvinVector::kelvinVectorToSymmetricTensor(sigma_avg);
}

template <typename ShapeFunction, int DisplacementDim>
std::vector<double> const& SmallDeformationLocalAssemblerMatrixNearFracture<
    ShapeFunction, DisplacementDim>::
    getIntPtSigma(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const
{
    return ProcessLib::getIntegrationPointKelvinVectorData<DisplacementDim>(
        _ip_data, &IntegrationPointDataType::_sigma, cache);
}
template <typename ShapeFunction, int DisplacementDim>
std::vector<double> const& SmallDeformationLocalAssemblerMatrixNearFracture<
    ShapeFunction, DisplacementDim>::
    getIntPtEpsilon(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const
{
    return ProcessLib::getIntegrationPointKelvinVectorData<DisplacementDim>(
        _ip_data, &IntegrationPointDataType::_eps, cache);
}

}  // namespace SmallDeformation
}  // namespace LIE
}  // namespace ProcessLib
