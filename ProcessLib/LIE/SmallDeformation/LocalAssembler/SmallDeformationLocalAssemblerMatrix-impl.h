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
#include <valarray>
#include <vector>

#include "IntegrationPointDataMatrix.h"
#include "MaterialLib/PhysicalConstant.h"
#include "MaterialLib/SolidModels/SelectSolidConstitutiveRelation.h"
#include "MathLib/KelvinVector.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "NumLib/Fem/InitShapeMatrices.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "ProcessLib/Deformation/BMatrixPolicy.h"
#include "ProcessLib/Deformation/LinearBMatrix.h"
#include "ProcessLib/LIE/SmallDeformation/SmallDeformationProcessData.h"
#include "ProcessLib/Utils/SetOrGetIntegrationPointData.h"
#include "SecondaryData.h"
#include "SmallDeformationLocalAssemblerInterface.h"
#include "SmallDeformationLocalAssemblerMatrix.h"

namespace ProcessLib
{
namespace LIE
{
namespace SmallDeformation
{
template <typename ShapeFunction, int DisplacementDim>
SmallDeformationLocalAssemblerMatrix<ShapeFunction, DisplacementDim>::
    SmallDeformationLocalAssemblerMatrix(
        MeshLib::Element const& e,
        std::size_t const /*local_matrix_size*/,
        NumLib::GenericIntegrationMethod const& integration_method,
        bool const is_axially_symmetric,
        SmallDeformationProcessData<DisplacementDim>& process_data)
    : _process_data(process_data),
      _integration_method(integration_method),
      _element(e),
      _is_axially_symmetric(is_axially_symmetric)
{
    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    _ip_data.reserve(n_integration_points);
    _secondary_data.N.resize(n_integration_points);

    auto const shape_matrices =
        NumLib::initShapeMatrices<ShapeFunction, ShapeMatricesType,
                                  DisplacementDim>(e, is_axially_symmetric,
                                                   _integration_method);

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
}

template <typename ShapeFunction, int DisplacementDim>
void SmallDeformationLocalAssemblerMatrix<ShapeFunction, DisplacementDim>::
    assembleWithJacobian(double const t, double const dt,
                         std::vector<double> const& local_x,
                         std::vector<double> const& /*local_x_prev*/,
                         std::vector<double>& /*local_M_data*/,
                         std::vector<double>& /*local_K_data*/,
                         std::vector<double>& local_b_data,
                         std::vector<double>& local_Jac_data)
{
    assert(_element.getDimension() == DisplacementDim);

    auto const local_matrix_size = local_x.size();

    auto local_Jac = MathLib::createZeroedMatrix<StiffnessMatrixType>(
        local_Jac_data, local_matrix_size, local_matrix_size);

    auto local_b = MathLib::createZeroedVector<NodalDisplacementVectorType>(
        local_b_data, local_matrix_size);

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    MPL::VariableArray variables;
    MPL::VariableArray variables_prev;
    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(_element.getID());

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);
        auto const& w = _ip_data[ip].integration_weight;

        auto const& N = _ip_data[ip].N;
        auto const& dNdx = _ip_data[ip].dNdx;
        auto const x_coord =
            NumLib::interpolateXCoordinate<ShapeFunction, ShapeMatricesType>(
                _element, N);
        auto const B =
            LinearBMatrix::computeBMatrix<DisplacementDim,
                                          ShapeFunction::NPOINTS,
                                          typename BMatricesType::BMatrixType>(
                dNdx, N, x_coord, _is_axially_symmetric);

        auto const& eps_prev = _ip_data[ip]._eps_prev;
        auto const& sigma_prev = _ip_data[ip]._sigma_prev;

        auto& eps = _ip_data[ip]._eps;
        auto& sigma = _ip_data[ip]._sigma;
        auto& state = _ip_data[ip]._material_state_variables;

        eps.noalias() =
            B * Eigen::Map<typename BMatricesType::NodalForceVectorType const>(
                    local_x.data(), ShapeFunction::NPOINTS * DisplacementDim);

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

        local_b.noalias() -= B.transpose() * sigma * w;
        local_Jac.noalias() += B.transpose() * C * B * w;
    }
}

template <typename ShapeFunction, int DisplacementDim>
void SmallDeformationLocalAssemblerMatrix<ShapeFunction, DisplacementDim>::
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
std::vector<double> const&
SmallDeformationLocalAssemblerMatrix<ShapeFunction, DisplacementDim>::
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
std::vector<double> const&
SmallDeformationLocalAssemblerMatrix<ShapeFunction, DisplacementDim>::
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
