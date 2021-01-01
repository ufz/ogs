/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <vector>

#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Utils/FormEigenTensor.h"
#include "MaterialLib/MPL/VariableType.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "NumLib/DOF/DOFTableUtil.h"
#include "NumLib/Extrapolation/ExtrapolatableElement.h"
#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/InitShapeMatrices.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "NumLib/Function/Interpolation.h"
#include "ParameterLib/Parameter.h"
#include "ProcessLib/LocalAssemblerInterface.h"
#include "ProcessLib/LocalAssemblerTraits.h"
#include "SteadyStateDiffusionData.h"

namespace ProcessLib
{
namespace SteadyStateDiffusion
{
const unsigned NUM_NODAL_DOF = 1;

class SteadyStateDiffusionLocalAssemblerInterface
    : public ProcessLib::LocalAssemblerInterface,
      public NumLib::ExtrapolatableElement
{
public:
    virtual std::vector<double> const& getIntPtDarcyVelocity(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const = 0;
};

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
class LocalAssemblerData : public SteadyStateDiffusionLocalAssemblerInterface
{
    using ShapeMatricesType = ShapeMatrixPolicyType<ShapeFunction, GlobalDim>;
    using ShapeMatrices = typename ShapeMatricesType::ShapeMatrices;

    using LocalAssemblerTraits = ProcessLib::LocalAssemblerTraits<
        ShapeMatricesType, ShapeFunction::NPOINTS, NUM_NODAL_DOF, GlobalDim>;

    using NodalMatrixType = typename LocalAssemblerTraits::LocalMatrix;
    using NodalVectorType = typename LocalAssemblerTraits::LocalVector;
    using GlobalDimVectorType = typename ShapeMatricesType::GlobalDimVectorType;

public:
    /// The hydraulic_conductivity factor is directly integrated into the local
    /// element matrix.
    LocalAssemblerData(MeshLib::Element const& element,
                       std::size_t const /*local_matrix_size*/,
                       bool is_axially_symmetric,
                       unsigned const integration_order,
                       SteadyStateDiffusionData const& process_data)
        : _element(element),
          _process_data(process_data),
          _integration_method(integration_order),
          _shape_matrices(
              NumLib::initShapeMatrices<ShapeFunction, ShapeMatricesType,
                                        GlobalDim>(
                  element, is_axially_symmetric, _integration_method))
    {
    }

    void assemble(double const t, double const dt,
                  std::vector<double> const& local_x,
                  std::vector<double> const& /*local_xdot*/,
                  std::vector<double>& /*local_M_data*/,
                  std::vector<double>& local_K_data,
                  std::vector<double>& /*local_b_data*/) override
    {
        auto const local_matrix_size = local_x.size();
        // This assertion is valid only if all nodal d.o.f. use the same shape
        // matrices.
        assert(local_matrix_size == ShapeFunction::NPOINTS * NUM_NODAL_DOF);

        auto local_K = MathLib::createZeroedMatrix<NodalMatrixType>(
            local_K_data, local_matrix_size, local_matrix_size);

        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        ParameterLib::SpatialPosition pos;
        pos.setElementID(_element.getID());

        auto const& medium =
            *_process_data.media_map->getMedium(_element.getID());
        MaterialPropertyLib::VariableArray vars;
        vars[static_cast<int>(MaterialPropertyLib::Variable::temperature)] =
            medium
                .property(
                    MaterialPropertyLib::PropertyType::reference_temperature)
                .template value<double>(vars, pos, t, dt);

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            pos.setIntegrationPoint(ip);
            auto const& sm = _shape_matrices[ip];
            auto const& wp = _integration_method.getWeightedPoint(ip);

            double p_int_pt = 0.0;
            NumLib::shapeFunctionInterpolate(local_x, sm.N, p_int_pt);
            vars[static_cast<int>(
                MaterialPropertyLib::Variable::phase_pressure)] = p_int_pt;
            auto const k = MaterialPropertyLib::formEigenTensor<GlobalDim>(
                medium.property(MaterialPropertyLib::PropertyType::diffusion)
                    .value(vars, pos, t, dt));

            local_K.noalias() += sm.dNdx.transpose() * k * sm.dNdx * sm.detJ *
                                 sm.integralMeasure * wp.getWeight();
        }
    }

    /// Computes the flux in the point \c p_local_coords that is given in local
    /// coordinates using the values from \c local_x.
    Eigen::Vector3d getFlux(MathLib::Point3d const& p_local_coords,
                            double const t,
                            std::vector<double> const& local_x) const override
    {
        // TODO (tf) Temporary value not used by current material models. Need
        // extension of getFlux interface
        double const dt = std::numeric_limits<double>::quiet_NaN();

        // Eval shape matrices at given point
        // Note: Axial symmetry is set to false here, because we only need dNdx
        // here, which is not affected by axial symmetry.
        auto const shape_matrices =
            NumLib::computeShapeMatrices<ShapeFunction, ShapeMatricesType,
                                         GlobalDim>(
                _element, false /*is_axially_symmetric*/,
                std::array{p_local_coords})[0];

        // fetch hydraulic conductivity
        ParameterLib::SpatialPosition pos;
        pos.setElementID(_element.getID());
        auto const& medium =
            *_process_data.media_map->getMedium(_element.getID());

        MaterialPropertyLib::VariableArray vars;
        vars[static_cast<int>(MaterialPropertyLib::Variable::temperature)] =
            medium
                .property(
                    MaterialPropertyLib::PropertyType::reference_temperature)
                .template value<double>(vars, pos, t, dt);
        double pressure = 0.0;
        NumLib::shapeFunctionInterpolate(local_x, shape_matrices.N, pressure);
        vars[static_cast<int>(MaterialPropertyLib::Variable::phase_pressure)] =
            pressure;

        auto const k = MaterialPropertyLib::formEigenTensor<GlobalDim>(
            medium.property(MaterialPropertyLib::PropertyType::diffusion)
                .value(vars, pos, t, dt));

        Eigen::Vector3d flux(0.0, 0.0, 0.0);
        flux.head<GlobalDim>() =
            -k * shape_matrices.dNdx *
            Eigen::Map<const NodalVectorType>(local_x.data(), local_x.size());

        return flux;
    }

    Eigen::Map<const Eigen::RowVectorXd> getShapeMatrix(
        const unsigned integration_point) const override
    {
        auto const& N = _shape_matrices[integration_point].N;

        // assumes N is stored contiguously in memory
        return Eigen::Map<const Eigen::RowVectorXd>(N.data(), N.size());
    }

    std::vector<double> const& getIntPtDarcyVelocity(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const override
    {
        // TODO (tf) Temporary value not used by current material models. Need
        // extension of secondary variable interface.
        double const dt = std::numeric_limits<double>::quiet_NaN();

        auto const n_integration_points =
            _integration_method.getNumberOfPoints();

        int const process_id = 0;  // monolithic scheme
        auto const indices =
            NumLib::getIndices(_element.getID(), *dof_table[process_id]);
        assert(!indices.empty());
        auto const local_x = x[process_id]->get(indices);
        auto const local_x_vec =
            MathLib::toVector<Eigen::Matrix<double, ShapeFunction::NPOINTS, 1>>(
                local_x, ShapeFunction::NPOINTS);

        cache.clear();
        auto cache_mat = MathLib::createZeroedMatrix<
            Eigen::Matrix<double, GlobalDim, Eigen::Dynamic, Eigen::RowMajor>>(
            cache, GlobalDim, n_integration_points);

        ParameterLib::SpatialPosition pos;
        pos.setElementID(_element.getID());

        auto const& medium =
            *_process_data.media_map->getMedium(_element.getID());

        MaterialPropertyLib::VariableArray vars;
        vars[static_cast<int>(MaterialPropertyLib::Variable::temperature)] =
            medium
                .property(
                    MaterialPropertyLib::PropertyType::reference_temperature)
                .template value<double>(vars, pos, t, dt);
        double pressure = 0.0;
        for (unsigned i = 0; i < n_integration_points; ++i)
        {
            pos.setIntegrationPoint(i);
            NumLib::shapeFunctionInterpolate(local_x, _shape_matrices[i].N,
                                             pressure);
            vars[static_cast<int>(
                MaterialPropertyLib::Variable::phase_pressure)] = pressure;

            auto const k = MaterialPropertyLib::formEigenTensor<GlobalDim>(
                medium.property(MaterialPropertyLib::PropertyType::diffusion)
                    .value(vars, pos, t, dt));
            // dimensions: (d x 1) = (d x n) * (n x 1)
            cache_mat.col(i).noalias() =
                -k * _shape_matrices[i].dNdx * local_x_vec;
        }

        return cache;
    }

private:
    MeshLib::Element const& _element;
    SteadyStateDiffusionData const& _process_data;

    IntegrationMethod const _integration_method;
    std::vector<ShapeMatrices, Eigen::aligned_allocator<ShapeMatrices>>
        _shape_matrices;
};

}  // namespace SteadyStateDiffusion
}  // namespace ProcessLib
