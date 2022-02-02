/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <vector>

#include "HeatConductionProcessData.h"
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
#include "ProcessLib/LocalAssemblerInterface.h"
#include "ProcessLib/LocalAssemblerTraits.h"

namespace ProcessLib
{
namespace HeatConduction
{
const unsigned NUM_NODAL_DOF = 1;

class HeatConductionLocalAssemblerInterface
    : public ProcessLib::LocalAssemblerInterface,
      public NumLib::ExtrapolatableElement
{
public:
    virtual std::vector<double> const& getIntPtHeatFlux(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const = 0;
};

template <typename ShapeFunction, typename IntegrationMethod, int GlobalDim>
class LocalAssemblerData : public HeatConductionLocalAssemblerInterface
{
    using ShapeMatricesType = ShapeMatrixPolicyType<ShapeFunction, GlobalDim>;
    using ShapeMatrices = typename ShapeMatricesType::ShapeMatrices;

    using LocalAssemblerTraits = ProcessLib::LocalAssemblerTraits<
        ShapeMatricesType, ShapeFunction::NPOINTS, NUM_NODAL_DOF, GlobalDim>;

    using NodalMatrixType = typename LocalAssemblerTraits::LocalMatrix;
    using NodalVectorType = typename LocalAssemblerTraits::LocalVector;
    using GlobalDimVectorType = typename ShapeMatricesType::GlobalDimVectorType;

public:
    /// The thermal_conductivity factor is directly integrated into the local
    /// element matrix.
    LocalAssemblerData(MeshLib::Element const& element,
                       std::size_t const local_matrix_size,
                       bool is_axially_symmetric,
                       unsigned const integration_order,
                       HeatConductionProcessData const& process_data)
        : _element(element),
          _process_data(process_data),
          _integration_method(integration_order),
          _shape_matrices(
              NumLib::initShapeMatrices<ShapeFunction, ShapeMatricesType,
                                        GlobalDim>(
                  element, is_axially_symmetric, _integration_method)),
          _heat_fluxes(
              GlobalDim,
              std::vector<double>(_integration_method.getNumberOfPoints()))
    {
        // This assertion is valid only if all nodal d.o.f. use the same shape
        // matrices.
        assert(local_matrix_size == ShapeFunction::NPOINTS * NUM_NODAL_DOF);
        (void)local_matrix_size;
    }

    void assemble(double const t, double const dt,
                  std::vector<double> const& local_x,
                  std::vector<double> const& /*local_xdot*/,
                  std::vector<double>& local_M_data,
                  std::vector<double>& local_K_data,
                  std::vector<double>& /*local_b_data*/) override
    {
        auto const local_matrix_size = local_x.size();
        // This assertion is valid only if all nodal d.o.f. use the same shape
        // matrices.
        assert(local_matrix_size == ShapeFunction::NPOINTS * NUM_NODAL_DOF);

        auto local_M = MathLib::createZeroedMatrix<NodalMatrixType>(
            local_M_data, local_matrix_size, local_matrix_size);
        auto local_K = MathLib::createZeroedMatrix<NodalMatrixType>(
            local_K_data, local_matrix_size, local_matrix_size);

        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        ParameterLib::SpatialPosition pos;
        pos.setElementID(_element.getID());

        auto const& medium =
            *_process_data.media_map->getMedium(_element.getID());
        MaterialPropertyLib::VariableArray vars;

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            pos.setIntegrationPoint(ip);
            auto const& sm = _shape_matrices[ip];
            auto const& wp = _integration_method.getWeightedPoint(ip);

            // get the local temperature and put it in the variable array for
            // access in MPL
            double T_int_pt = 0.0;
            NumLib::shapeFunctionInterpolate(local_x, sm.N, T_int_pt);
            vars[static_cast<int>(MaterialPropertyLib::Variable::temperature)] =
                T_int_pt;

            auto const k = MaterialPropertyLib::formEigenTensor<GlobalDim>(
                medium
                    .property(
                        MaterialPropertyLib::PropertyType::thermal_conductivity)
                    .value(vars, pos, t, dt));
            auto const specific_heat_capacity =
                medium
                    .property(MaterialPropertyLib::PropertyType::
                                  specific_heat_capacity)
                    .template value<double>(vars, pos, t, dt);
            auto const density =
                medium.property(MaterialPropertyLib::PropertyType::density)
                    .template value<double>(vars, pos, t, dt);

            local_K.noalias() += sm.dNdx.transpose() * k * sm.dNdx * sm.detJ *
                                 wp.getWeight() * sm.integralMeasure;
            local_M.noalias() += sm.N.transpose() * density *
                                 specific_heat_capacity * sm.N * sm.detJ *
                                 wp.getWeight() * sm.integralMeasure;
        }
        if (_process_data.mass_lumping)
        {
            local_M = local_M.colwise().sum().eval().asDiagonal();
        }
    }

    void assembleWithJacobian(double const t, double const dt,
                              std::vector<double> const& local_x,
                              std::vector<double> const& local_xdot,
                              std::vector<double>& /*local_M_data*/,
                              std::vector<double>& /*local_K_data*/,
                              std::vector<double>& local_rhs_data,
                              std::vector<double>& local_Jac_data) override
    {
        auto const local_matrix_size = local_x.size();
        // This assertion is valid only if all nodal d.o.f. use the same shape
        // matrices.
        assert(local_matrix_size == ShapeFunction::NPOINTS * NUM_NODAL_DOF);

        auto x = Eigen::Map<NodalVectorType const>(local_x.data(),
                                                   local_matrix_size);

        auto x_dot = Eigen::Map<NodalVectorType const>(local_xdot.data(),
                                                       local_matrix_size);

        auto local_Jac = MathLib::createZeroedMatrix<NodalMatrixType>(
            local_Jac_data, local_matrix_size, local_matrix_size);
        auto local_rhs = MathLib::createZeroedVector<NodalVectorType>(
            local_rhs_data, local_matrix_size);

        NodalMatrixType laplace =
            NodalMatrixType::Zero(local_matrix_size, local_matrix_size);
        NodalMatrixType storage =
            NodalMatrixType::Zero(local_matrix_size, local_matrix_size);

        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        ParameterLib::SpatialPosition pos;
        pos.setElementID(_element.getID());

        auto const& medium =
            *_process_data.media_map->getMedium(_element.getID());
        MaterialPropertyLib::VariableArray vars;

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            pos.setIntegrationPoint(ip);
            auto const& sm = _shape_matrices[ip];
            double const w =
                _integration_method.getWeightedPoint(ip).getWeight() * sm.detJ *
                sm.integralMeasure;

            // get the local temperature and put it in the variable array for
            // access in MPL
            double T_int_pt = 0.0;
            NumLib::shapeFunctionInterpolate(local_x, sm.N, T_int_pt);
            vars[static_cast<int>(MaterialPropertyLib::Variable::temperature)] =
                T_int_pt;

            auto const k = MaterialPropertyLib::formEigenTensor<GlobalDim>(
                medium
                    .property(
                        MaterialPropertyLib::PropertyType::thermal_conductivity)
                    .value(vars, pos, t, dt));
            auto const specific_heat_capacity =
                medium
                    .property(MaterialPropertyLib::PropertyType::
                                  specific_heat_capacity)
                    .template value<double>(vars, pos, t, dt);
            auto const density =
                medium.property(MaterialPropertyLib::PropertyType::density)
                    .template value<double>(vars, pos, t, dt);

            laplace.noalias() += sm.dNdx.transpose() * k * sm.dNdx * w;
            storage.noalias() +=
                sm.N.transpose() * density * specific_heat_capacity * sm.N * w;
        }
        if (_process_data.mass_lumping)
        {
            storage = storage.colwise().sum().eval().asDiagonal();
        }

        local_Jac.noalias() += laplace + storage / dt;
        local_rhs.noalias() -= laplace * x + storage * x_dot;
    }

    void computeSecondaryVariableConcrete(
        double const t, double const dt, Eigen::VectorXd const& local_x,
        Eigen::VectorXd const& /*local_x_dot*/) override
    {
        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        ParameterLib::SpatialPosition pos;
        pos.setElementID(_element.getID());

        auto const& medium =
            *_process_data.media_map->getMedium(_element.getID());
        MaterialPropertyLib::VariableArray vars;

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            pos.setIntegrationPoint(ip);
            auto const& sm = _shape_matrices[ip];
            // get the local temperature and put it in the variable array for
            // access in MPL
            double T_int_pt = 0.0;
            NumLib::shapeFunctionInterpolate(local_x, sm.N, T_int_pt);
            vars[static_cast<int>(MaterialPropertyLib::Variable::temperature)] =
                T_int_pt;

            auto const k = MaterialPropertyLib::formEigenTensor<GlobalDim>(
                medium
                    .property(
                        MaterialPropertyLib::PropertyType::thermal_conductivity)
                    .value(vars, pos, t, dt));
            // heat flux only computed for output.
            GlobalDimVectorType const heat_flux = -k * sm.dNdx * local_x;

            for (int d = 0; d < GlobalDim; ++d)
            {
                _heat_fluxes[d][ip] = heat_flux[d];
            }
        }
    }

    Eigen::Map<const Eigen::RowVectorXd> getShapeMatrix(
        const unsigned integration_point) const override
    {
        auto const& N = _shape_matrices[integration_point].N;

        // assumes N is stored contiguously in memory
        return Eigen::Map<const Eigen::RowVectorXd>(N.data(), N.size());
    }

    std::vector<double> const& getIntPtHeatFlux(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const override
    {
        int const process_id = 0;  // monolithic case.
        auto const indices =
            NumLib::getIndices(_element.getID(), *dof_table[process_id]);
        assert(!indices.empty());
        auto const& local_x = x[process_id]->get(indices);

        auto const T_nodal_values = Eigen::Map<const NodalVectorType>(
            local_x.data(), ShapeFunction::NPOINTS);

        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        ParameterLib::SpatialPosition pos;
        pos.setElementID(_element.getID());

        auto const& medium =
            *_process_data.media_map->getMedium(_element.getID());
        MaterialPropertyLib::VariableArray vars;

        double const dt = std::numeric_limits<double>::quiet_NaN();
        cache.clear();
        auto cache_mat = MathLib::createZeroedMatrix<
            Eigen::Matrix<double, GlobalDim, Eigen::Dynamic, Eigen::RowMajor>>(
            cache, GlobalDim, n_integration_points);

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            pos.setIntegrationPoint(ip);
            auto const& sm = _shape_matrices[ip];
            // get the local temperature and put it in the variable array for
            // access in MPL
            vars[static_cast<int>(MaterialPropertyLib::Variable::temperature)] =
                sm.N.dot(T_nodal_values);

            auto const k = MaterialPropertyLib::formEigenTensor<GlobalDim>(
                medium
                    .property(
                        MaterialPropertyLib::PropertyType::thermal_conductivity)
                    .value(vars, pos, t, dt));

            // heat flux only computed for output.
            cache_mat.col(ip).noalias() = -k * sm.dNdx * T_nodal_values;
        }

        return cache;
    }

private:
    MeshLib::Element const& _element;
    HeatConductionProcessData const& _process_data;

    IntegrationMethod const _integration_method;
    std::vector<ShapeMatrices, Eigen::aligned_allocator<ShapeMatrices>>
        _shape_matrices;

    std::vector<std::vector<double>> _heat_fluxes;
};

}  // namespace HeatConduction
}  // namespace ProcessLib
