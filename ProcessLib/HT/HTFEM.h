/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_HTPROCESS_FEM_H_
#define PROCESS_LIB_HTPROCESS_FEM_H_

#include <iostream>
#include <Eigen/Dense>
#include <vector>


#include "HTProcessData.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "NumLib/Extrapolation/ExtrapolatableElement.h"
#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "NumLib/Function/Interpolation.h"
#include "ProcessLib/LocalAssemblerInterface.h"
#include "ProcessLib/LocalAssemblerTraits.h"
#include "ProcessLib/Parameter/Parameter.h"
#include "ProcessLib/Utils/InitShapeMatrices.h"

namespace ProcessLib
{
namespace HT
{
const unsigned NUM_NODAL_DOF = 2;

class HTLocalAssemblerInterface
    : public ProcessLib::LocalAssemblerInterface,
      public NumLib::ExtrapolatableElement
{
public:
    virtual std::vector<double> const& getIntPtDarcyVelocityX(
        std::vector<double>& /*cache*/) const = 0;

    virtual std::vector<double> const& getIntPtDarcyVelocityY(
        std::vector<double>& /*cache*/) const = 0;

    virtual std::vector<double> const& getIntPtDarcyVelocityZ(
        std::vector<double>& /*cache*/) const = 0;
};

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
class LocalAssemblerData : public HTLocalAssemblerInterface
{
    using ShapeMatricesType = ShapeMatrixPolicyType<ShapeFunction, GlobalDim>;
    using ShapeMatrices = typename ShapeMatricesType::ShapeMatrices;

    using LocalAssemblerTraits = ProcessLib::LocalAssemblerTraits<
        ShapeMatricesType, ShapeFunction::NPOINTS, NUM_NODAL_DOF, GlobalDim>;
    using NodalMatrixType = typename LocalAssemblerTraits::LocalMatrix;
    using NodalVectorType = typename LocalAssemblerTraits::LocalVector;
    using GlobalDimVectorType = typename ShapeMatricesType::GlobalDimVectorType;

public:
    LocalAssemblerData(MeshLib::Element const& element,
                       std::size_t const local_matrix_size,
                       bool is_axially_symmetric,
                       unsigned const integration_order,
                       HTProcessData const& process_data)
        : _element(element),
          _process_data(process_data),
          _integration_method(integration_order),
          _shape_matrices(initShapeMatrices<ShapeFunction, ShapeMatricesType,
                                            IntegrationMethod, GlobalDim>(
              element, is_axially_symmetric, _integration_method)),
          _darcy_velocities(
              GlobalDim,
              std::vector<double>(_integration_method.getNumberOfPoints()))
    {
        // This assertion is valid only if all nodal d.o.f. use the same shape
        // matrices.
        assert(local_matrix_size == ShapeFunction::NPOINTS * NUM_NODAL_DOF);
        (void)local_matrix_size;
    }

    void assemble(double const t, std::vector<double> const& local_x,
                  std::vector<double>& local_M_data,
                  std::vector<double>& local_K_data,
                  std::vector<double>& local_b_data) override
    {
        auto const local_matrix_size = local_x.size();
        // This assertion is valid only if all nodal d.o.f. use the same shape
        // matrices.
        assert(local_matrix_size == ShapeFunction::NPOINTS * NUM_NODAL_DOF);

        auto local_M = MathLib::createZeroedMatrix<NodalMatrixType>(
            local_M_data, local_matrix_size, local_matrix_size);
        auto local_K = MathLib::createZeroedMatrix<NodalMatrixType>(
            local_K_data, local_matrix_size, local_matrix_size);
        auto local_b = MathLib::createZeroedVector<NodalVectorType>(
            local_b_data, local_matrix_size);

        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        SpatialPosition pos;
        pos.setElementID(_element.getID());

        auto const num_nodes = ShapeFunction::NPOINTS;
        auto p_nodal_values =
            Eigen::Map<const Eigen::VectorXd>(&local_x[num_nodes], num_nodes);

        auto const & b = _process_data.specific_body_force.head(GlobalDim);

        for (std::size_t ip(0); ip < n_integration_points; ip++)
        {
            pos.setIntegrationPoint(ip);
            auto const temperature0 =
                _process_data.reference_temperature_fluid_density_model(t,
                                                                        pos)[0];
            auto const beta =
                _process_data.thermal_expansion_coefficient(t, pos)[0];

            auto const density_fluid = _process_data.density_fluid(t, pos)[0];
            auto const density_solid = _process_data.density_solid(t, pos)[0];
            auto const specific_storage =
                _process_data.specific_storage(t, pos)[0];
            auto const intrinsic_permeability =
                _process_data.intrinsic_permeability(t, pos)[0];
            auto const viscosity0 = _process_data.viscosity(t, pos)[0];

            auto const porosity = _process_data.porosity(t, pos)[0];

            auto const specific_heat_capacity_solid =
                _process_data.specific_heat_capacity_solid(t, pos)[0];
            auto const specific_heat_capacity_fluid =
                _process_data.specific_heat_capacity_fluid(t, pos)[0];
            double const heat_capacity =
                density_solid * specific_heat_capacity_solid * (1 - porosity) +
                density_fluid * specific_heat_capacity_fluid * porosity;

            auto const thermal_conductivity_solid =
                _process_data.thermal_conductivity_solid(t, pos)[0];
            auto const thermal_conductivity_fluid =
                _process_data.thermal_conductivity_fluid(t, pos)[0];
            double const thermal_conductivity =
                thermal_conductivity_solid * (1 - porosity) +
                thermal_conductivity_fluid * porosity;

            auto const& sm = _shape_matrices[ip];
            auto const& wp = _integration_method.getWeightedPoint(ip);
            auto Ktt = local_K.template block<num_nodes, num_nodes>(0, 0);
            auto Mtt = local_M.template block<num_nodes, num_nodes>(0, 0);
            auto Kpp = local_K.template block<num_nodes, num_nodes>(num_nodes,
                                                                    num_nodes);
            auto Mpp = local_M.template block<num_nodes, num_nodes>(num_nodes,
                                                                    num_nodes);
            auto Bp = local_b.template block<num_nodes, 1>(num_nodes, 0);

            double T_int_pt = 0.0;
            double p_int_pt = 0.0;
            // Order matters: First T, then P!
            NumLib::shapeFunctionInterpolate(local_x, sm.N, T_int_pt, p_int_pt);

            double const delta_T(T_int_pt - temperature0);
            // TODO include this via material lib
            double density_water_T = density_fluid * (1 - beta * delta_T);
            // TODO include viscosity computations via material lib
            // double const viscosity = viscosity0 * std::exp(- delta_T/75.0);
            double const perm_visc =
                intrinsic_permeability / viscosity0;

            Eigen::Matrix<double, -1, 1, 0, -1, 1> const velocity =
                -perm_visc * (sm.dNdx * p_nodal_values - density_water_T * b);

            auto const integral_term =
                sm.integralMeasure * sm.detJ * wp.getWeight();
            // matrix assembly
            Ktt.noalias() +=
                integral_term *
                (sm.dNdx.transpose() * thermal_conductivity * sm.dNdx +
                 sm.N.transpose() * velocity.transpose() * sm.dNdx *
                     density_fluid * specific_heat_capacity_fluid);
            Kpp.noalias() +=
                integral_term * sm.dNdx.transpose() * perm_visc * sm.dNdx;
            Mtt.noalias() +=
                integral_term * sm.N.transpose() * heat_capacity * sm.N;
            Mpp.noalias() +=
                integral_term * sm.N.transpose() * specific_storage * sm.N;
            Bp += perm_visc * integral_term * sm.dNdx.transpose() * b *
                  density_water_T;
            /* with Oberbeck-Boussing assumption density difference only exists
             * in buoyancy effects */

            // velocity computed for output.
            for (unsigned d = 0; d < GlobalDim; ++d)
            {
                _darcy_velocities[d][ip] = velocity[d];
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

    std::vector<double> const& getIntPtDarcyVelocityX(
        std::vector<double>& /*cache*/) const override
    {
        assert(_darcy_velocities.size() > 0);
        return _darcy_velocities[0];
    }

    std::vector<double> const& getIntPtDarcyVelocityY(
        std::vector<double>& /*cache*/) const override
    {
        assert(_darcy_velocities.size() > 1);
        return _darcy_velocities[1];
    }

    std::vector<double> const& getIntPtDarcyVelocityZ(
        std::vector<double>& /*cache*/) const override
    {
        assert(_darcy_velocities.size() > 2);
        return _darcy_velocities[2];
    }

private:
    MeshLib::Element const& _element;
    HTProcessData const& _process_data;

    IntegrationMethod const _integration_method;
    std::vector<ShapeMatrices> _shape_matrices;
    std::vector<std::vector<double>> _darcy_velocities;
};

}  // namespace HT
}  // namespace ProcessLib

#endif  // PROCESS_LIB_HT_FEM_H_
