/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

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

    using NodalMatrixType = typename ShapeMatricesType::NodalMatrixType;
    using NodalVectorType = typename ShapeMatricesType::NodalVectorType;
        using NodalRowVectorType = typename ShapeMatricesType::NodalRowVectorType;
    using GlobalDimVectorType = typename ShapeMatricesType::GlobalDimVectorType;
    using GlobalDimMatrixType = typename ShapeMatricesType::GlobalDimMatrixType;

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
            Eigen::Map<const NodalVectorType>(&local_x[num_nodes], num_nodes);

        auto const & b = _process_data.specific_body_force.head(GlobalDim);

        MaterialLib::Fluid::FluidProperty::ArrayType vars;

        for (std::size_t ip(0); ip < n_integration_points; ip++)
        {
            pos.setIntegrationPoint(ip);

            auto const fluid_reference_density =
                _process_data.fluid_reference_density(t, pos)[0];

            auto const density_solid = _process_data.density_solid(t, pos)[0];
            // TODO the argument to getValue() has to be changed for non
            // constant storage model
            auto const specific_storage =
                _process_data.specific_storage_model->getValue(0.0);
            auto const& intrinsic_permeability =
                _process_data.permeability_model;

            auto const thermal_conductivity_solid =
                _process_data.thermal_conductivity_solid(t, pos)[0];
            auto const thermal_conductivity_fluid =
                _process_data.thermal_conductivity_fluid(t, pos)[0];

            auto const& sm = _shape_matrices[ip];
            double T_int_pt = 0.0;
            double p_int_pt = 0.0;
            // Order matters: First T, then P!
            NumLib::shapeFunctionInterpolate(local_x, sm.N, T_int_pt, p_int_pt);

            // TODO the first argument has to be changed for non constant
            // porosity model
            auto const porosity =
                _process_data.porosity_model->getValue(0.0, T_int_pt);

            double const thermal_conductivity =
                thermal_conductivity_solid * (1 - porosity) +
                 thermal_conductivity_fluid * porosity;

            auto const specific_heat_capacity_solid =
                _process_data.specific_heat_capacity_solid(t, pos)[0];
            auto const specific_heat_capacity_fluid =
                _process_data.specific_heat_capacity_fluid(t, pos)[0];

            auto const thermal_dispersivity_longitudinal =
                _process_data.thermal_dispersivity_longitudinal(t, pos)[0];
            auto const thermal_dispersivity_transversal =
                _process_data.thermal_dispersivity_transversal(t, pos)[0];

            auto const& wp = _integration_method.getWeightedPoint(ip);
            auto Ktt = local_K.template block<num_nodes, num_nodes>(0, 0);
            auto Mtt = local_M.template block<num_nodes, num_nodes>(0, 0);
            auto Kpp = local_K.template block<num_nodes, num_nodes>(num_nodes,
                                                                    num_nodes);
            auto Mpp = local_M.template block<num_nodes, num_nodes>(num_nodes,
                                                                    num_nodes);
            auto Bp = local_b.template block<num_nodes, 1>(num_nodes, 0);

            // Use the fluid density model to compute the density
            vars[static_cast<int>(
                MaterialLib::Fluid::PropertyVariableType::T)] = T_int_pt;
            vars[static_cast<int>(
                MaterialLib::Fluid::PropertyVariableType::p)] = p_int_pt;
            auto const density_water_T =
                _process_data.fluid_density->getValue(vars);

            // Use the viscosity model to compute the viscosity
            auto const viscosity =
                _process_data.viscosity_model->getValue(vars);
            GlobalDimMatrixType perm_over_visc =
                intrinsic_permeability / viscosity;

            GlobalDimVectorType const velocity =
                -perm_over_visc *
                (sm.dNdx * p_nodal_values - density_water_T * b);

            double const velocity_magnitude = velocity.norm();
            GlobalDimMatrixType const& I(
                GlobalDimMatrixType::Identity(GlobalDim, GlobalDim));
            GlobalDimMatrixType thermal_dispersivity =
                fluid_reference_density * specific_heat_capacity_fluid *
                (thermal_dispersivity_transversal * velocity_magnitude *
                     I +
                 (thermal_dispersivity_longitudinal -
                  thermal_dispersivity_transversal) /
                     velocity_magnitude * velocity * velocity.transpose());

            auto const hydrodynamic_thermodispersion =
                thermal_conductivity * I + thermal_dispersivity;

            double const heat_capacity =
                density_solid * specific_heat_capacity_solid * (1 - porosity) +
                fluid_reference_density * specific_heat_capacity_fluid * porosity;

            auto const integral_term =
                sm.integralMeasure * sm.detJ * wp.getWeight();
            // matrix assembly
            Ktt.noalias() +=
                integral_term *
                (sm.dNdx.transpose() * hydrodynamic_thermodispersion * sm.dNdx +
                 sm.N.transpose() * velocity.transpose() * sm.dNdx *
                     fluid_reference_density * specific_heat_capacity_fluid);
            Kpp.noalias() +=
                integral_term * sm.dNdx.transpose() * perm_over_visc * sm.dNdx;
            Mtt.noalias() +=
                integral_term * sm.N.transpose() * heat_capacity * sm.N;
            Mpp.noalias() +=
                integral_term * sm.N.transpose() * specific_storage * sm.N;
            Bp += integral_term * density_water_T * sm.dNdx.transpose() *
                  perm_over_visc * b;
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
