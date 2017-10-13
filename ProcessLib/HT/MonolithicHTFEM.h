/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file   MonolithicHTFEM.h
 *  Created on October 11, 2017, 2:33 PM
 */

#pragma once

#include <Eigen/Dense>
#include <vector>

#include "HTMaterialProperties.h"
#include "NumLib/DOF/DOFTableUtil.h"
#include "NumLib/Extrapolation/ExtrapolatableElement.h"
#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "ProcessLib/Parameter/Parameter.h"
#include "ProcessLib/Utils/InitShapeMatrices.h"

#include "HTFEM.h"

namespace ProcessLib
{
namespace HT
{
const unsigned NUM_NODAL_DOF = 2;

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
class MonolithicHTFEM
    : public HTFEM<ShapeFunction, IntegrationMethod, GlobalDim>
{
    using ShapeMatricesType = ShapeMatrixPolicyType<ShapeFunction, GlobalDim>;
    using ShapeMatrices = typename ShapeMatricesType::ShapeMatrices;

    using LocalMatrixType = typename ShapeMatricesType::template MatrixType<
        NUM_NODAL_DOF * ShapeFunction::NPOINTS,
        NUM_NODAL_DOF * ShapeFunction::NPOINTS>;
    using LocalVectorType =
        typename ShapeMatricesType::template VectorType<NUM_NODAL_DOF *
                                                        ShapeFunction::NPOINTS>;

    using NodalVectorType = typename ShapeMatricesType::NodalVectorType;
    using NodalRowVectorType = typename ShapeMatricesType::NodalRowVectorType;

    using GlobalDimVectorType = typename ShapeMatricesType::GlobalDimVectorType;
    using GlobalDimMatrixType = typename ShapeMatricesType::GlobalDimMatrixType;

public:
    MonolithicHTFEM(MeshLib::Element const& element,
                    std::size_t const local_matrix_size,
                    bool is_axially_symmetric,
                    unsigned const integration_order,
                    HTMaterialProperties const& material_properties)
        : HTFEM<ShapeFunction, IntegrationMethod, GlobalDim>(
              element, local_matrix_size, is_axially_symmetric,
              integration_order, material_properties, NUM_NODAL_DOF)
    {
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

        auto local_M = MathLib::createZeroedMatrix<LocalMatrixType>(
            local_M_data, local_matrix_size, local_matrix_size);
        auto local_K = MathLib::createZeroedMatrix<LocalMatrixType>(
            local_K_data, local_matrix_size, local_matrix_size);
        auto local_b = MathLib::createZeroedVector<LocalVectorType>(
            local_b_data, local_matrix_size);

        auto const num_nodes = ShapeFunction::NPOINTS;

        auto Ktt = local_K.template block<num_nodes, num_nodes>(0, 0);
        auto Mtt = local_M.template block<num_nodes, num_nodes>(0, 0);
        auto Kpp =
            local_K.template block<num_nodes, num_nodes>(num_nodes, num_nodes);
        auto Mpp =
            local_M.template block<num_nodes, num_nodes>(num_nodes, num_nodes);
        auto Bp = local_b.template block<num_nodes, 1>(num_nodes, 0);

        SpatialPosition pos;
        pos.setElementID(this->_element.getID());

        auto p_nodal_values =
            Eigen::Map<const NodalVectorType>(&local_x[num_nodes], num_nodes);

        auto const& process_data = this->_material_properties;

        auto const& b = process_data.specific_body_force;

        GlobalDimMatrixType const& I(
            GlobalDimMatrixType::Identity(GlobalDim, GlobalDim));

        MaterialLib::Fluid::FluidProperty::ArrayType vars;

        unsigned const n_integration_points =
            this->_integration_method.getNumberOfPoints();

        for (std::size_t ip(0); ip < n_integration_points; ip++)
        {
            pos.setIntegrationPoint(ip);

            auto const fluid_reference_density =
                process_data.fluid_reference_density(t, pos)[0];

            auto const density_solid = process_data.density_solid(t, pos)[0];
            // \todo the argument to getValue() has to be changed for non
            // constant storage model
            auto const specific_storage =
                process_data.porous_media_properties.getSpecificStorage(t, pos)
                    .getValue(0.0);

            auto const thermal_conductivity_solid =
                process_data.thermal_conductivity_solid(t, pos)[0];
            auto const thermal_conductivity_fluid =
                process_data.thermal_conductivity_fluid(t, pos)[0];

            auto const& ip_data = this->_ip_data[ip];
            auto const& N = ip_data.N;
            auto const& dNdx = ip_data.dNdx;
            auto const& w = ip_data.integration_weight;

            double T_int_pt = 0.0;
            double p_int_pt = 0.0;
            // Order matters: First T, then P!
            NumLib::shapeFunctionInterpolate(local_x, N, T_int_pt, p_int_pt);

            // \todo the first argument has to be changed for non constant
            // porosity model
            auto const porosity =
                process_data.porous_media_properties.getPorosity(t, pos)
                    .getValue(t, pos, 0.0, T_int_pt);
            auto const intrinsic_permeability =
                process_data.porous_media_properties.getIntrinsicPermeability(
                    t, pos).getValue(t, pos, 0.0, T_int_pt);

            double const thermal_conductivity =
                thermal_conductivity_solid * (1 - porosity) +
                thermal_conductivity_fluid * porosity;

            auto const specific_heat_capacity_solid =
                process_data.specific_heat_capacity_solid(t, pos)[0];

            vars[static_cast<int>(
                MaterialLib::Fluid::PropertyVariableType::T)] = T_int_pt;
            vars[static_cast<int>(
                MaterialLib::Fluid::PropertyVariableType::p)] = p_int_pt;
            auto const specific_heat_capacity_fluid =
                process_data.fluid_properties->getValue(
                    MaterialLib::Fluid::FluidPropertyType::HeatCapacity, vars);

            auto const thermal_dispersivity_longitudinal =
                process_data.thermal_dispersivity_longitudinal(t, pos)[0];
            auto const thermal_dispersivity_transversal =
                process_data.thermal_dispersivity_transversal(t, pos)[0];

            // Use the fluid density model to compute the density
            auto const density = process_data.fluid_properties->getValue(
                MaterialLib::Fluid::FluidPropertyType::Density, vars);

            // Use the viscosity model to compute the viscosity
            auto const viscosity = process_data.fluid_properties->getValue(
                MaterialLib::Fluid::FluidPropertyType::Viscosity, vars);
            GlobalDimMatrixType K_over_mu = intrinsic_permeability / viscosity;

            GlobalDimVectorType const velocity =
                process_data.has_gravity
                    ? GlobalDimVectorType(-K_over_mu *
                                          (dNdx * p_nodal_values - density * b))
                    : GlobalDimVectorType(-K_over_mu * dNdx * p_nodal_values);

            double const velocity_magnitude = velocity.norm();
            GlobalDimMatrixType const thermal_dispersivity =
                fluid_reference_density * specific_heat_capacity_fluid *
                (thermal_dispersivity_transversal * velocity_magnitude * I +
                 (thermal_dispersivity_longitudinal -
                  thermal_dispersivity_transversal) /
                     velocity_magnitude * velocity * velocity.transpose());

            GlobalDimMatrixType const hydrodynamic_thermodispersion =
                thermal_conductivity * I + thermal_dispersivity;

            double const heat_capacity =
                density_solid * specific_heat_capacity_solid * (1 - porosity) +
                fluid_reference_density * specific_heat_capacity_fluid *
                    porosity;

            // matrix assembly
            Ktt.noalias() +=
                (dNdx.transpose() * hydrodynamic_thermodispersion * dNdx +
                 N.transpose() * velocity.transpose() * dNdx *
                     fluid_reference_density * specific_heat_capacity_fluid) *
                w;
            Kpp.noalias() += w * dNdx.transpose() * K_over_mu * dNdx;
            Mtt.noalias() += w * N.transpose() * heat_capacity * N;
            Mpp.noalias() += w * N.transpose() * specific_storage * N;
            if (process_data.has_gravity)
                Bp += w * density * dNdx.transpose() * K_over_mu * b;
            /* with Oberbeck-Boussing assumption density difference only exists
             * in buoyancy effects */
        }
    }
};

}  // namespace HT
}  // namespace ProcessLib
