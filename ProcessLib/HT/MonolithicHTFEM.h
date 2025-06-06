/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  Created on October 11, 2017, 2:33 PM
 */

#pragma once

#include <Eigen/Core>
#include <typeinfo>
#include <vector>

#include "HTFEM.h"
#include "HTProcessData.h"
#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Utils/FormEigenTensor.h"
#include "NumLib/DOF/DOFTableUtil.h"
#include "NumLib/Extrapolation/ExtrapolatableElement.h"
#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/InitShapeMatrices.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "NumLib/NumericalStability/AdvectionMatrixAssembler.h"
#include "NumLib/NumericalStability/HydrodynamicDispersion.h"
#include "ParameterLib/Parameter.h"

namespace ProcessLib
{
namespace HT
{
const unsigned NUM_NODAL_DOF = 2;

template <typename ShapeFunction, int GlobalDim>
class MonolithicHTFEM : public HTFEM<ShapeFunction, GlobalDim>
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
    using NodalMatrixType = typename ShapeMatricesType::NodalMatrixType;
    using NodalRowVectorType = typename ShapeMatricesType::NodalRowVectorType;

    using GlobalDimVectorType = typename ShapeMatricesType::GlobalDimVectorType;
    using GlobalDimMatrixType = typename ShapeMatricesType::GlobalDimMatrixType;

public:
    MonolithicHTFEM(MeshLib::Element const& element,
                    std::size_t const local_matrix_size,
                    NumLib::GenericIntegrationMethod const& integration_method,
                    bool is_axially_symmetric,
                    HTProcessData const& process_data)
        : HTFEM<ShapeFunction, GlobalDim>(
              element, local_matrix_size, integration_method,
              is_axially_symmetric, process_data, NUM_NODAL_DOF)
    {
    }

    void assemble(double const t, double const dt,
                  std::vector<double> const& local_x,
                  std::vector<double> const& /*local_x_prev*/,
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

        auto KTT = local_K.template block<temperature_size, temperature_size>(
            temperature_index, temperature_index);
        auto MTT = local_M.template block<temperature_size, temperature_size>(
            temperature_index, temperature_index);
        auto Kpp = local_K.template block<pressure_size, pressure_size>(
            pressure_index, pressure_index);
        auto Mpp = local_M.template block<pressure_size, pressure_size>(
            pressure_index, pressure_index);
        auto Bp = local_b.template block<pressure_size, 1>(pressure_index, 0);

        auto const& process_data = this->_process_data;

        auto p_nodal_values = Eigen::Map<const NodalVectorType>(
            &local_x[pressure_index], pressure_size);

        auto const& medium =
            *process_data.media_map.getMedium(this->_element.getID());
        auto const& liquid_phase = medium.phase("AqueousLiquid");
        auto const& solid_phase = medium.phase("Solid");

        auto const& b =
            process_data
                .projected_specific_body_force_vectors[this->_element.getID()];

        MaterialPropertyLib::VariableArray vars;

        unsigned const n_integration_points =
            this->_integration_method.getNumberOfPoints();

        std::vector<GlobalDimVectorType> ip_flux_vector;
        double average_velocity_norm = 0.0;
        ip_flux_vector.reserve(n_integration_points);

        auto const& Ns =
            process_data.shape_matrix_cache
                .template NsHigherOrder<typename ShapeFunction::MeshElement>();

        for (unsigned ip(0); ip < n_integration_points; ip++)
        {
            auto const& ip_data = this->_ip_data[ip];
            auto const& dNdx = ip_data.dNdx;
            auto const& N = Ns[ip];
            auto const& w = ip_data.integration_weight;

            ParameterLib::SpatialPosition const pos{
                std::nullopt, this->_element.getID(),
                MathLib::Point3d(
                    NumLib::interpolateCoordinates<ShapeFunction,
                                                   ShapeMatricesType>(
                        this->_element, N))};

            double T_int_pt = 0.0;
            double p_int_pt = 0.0;
            // Order matters: First T, then P!
            NumLib::shapeFunctionInterpolate(local_x, N, T_int_pt, p_int_pt);

            vars.temperature = T_int_pt;
            vars.liquid_phase_pressure = p_int_pt;

            vars.liquid_saturation = 1.0;
            // \todo the argument to getValue() has to be changed for non
            // constant storage model
            auto const specific_storage =
                solid_phase.property(MaterialPropertyLib::PropertyType::storage)
                    .template value<double>(vars, pos, t, dt);

            auto const porosity =
                medium.property(MaterialPropertyLib::PropertyType::porosity)
                    .template value<double>(vars, pos, t, dt);
            vars.porosity = porosity;

            auto const intrinsic_permeability =
                MaterialPropertyLib::formEigenTensor<GlobalDim>(
                    medium
                        .property(
                            MaterialPropertyLib::PropertyType::permeability)
                        .value(vars, pos, t, dt));

            auto const specific_heat_capacity_fluid =
                liquid_phase
                    .property(MaterialPropertyLib::specific_heat_capacity)
                    .template value<double>(vars, pos, t, dt);

            // Use the fluid density model to compute the density
            auto const fluid_density =
                liquid_phase
                    .property(MaterialPropertyLib::PropertyType::density)
                    .template value<double>(vars, pos, t, dt);

            vars.density = fluid_density;
            // Use the viscosity model to compute the viscosity
            auto const viscosity =
                liquid_phase
                    .property(MaterialPropertyLib::PropertyType::viscosity)
                    .template value<double>(vars, pos, t, dt);
            GlobalDimMatrixType K_over_mu = intrinsic_permeability / viscosity;

            GlobalDimVectorType const velocity =
                process_data.has_gravity
                    ? GlobalDimVectorType(-K_over_mu * (dNdx * p_nodal_values -
                                                        fluid_density * b))
                    : GlobalDimVectorType(-K_over_mu * dNdx * p_nodal_values);

            // matrix assembly
            GlobalDimMatrixType const thermal_conductivity_dispersivity =
                this->getThermalConductivityDispersivity(
                    vars, fluid_density, specific_heat_capacity_fluid, velocity,
                    pos, t, dt);

            KTT.noalias() +=
                dNdx.transpose() * thermal_conductivity_dispersivity * dNdx * w;

            ip_flux_vector.emplace_back(velocity * fluid_density *
                                        specific_heat_capacity_fluid);
            average_velocity_norm += velocity.norm();

            Kpp.noalias() += w * dNdx.transpose() * K_over_mu * dNdx;
            MTT.noalias() += w *
                             this->getHeatEnergyCoefficient(
                                 vars, porosity, fluid_density,
                                 specific_heat_capacity_fluid, pos, t, dt) *
                             N.transpose() * N;
            Mpp.noalias() += w * N.transpose() * specific_storage * N;
            if (process_data.has_gravity)
            {
                Bp += w * fluid_density * dNdx.transpose() * K_over_mu * b;
            }
            /* with Oberbeck-Boussing assumption density difference only exists
             * in buoyancy effects */
        }

        NumLib::assembleAdvectionMatrix<typename ShapeFunction::MeshElement>(
            process_data.stabilizer, this->_ip_data,
            process_data.shape_matrix_cache, ip_flux_vector,
            average_velocity_norm / static_cast<double>(n_integration_points),
            KTT);
    }

    std::vector<double> const& getIntPtDarcyVelocity(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const override
    {
        int const process_id = 0;  // monolithic case.
        auto const indices =
            NumLib::getIndices(this->_element.getID(), *dof_table[process_id]);
        assert(!indices.empty());
        auto const& local_x = x[process_id]->get(indices);

        return this->getIntPtDarcyVelocityLocal(t, local_x, cache);
    }

private:
    using HTFEM<ShapeFunction, GlobalDim>::pressure_index;
    using HTFEM<ShapeFunction, GlobalDim>::pressure_size;
    using HTFEM<ShapeFunction, GlobalDim>::temperature_index;
    using HTFEM<ShapeFunction, GlobalDim>::temperature_size;
};

}  // namespace HT
}  // namespace ProcessLib
