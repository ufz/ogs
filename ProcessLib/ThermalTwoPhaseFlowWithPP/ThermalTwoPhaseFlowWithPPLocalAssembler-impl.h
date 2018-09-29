/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/**
* Common convenitions for naming:
* X_gas_nonwet           mass fraction of gas component(e.g air) in nonwetting phase
* (gas component doesn't include water vapor, same for the following)
* x_gas_nonwet           molar fraction of gas component in nonwetting phase
* x_vapor_nonwet         molar fraction of vapor in nonwetting phase
* p_vapor_nonwet         water vapor pressure
* p_gas_nonwet           partial pressure of gas component
* mol_density_nonwet         molar density of nonwetting phase
* mol_density_water          molar density of water
* density_water              mass density of water
* density_nonwet_gas         mass density of gas component in the nonwetting phase
* density_nonwet_vapor       mass density of vapor in the nonwetting phase
* density_nonwet             mass density of the nonwetting phase
* density_wet                mass density of wetting pahse
* density_solid              mass density of the solid phase
* velocity_nonwet              velocity of nonwetting phase
* velocity_wet                 velocity of wetting phase
* heat_capacity_dry_gas        heat capacity of dry gas
* heat_capacity_water_vapor    heat capacity of water vapor
* heat_capacity_water          heat capacity of liquid water
* heat_capacity_solid          heat capacity of soil grain
* latent_heat_evaporation      latent heat for evaporation(water to vapor)
* enthalpy_nonwet_gas          enthalpy of gas component in the nonwetting phase
* enthalpy_nonwet_vapor        enthalpy of water vapor in the nonwetting phase
* enthalpy_wet                 enthalpy of wetting phase
* enthalpy_nonwet                 enthalpy of the nonwetting phase
* internal_energy_nonwet        specific internal energy for the nonwetting phase
* internal_energy_wet           specific internal energy for the wetting phase
* heat_conductivity_dry_solid   heat conductivity of the dry porous medium
* heat_conductivity_wet_solid   heat conductivity of the fully saturated porous medium
* heat_conductivity_unsaturated   heat conductivity of the unsaturated porous medium
*/
#pragma once

#include "ThermalTwoPhaseFlowWithPPLocalAssembler.h"

#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
#include "NumLib/Function/Interpolation.h"

#include "ThermalTwoPhaseFlowWithPPProcessData.h"

namespace ProcessLib
{
namespace ThermalTwoPhaseFlowWithPP
{
template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
void ThermalTwoPhaseFlowWithPPLocalAssembler<
    ShapeFunction, IntegrationMethod,
    GlobalDim>::assemble(double const t, std::vector<double> const& local_x,
                         std::vector<double>& local_M_data,
                         std::vector<double>& local_K_data,
                         std::vector<double>& local_b_data)
{
    using MaterialLib::PhysicalConstant::IdealGasConstant;
    auto const& water_mol_mass =
        MaterialLib::PhysicalConstant::MolarMass::Water;
    auto const& air_mol_mass = MaterialLib::PhysicalConstant::MolarMass::Air;

    auto const local_matrix_size = local_x.size();

    assert(local_matrix_size == ShapeFunction::NPOINTS * NUM_NODAL_DOF);

    auto local_M = MathLib::createZeroedMatrix<LocalMatrixType>(
        local_M_data, local_matrix_size, local_matrix_size);
    auto local_K = MathLib::createZeroedMatrix<LocalMatrixType>(
        local_K_data, local_matrix_size, local_matrix_size);
    auto local_b = MathLib::createZeroedVector<LocalVectorType>(
        local_b_data, local_matrix_size);

    auto Mgp =
        local_M.template block<nonwet_pressure_size, nonwet_pressure_size>(
            nonwet_pressure_matrix_index, nonwet_pressure_matrix_index);
    auto Mgpc = local_M.template block<nonwet_pressure_size, cap_pressure_size>(
        nonwet_pressure_matrix_index, cap_pressure_matrix_index);
    auto Mgt = local_M.template block<nonwet_pressure_size, temperature_size>(
        nonwet_pressure_matrix_index, temperature_matrix_index);

    auto Mlp = local_M.template block<cap_pressure_size, nonwet_pressure_size>(
        cap_pressure_matrix_index, nonwet_pressure_matrix_index);
    auto Mlpc = local_M.template block<cap_pressure_size, cap_pressure_size>(
        cap_pressure_matrix_index, cap_pressure_matrix_index);
    auto Mlt = local_M.template block<cap_pressure_size, temperature_size>(
        cap_pressure_matrix_index, temperature_matrix_index);

    auto Mep = local_M.template block<temperature_size, nonwet_pressure_size>(
        temperature_matrix_index, nonwet_pressure_matrix_index);
    auto Mepc = local_M.template block<temperature_size, cap_pressure_size>(
        temperature_matrix_index, cap_pressure_matrix_index);
    auto Met = local_M.template block<temperature_size, temperature_size>(
        temperature_matrix_index, temperature_matrix_index);

    NodalMatrixType laplace_operator =
        NodalMatrixType::Zero(ShapeFunction::NPOINTS, ShapeFunction::NPOINTS);

    auto Kgp =
        local_K.template block<nonwet_pressure_size, nonwet_pressure_size>(
            nonwet_pressure_matrix_index, nonwet_pressure_matrix_index);
    auto Kgpc = local_K.template block<nonwet_pressure_size, cap_pressure_size>(
        nonwet_pressure_matrix_index, cap_pressure_matrix_index);
    auto Kgt = local_K.template block<nonwet_pressure_size, temperature_size>(
        nonwet_pressure_matrix_index, temperature_matrix_index);

    auto Klp = local_K.template block<cap_pressure_size, nonwet_pressure_size>(
        cap_pressure_matrix_index, nonwet_pressure_matrix_index);
    auto Klpc = local_K.template block<cap_pressure_size, cap_pressure_size>(
        cap_pressure_matrix_index, cap_pressure_matrix_index);
    auto Klt = local_K.template block<cap_pressure_size, temperature_size>(
        cap_pressure_matrix_index, temperature_matrix_index);

    auto Kep = local_K.template block<temperature_size, nonwet_pressure_size>(
        temperature_matrix_index, nonwet_pressure_matrix_index);
    auto Kepc = local_K.template block<temperature_size, cap_pressure_size>(
        temperature_matrix_index, cap_pressure_matrix_index);
    auto Ket = local_K.template block<temperature_size, temperature_size>(
        temperature_matrix_index, temperature_matrix_index);

    auto Bg = local_b.template segment<nonwet_pressure_size>(
        nonwet_pressure_matrix_index);
    auto Bl =
        local_b.template segment<cap_pressure_size>(cap_pressure_matrix_index);
    auto Be =
        local_b.template segment<temperature_size>(temperature_matrix_index);

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    SpatialPosition pos;
    pos.setElementID(_element.getID());
    auto const& two_phase_material_model =
        _process_data.material->getTwoPhaseMaterialModel();
    const int material_id =
        two_phase_material_model.getMaterialID(pos.getElementID().get());

    auto const num_nodes = ShapeFunction::NPOINTS;
    auto const pg_nodal_values =
        Eigen::Map<const NodalVectorType>(&local_x[0], num_nodes);
    auto const pc_nodal_values =
        Eigen::Map<const NodalVectorType>(&local_x[num_nodes], num_nodes);

    const Eigen::MatrixXd& perm = two_phase_material_model.getPermeability(
        material_id, t, pos, _element.getDimension());
    assert(perm.rows() == _element.getDimension() || perm.rows() == 1);
    GlobalDimMatrixType permeability = GlobalDimMatrixType::Zero(
        _element.getDimension(), _element.getDimension());
    if (perm.rows() == _element.getDimension())
        permeability = perm;
    else if (perm.rows() == 1)
        permeability.diagonal().setConstant(perm(0, 0));

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        double pg_int_pt = 0.;
        double pc_int_pt = 0.;
        double T_int_pt = 0.0;
        NumLib::shapeFunctionInterpolate(local_x, _ip_data[ip].N, pg_int_pt,
                                         pc_int_pt, T_int_pt);

        double const density_water =
            two_phase_material_model.getLiquidDensity(pg_int_pt, T_int_pt);

        double const Sw = two_phase_material_model.getSaturation(
            material_id, t, pos, pg_int_pt, T_int_pt, pc_int_pt);

        _saturation[ip] = Sw;

        double dSwdpc =
            (pc_int_pt > two_phase_material_model.getCapillaryPressure(
                             material_id, t, pos, pg_int_pt, T_int_pt, 0.0))
                ? 0.0
                : two_phase_material_model.getSaturationDerivative(
                      material_id, t, pos, pg_int_pt, T_int_pt, Sw);

        double const p_vapor_nonwet =
            _process_data.material->calculateVaporPressureNonwet(
                pc_int_pt, T_int_pt, density_water);
        // partial pressure of gas component
        double const p_gas_nonwet = pg_int_pt - p_vapor_nonwet;
        // molar fraction of gas component in nonwet phase
        double const x_gas_nonwet = p_gas_nonwet / pg_int_pt;
        // molar fraction of water vapor in nonwet phase
        double const x_vapor_nonwet = p_vapor_nonwet / pg_int_pt;
        // mass fraction of gas component in the nonwet phase
        double const X_gas_nonwet =
            x_gas_nonwet /
            (x_gas_nonwet + x_vapor_nonwet * water_mol_mass / air_mol_mass);
        double const mol_density_nonwet = pg_int_pt / IdealGasConstant / T_int_pt;
        double const mol_density_water = density_water / water_mol_mass;

        double const d_mol_density_nonwet_d_pg = 1 / IdealGasConstant / T_int_pt;
        double const d_p_vapor_nonwet_d_T =
            _process_data.material->calculateDerivativedPgwdT(
                pc_int_pt, T_int_pt, density_water);
        double const d_p_vapor_nonwet_d_pc =
            _process_data.material->calculateDerivativedPgwdPC(
                pc_int_pt, T_int_pt, density_water);
        double const d_mol_density_nonwet_d_T =
            -pg_int_pt / IdealGasConstant / T_int_pt / T_int_pt;
        double const d_x_gas_nonwet_d_pg =
            p_vapor_nonwet / pg_int_pt / pg_int_pt;
        double const d_x_gas_nonwet_d_pc = -d_p_vapor_nonwet_d_pc / pg_int_pt;
        double const d_x_gas_nonwet_d_T = -d_p_vapor_nonwet_d_T / pg_int_pt;

        double const density_nonwet_gas =
            p_gas_nonwet * air_mol_mass / IdealGasConstant / T_int_pt;
        double const density_nonwet_vapor =
            p_vapor_nonwet * water_mol_mass / IdealGasConstant / T_int_pt;
        double const density_nonwet = density_nonwet_gas + density_nonwet_vapor;
        double const density_wet = density_water;
        double const density_solid = _process_data.density_solid(t, pos)[0];
        // Derivative of nonwet phase density in terms of T
        double const d_density_nonwet_d_T =
            _process_data.material->calculatedDensityNonwetdT (
                p_gas_nonwet, p_vapor_nonwet, pc_int_pt, T_int_pt, density_water);

        _pressure_wetting[ip] = pg_int_pt - pc_int_pt;
        // heat capacity of nonwet phase
        double const heat_capacity_dry_gas =
            _process_data.material->getSpecificHeatCapacityAir(pg_int_pt,
                                                               T_int_pt);
        const double heat_capacity_water_vapor =
            _process_data.material->getSpecificHeatCapacityVapor(pg_int_pt,
                                                                 T_int_pt);

        double const heat_capacity_water =
            _process_data.material->getSpecificHeatCapacityWater(pg_int_pt,
                                                                 T_int_pt);
        double const heat_capacity_solid =
            _process_data.material->getSpecificHeatCapacitySolid(pg_int_pt,
                                                                 T_int_pt);
        double const latent_heat_evaporation =
            _process_data.latent_heat_evaporation(t, pos)[0];

        double const enthalpy_nonwet_gas =
            _process_data.material->getAirEnthalpySimple(
                T_int_pt, heat_capacity_dry_gas, pg_int_pt);

        double const enthalpy_wet =
            _process_data.material->getLiquidWaterEnthalpySimple(
                T_int_pt, heat_capacity_water, _pressure_wetting[ip]);

        double const enthalpy_nonwet_vapor =
            _process_data.material->getWaterVaporEnthalpySimple(
                T_int_pt, heat_capacity_water_vapor, pg_int_pt,
                latent_heat_evaporation);
        double const enthalpy_nonwet =
            enthalpy_nonwet_gas * X_gas_nonwet +
            enthalpy_nonwet_vapor * (1 - X_gas_nonwet);
        double const internal_energy_nonwet =
            enthalpy_nonwet - pg_int_pt / density_nonwet;
        double const internal_energy_wet = enthalpy_wet;
        /// Derivative
        double const d_enthalpy_gas_nonwet_d_T =
            heat_capacity_dry_gas + IdealGasConstant / air_mol_mass;
        double const d_enthalpy_nonwet_d_T =
            heat_capacity_water * (1 - X_gas_nonwet) +
            d_enthalpy_gas_nonwet_d_T * X_gas_nonwet;
        // Assemble M matrix
        // nonwetting
        double const porosity = two_phase_material_model.getPorosity(
            material_id, t, pos, pg_int_pt, T_int_pt, 0);

        Mgp.noalias() += porosity *
                         ((1 - Sw) * (mol_density_nonwet * d_x_gas_nonwet_d_pg +
                                      x_gas_nonwet * d_mol_density_nonwet_d_pg)) *
                         _ip_data[ip].mass_operator;
        Mgpc.noalias() += porosity *
                          ((1 - Sw) * mol_density_nonwet * d_x_gas_nonwet_d_pc -
                           mol_density_nonwet * x_gas_nonwet * dSwdpc) *
                          _ip_data[ip].mass_operator;
        Mgt.noalias() += porosity *
                         ((1 - Sw) * (mol_density_nonwet * d_x_gas_nonwet_d_T +
                                      x_gas_nonwet * d_mol_density_nonwet_d_T)) *
                         _ip_data[ip].mass_operator;

        Mlpc.noalias() +=
            porosity *
            ((1 - Sw) * d_p_vapor_nonwet_d_pc / IdealGasConstant / T_int_pt +
             mol_density_nonwet * x_vapor_nonwet * (-dSwdpc) +
             dSwdpc * mol_density_water) *
            _ip_data[ip].mass_operator;
        Mlt.noalias() +=
            porosity *
            ((1 - Sw) *
             (d_p_vapor_nonwet_d_T / IdealGasConstant / T_int_pt -
              p_vapor_nonwet / IdealGasConstant / T_int_pt / T_int_pt)) *
            _ip_data[ip].mass_operator;

        Mep.noalias() +=
            porosity *
            ((x_gas_nonwet * air_mol_mass + x_vapor_nonwet * water_mol_mass) *
                 d_mol_density_nonwet_d_pg * enthalpy_nonwet -
             mol_density_nonwet * (water_mol_mass - air_mol_mass) *
                 d_x_gas_nonwet_d_pg * enthalpy_nonwet -
             1) *
            (1 - Sw) * _ip_data[ip].mass_operator;
        Mepc.noalias() +=
            porosity * (density_wet * internal_energy_wet -
                        density_nonwet * internal_energy_nonwet) *
                dSwdpc * _ip_data[ip].mass_operator +
            porosity * ((water_mol_mass - air_mol_mass) * enthalpy_nonwet /
                        IdealGasConstant / T_int_pt) *
                (1 - Sw) * d_p_vapor_nonwet_d_pc * _ip_data[ip].mass_operator;
        Met.noalias() +=
            ((1 - porosity) * density_solid * heat_capacity_solid +
             porosity * ((1 - Sw) * (d_density_nonwet_d_T * enthalpy_nonwet +
                                     density_nonwet * d_enthalpy_nonwet_d_T) +
                         Sw * density_wet * heat_capacity_water)) *
            _ip_data[ip].mass_operator;

        // nonwet
        double const k_rel_nonwet =
            two_phase_material_model.getNonwetRelativePermeability(
                t, pos, _pressure_wetting[ip], T_int_pt, Sw);
        double const mu_nonwet = two_phase_material_model.getGasViscosity(
            _pressure_wetting[ip], T_int_pt);
        double const lambda_nonwet = k_rel_nonwet / mu_nonwet;
        double const diffusion_coeff_component_gas =
            _process_data.diffusion_coeff_component_b(t, pos)[0];

        // wet
        double const k_rel_wet =
            two_phase_material_model.getWetRelativePermeability(
                t, pos, pg_int_pt, T_int_pt, Sw);
        double const mu_wet =
            two_phase_material_model.getLiquidViscosity(pg_int_pt, T_int_pt);
        double const lambda_wet = k_rel_wet / mu_wet;

        GlobalDimVectorType const velocity_nonwet =
            -lambda_nonwet * permeability *
            (_ip_data[ip].dNdx * pg_nodal_values);
        GlobalDimVectorType const velocity_wet =
            -lambda_wet * permeability *
            (_ip_data[ip].dNdx * (pg_nodal_values - pc_nodal_values));

        laplace_operator.noalias() = _ip_data[ip].dNdx.transpose() *
                                     permeability * _ip_data[ip].dNdx *
                                     _ip_data[ip].integration_weight;

        Ket.noalias() +=
            _ip_data[ip].integration_weight * _ip_data[ip].N.transpose() *
                (d_density_nonwet_d_T * enthalpy_nonwet +
                 density_nonwet * d_enthalpy_nonwet_d_T) *
                velocity_nonwet.transpose() * _ip_data[ip].dNdx +
            _ip_data[ip].integration_weight * _ip_data[ip].N.transpose() *
                heat_capacity_water * density_water * velocity_wet.transpose() *
                _ip_data[ip].dNdx;

        double const heat_conductivity_dry_solid =
            _process_data.material->getThermalConductivityDrySolid(pg_int_pt,
                                                                   T_int_pt);
        double const heat_conductivity_wet_solid =
            _process_data.material->getThermalConductivityWetSolid(pg_int_pt,
                                                                   T_int_pt);
        double const heat_conductivity_unsaturated =
            _process_data.material->calculateUnsatHeatConductivity(
                t, pos, Sw, heat_conductivity_dry_solid,
                heat_conductivity_wet_solid);
        // Laplace
        Kgp.noalias() +=
            (mol_density_nonwet * x_gas_nonwet * lambda_nonwet) * laplace_operator +
            ((1 - Sw) * porosity * diffusion_coeff_component_gas *
             mol_density_nonwet * d_x_gas_nonwet_d_pg) *
                _ip_data[ip].diffusion_operator;
        Kgpc.noalias() += ((1 - Sw) * porosity * diffusion_coeff_component_gas *
                           mol_density_nonwet * d_x_gas_nonwet_d_pc) *
                          _ip_data[ip].diffusion_operator;
        Kgt.noalias() += ((1 - Sw) * porosity * diffusion_coeff_component_gas *
                          mol_density_nonwet * d_x_gas_nonwet_d_T) *
                         _ip_data[ip].diffusion_operator;

        Klp.noalias() += (mol_density_nonwet * x_vapor_nonwet * lambda_nonwet) *
                             laplace_operator +
                         mol_density_water * lambda_wet * laplace_operator -
                         ((1 - Sw) * porosity * diffusion_coeff_component_gas *
                          mol_density_nonwet * d_x_gas_nonwet_d_pg) *
                             _ip_data[ip].diffusion_operator;
        Klpc.noalias() += (-mol_density_water * lambda_wet * laplace_operator) -
                          ((1 - Sw) * porosity * diffusion_coeff_component_gas *
                           mol_density_nonwet * d_x_gas_nonwet_d_pc) *
                              _ip_data[ip].diffusion_operator;
        Klt.noalias() += -((1 - Sw) * porosity * diffusion_coeff_component_gas *
                           mol_density_nonwet * d_x_gas_nonwet_d_T) *
                         _ip_data[ip].diffusion_operator;

        Kep.noalias() +=
            (lambda_nonwet * density_nonwet * enthalpy_nonwet +
             lambda_wet * density_wet * enthalpy_wet) *
                laplace_operator +
            (1 - Sw) * porosity * diffusion_coeff_component_gas *
                mol_density_nonwet * (air_mol_mass * enthalpy_nonwet_gas -
                                  water_mol_mass * enthalpy_nonwet_vapor) *
                d_x_gas_nonwet_d_pg * _ip_data[ip].diffusion_operator;
        Kepc.noalias() +=
            -lambda_wet * enthalpy_wet * density_wet * laplace_operator +
            (1 - Sw) * porosity * diffusion_coeff_component_gas *
                mol_density_nonwet * (air_mol_mass * enthalpy_nonwet_gas -
                                  water_mol_mass * enthalpy_nonwet_vapor) *
                d_x_gas_nonwet_d_pc * _ip_data[ip].diffusion_operator;
        Ket.noalias() +=
            _ip_data[ip].dNdx.transpose() * heat_conductivity_unsaturated *
                _ip_data[ip].dNdx * _ip_data[ip].integration_weight +
            (1 - Sw) * porosity * diffusion_coeff_component_gas *
                mol_density_nonwet * (air_mol_mass * enthalpy_nonwet_gas -
                                  water_mol_mass * enthalpy_nonwet_vapor) *
                d_x_gas_nonwet_d_T * _ip_data[ip].diffusion_operator;

        if (_process_data.has_gravity)
        {
            auto const& b = _process_data.specific_body_force;
            NodalVectorType gravity_operator = _ip_data[ip].dNdx.transpose() *
                                               permeability * b *
                                               _ip_data[ip].integration_weight;
            Bg.noalias() +=
                (mol_density_nonwet * x_gas_nonwet * lambda_nonwet * density_nonwet) *
                gravity_operator;
            Bl.noalias() +=
                (mol_density_water * lambda_wet * density_wet +
                 mol_density_nonwet * x_vapor_nonwet * lambda_nonwet * density_nonwet) *
                gravity_operator;
            Be.noalias() +=
                (lambda_nonwet * density_nonwet * density_nonwet * enthalpy_nonwet +
                 lambda_wet * density_wet * density_wet * enthalpy_wet) *
                gravity_operator;
        }  // end of has gravity
    }
    if (_process_data.has_mass_lumping)
    {
        for (unsigned row = 0; row < Mgp.cols(); row++)
        {
            for (unsigned column = 0; column < Mgp.cols(); column++)
            {
                if (row != column)
                {
                    Mgp(row, row) += Mgp(row, column);
                    Mgp(row, column) = 0.0;
                    Mgpc(row, row) += Mgpc(row, column);
                    Mgpc(row, column) = 0.0;
                    Mgt(row, row) += Mgt(row, column);
                    Mgt(row, column) = 0.0;
                    Mlp(row, row) += Mlp(row, column);
                    Mlp(row, column) = 0.0;
                    Mlpc(row, row) += Mlpc(row, column);
                    Mlpc(row, column) = 0.0;
                    Mlt(row, row) += Mlt(row, column);
                    Mlt(row, column) = 0.0;
                    Mep(row, row) += Mep(row, column);
                    Mep(row, column) = 0.0;
                    Mepc(row, row) += Mepc(row, column);
                    Mepc(row, column) = 0.0;
                    Met(row, row) += Met(row, column);
                    Met(row, column) = 0.0;
                }
            }
        }
    }  // end of mass-lumping
}

}  // end of namespace
}  // end of namespace
