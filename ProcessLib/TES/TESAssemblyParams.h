/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <Eigen/Eigen>
#include <Eigen/Sparse>

#include "MaterialLib/PhysicalConstant.h"
#include "MaterialLib/Adsorption/Reaction.h"

#include "ProcessLib/VariableTransformation.h"

namespace ProcessLib
{
namespace TES
{
const unsigned NODAL_DOF = 3;
const unsigned COMPONENT_ID_PRESSURE = 0;
const unsigned COMPONENT_ID_TEMPERATURE = 1;
const unsigned COMPONENT_ID_MASS_FRACTION = 2;

const double M_N2 = 0.028013;
const double M_H2O = 0.018016;

struct AssemblyParams
{
    Trafo trafo_p{1.0};
    Trafo trafo_T{1.0};
    Trafo trafo_x{1.0};

    std::unique_ptr<Adsorption::Reaction> react_sys;

    double fluid_specific_heat_source =
        std::numeric_limits<double>::quiet_NaN();
    double cpG = std::numeric_limits<
        double>::quiet_NaN();  // specific isobaric fluid heat capacity

    Eigen::MatrixXd solid_perm_tensor = Eigen::MatrixXd::Constant(
        3, 3, std::numeric_limits<double>::quiet_NaN());  // TODO get dimensions
    double solid_specific_heat_source =
        std::numeric_limits<double>::quiet_NaN();
    double solid_heat_cond = std::numeric_limits<double>::quiet_NaN();
    double cpS = std::numeric_limits<
        double>::quiet_NaN();  // specific isobaric solid heat capacity

    double tortuosity = std::numeric_limits<double>::quiet_NaN();
    double diffusion_coefficient_component =
        std::numeric_limits<double>::quiet_NaN();

    double poro = std::numeric_limits<double>::quiet_NaN();

    double rho_SR_dry = std::numeric_limits<double>::quiet_NaN();

    const double M_inert = MaterialLib::PhysicalConstant::MolarMass::N2;
    const double M_react = MaterialLib::PhysicalConstant::MolarMass::Water;

    // TODO unify variable names
    double initial_solid_density = std::numeric_limits<double>::quiet_NaN();

    double delta_t = std::numeric_limits<double>::quiet_NaN();
    unsigned iteration_in_current_timestep = 0;

    bool output_element_matrices = false;

    unsigned number_of_try_of_iteration = 0;
    double current_time = std::numeric_limits<double>::quiet_NaN();

    //! Output global matrix/rhs after first iteration.
    std::size_t timestep = 0;
    std::size_t total_iteration = 0;
};

}  // namespace TES

}  // namespace ProcessLib
