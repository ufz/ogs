/**
 * \file
 * \copyright
 * Copyright (c) 2012-2023, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <Eigen/Core>
#include <memory>

#include "ChemistryLib/ChemicalSolverInterface.h"
#include "LookupTable.h"
#include "MaterialLib/MPL/MaterialSpatialDistributionMap.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "NumLib/NumericalStability/NumericalStabilization.h"
#include "ParameterLib/ConstantParameter.h"
#include "ParameterLib/Parameter.h"

namespace MaterialPropertyLib
{
class Medium;
}

namespace MeshLib
{
template <typename PROP_VAL_TYPE>
class PropertyVector;
}

namespace ProcessLib
{
namespace ComponentTransport
{
struct ComponentTransportProcessData
{
    MaterialPropertyLib::MaterialSpatialDistributionMap media_map;
    bool const has_gravity;
    bool const non_advective_form;
    /// This optional tag provides a simple means of considering the temperature
    /// effect on the solute transport process.
    ParameterLib::Parameter<double> const* const temperature;
    /**
     * When this optional tag is on, the feedback of chemical reactions on the
     * porosity will be counted. The change of porosity equals to the summation
     * over the changes in the volume fractions of solid constituents. The
     * change of the volume fraction, in terms of a solid constituent, results
     * from chemical reactions.
     *
     * \note In order to use this optional tag, the amount of solid constituents
     * should be given as volume fraction instead of molality. In addition,
     * an appropriate molar volume is required for each solid. The relationship
     * to calculate volume fractions of m solids from molalities is as follows:
     * \f[
     * b_i = \frac{n_i}{m^l} = \frac{\phi_i}{\rho^l \phi V_{m,i}}, i=1,...,m
     * \f]
     * where \f$b_i\f$ is the molality in mol/kg of water,
     * \f$n_i\f$ is the amount of solid in mol,
     * \f$m^l\f$ is the mass of water in kg,
     * \f$\phi_i\f$ is the volume fraction of solid i,
     * \f$\rho^l\f$ is the density of water in kg/m\f$^3\f$,
     * \f$\phi\f$ is the porosity,
     * \f$V_{m,i}\f$ is the molar volume of solid i in m\f$^3\f$/mol.
     */
    bool const chemically_induced_porosity_change;
    ChemistryLib::ChemicalSolverInterface* const chemical_solver_interface;
    std::unique_ptr<LookupTable> lookup_table;

    NumLib::NumericalStabilization stabilizer;

    /// Projected specific body force vector: R * R^T * b.
    std::vector<Eigen::VectorXd> const projected_specific_body_force_vectors;

    int const mesh_space_dimension;

    /// The aperture size is the thickness of 2D element or the
    /// cross section area of 1D element. For 3D element, the value is set to 1.
    ParameterLib::Parameter<double> const& aperture_size =
        ParameterLib::ConstantParameter<double>("constant_one", 1.0);

    static const int hydraulic_process_id = 0;
    // TODO (renchao-lu): This variable is used in the calculation of the
    // fluid's density and flux, indicating the transport process id. For now it
    // is assumed that these quantities depend on the first occurring transport
    // process only. The density and flux calculations have to be extended to
    // all processes.
    static const int first_transport_process_id = 1;

    MeshLib::PropertyVector<double>* mesh_prop_velocity = nullptr;
    MeshLib::PropertyVector<double>* mesh_prop_porosity = nullptr;
};

}  // namespace ComponentTransport
}  // namespace ProcessLib
