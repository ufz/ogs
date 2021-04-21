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

#include <memory>

#include "ChemistryLib/ChemicalSolverInterface.h"
#include "MaterialLib/MPL/CreateMaterialSpatialDistributionMap.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"

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
    std::unique_ptr<MaterialPropertyLib::MaterialSpatialDistributionMap>
        media_map;
    Eigen::VectorXd const specific_body_force;
    bool const has_gravity;
    bool const non_advective_form;
    /// When this optional tag is on, the feedback of chemical reactions on the
    /// porosity will be counted. The change of porosity equals to the summation
    /// over the changes in the volume fractions of solid constituents. The
    /// change of the volume fraction, in terms of a solid constituent, results
    /// from chemical reactions.
    bool const chemically_induced_porosity_change;
    ChemistryLib::ChemicalSolverInterface* const chemical_solver_interface;

    const int hydraulic_process_id;
    // TODO (renchao-lu): This variable is used in the calculation of the
    // fluid's density and flux, indicating the transport process id. For now it
    // is assumed that these quantities depend on the first occurring transport
    // process only. The density and flux calculations have to be extended to
    // all processes.
    const int first_transport_process_id;

    MeshLib::PropertyVector<double>* mesh_prop_velocity = nullptr;
};

}  // namespace ComponentTransport
}  // namespace ProcessLib
