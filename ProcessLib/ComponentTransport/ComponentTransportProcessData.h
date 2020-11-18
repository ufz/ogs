/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>

#include "ChemicalProcessData.h"
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
    std::unique_ptr<ChemicalProcessData> chemical_process_data;

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
