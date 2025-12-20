// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <unordered_map>

#include "MaterialLib/MPL/MaterialSpatialDistributionMap.h"
#include "MeshLib/PropertyVector.h"
#include "ProcessLib/BoundaryConditionAndSourceTerm/Python/BHEInflowPythonBoundaryConditionPythonSideInterface.h"
#include "ProcessLib/HeatTransportBHE/BHE/BHETypes.h"
namespace MeshLib
{
class Element;
}

namespace ProcessLib::HeatTransportBHE
{
struct AlgebraicBCSetting
{
    const bool _use_algebraic_bc;

    const double _weighting_factor;
};

struct HeatTransportBHEProcessData final
{
    HeatTransportBHEProcessData(
        MaterialPropertyLib::MaterialSpatialDistributionMap media_map_,
        std::vector<BHE::BHETypes>&& vec_BHEs_,
        BHEInflowPythonBoundaryConditionPythonSideInterface* py_bc_object_ =
            nullptr,
        const bool use_tespy = false,
        const bool use_server_communication = false,
        const bool mass_lumping = false,
        AlgebraicBCSetting algebraicBCSetting = {false, 100.0},
        const bool is_linear = false)
        : media_map(media_map_),
          _vec_BHE_property(std::move(vec_BHEs_)),
          py_bc_object(py_bc_object_),
          _use_tespy(use_tespy),
          _use_server_communication(use_server_communication),
          _mass_lumping(mass_lumping),
          _algebraic_BC_Setting(algebraicBCSetting),
          _is_linear(is_linear)
    {
    }
    MaterialPropertyLib::MaterialSpatialDistributionMap media_map;

    std::vector<BHE::BHETypes> _vec_BHE_property;

    MeshLib::PropertyVector<int> const* _mesh_prop_materialIDs = nullptr;
    std::unordered_map<int, int> _map_materialID_to_BHE_ID{};

    //! Python object computing BC values.
    BHEInflowPythonBoundaryConditionPythonSideInterface* py_bc_object;

    const bool _use_tespy;

    const bool _use_server_communication;

    const bool _mass_lumping;

    std::vector<bool> mass_lumping_soil_elements;

    AlgebraicBCSetting const _algebraic_BC_Setting;

    const bool _is_linear;
};
}  // namespace ProcessLib::HeatTransportBHE
