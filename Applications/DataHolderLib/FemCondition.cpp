// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "FemCondition.h"

namespace DataHolderLib
{
FemCondition::FemCondition(ProcessVariable const& process_var,
                           std::string const& param_name)
    : _process_var(process_var), _param_name(param_name){};

void FemCondition::setMesh(std::string const& mesh_name)
{
    _base_type = BaseObjType::MESH;
    _base_obj_name = mesh_name;
    _obj_name = "";
}

void FemCondition::setGeoObject(std::string const& geo_name,
                                std::string const& obj_name)
{
    _base_type = BaseObjType::GEOMETRY;
    _base_obj_name = geo_name;
    _obj_name = obj_name;
}

}  // namespace DataHolderLib
