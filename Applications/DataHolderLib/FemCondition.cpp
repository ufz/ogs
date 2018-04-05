/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 */

#include "FemCondition.h"

namespace DataHolderLib
{
FemCondition::FemCondition(ProcessVariable const& process_var,
                           std::string const& param_name)
    : _process_var(process_var), _param_name(param_name){};

void FemCondition::setMesh(std::string const mesh_name)
{
    _base_type = BaseObjType::MESH;
    _base_obj_name = mesh_name;
    _obj_name = "";
}

void FemCondition::setGeoObject(std::string geo_name, std::string obj_name)
{
    _base_type = BaseObjType::GEOMETRY;
    _base_obj_name = geo_name;
    _obj_name = obj_name;
}

}  // namespace
