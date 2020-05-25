/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
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
    : process_var_(process_var), param_name_(param_name){};

void FemCondition::setMesh(std::string const& mesh_name)
{
    base_type_ = BaseObjType::MESH;
    base_obj_name_ = mesh_name;
    obj_name_ = "";
}

void FemCondition::setGeoObject(std::string const& geo_name,
                                std::string const& obj_name)
{
    base_type_ = BaseObjType::GEOMETRY;
    base_obj_name_ = geo_name;
    obj_name_ = obj_name;
}

}  // namespace DataHolderLib
