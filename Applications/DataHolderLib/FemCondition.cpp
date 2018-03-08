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

FemCondition::FemCondition(std::string const process_var, std::string const param_name, ConditionType type)
: _type(type), _process_var(process_var), _param_name(param_name)
{};

FemCondition::FemCondition(FemCondition const& c)
: _base_type(c.getBaseObjType()), _type(c.getType()),
  _process_var(c.getProcessVarName()), _param_name(c.getParamName()),
  _base_obj_name(c.getBaseObjName()), _obj_name(c.getObjName())
{}

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

ConditionType FemCondition::convertStringToType(std::string const& str)
{
    if ( str == "Dirichlet")
        return ConditionType::DIRICHLET;
    else if (str == "NonlinearDirichlet")
        return ConditionType::NONLINEARDIRICHLET;
    else if (str == "Neumann")
        return ConditionType::NEUMANN;
    else if (str == "NonuniformNeumann")
        return ConditionType::NONLINEARNEUMANN;
    else if (str == "Robin")
        return ConditionType::ROBIN;

    return ConditionType::NONE;
}

std::string FemCondition::convertTypeToString(ConditionType type)
{
    if (type == ConditionType::DIRICHLET)
        return "Dirichlet";
    else if (type == ConditionType::NONLINEARDIRICHLET)
        return "NonlinearDirichlet";
    else if (type == ConditionType::NEUMANN)
        return "Neumann";
    else if (type == ConditionType::NONLINEARDIRICHLET)
        return "NonlinearNeumann";
    else if (type == ConditionType::ROBIN)
        return "Robin";

    return "";
}


} // namespace
