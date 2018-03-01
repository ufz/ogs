/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 */

#pragma once

#include <string>

namespace DataHolderLib
{

enum class BaseObjType
{
    MESH = 0,
    GEOMETRY = 1
};

enum class ConditionType
{
    NONE = 0,
    DIRICHLET,
    NONLINEARDIRICHLET,
    NEUMANN,
    NONLINEARNEUMANN,
    ROBIN
};

enum class ParameterType
{
    NONE,
    CONSTANT,
    FUNCTION
};

class Condition
{
public:
    Condition(std::string const process_var, std::string const param_name, ConditionType type);

    Condition(Condition const& c);

    ~Condition() {};

    /// Returns the name of the associated process variable
    std::string const getProcessVarName() const { return _process_var; }

    std::string const getParamName() const { return _param_name; }

    /// Is the condition set a geometry or on a mesh
    BaseObjType getBaseObjType() const { return _base_type; }

    ///Returns the name of the base object (i.e. geometry or mesh)
    std::string getBaseObjName() const { return _base_obj_name; }

    ///Returns the name of the geometric object
    std::string const getObjName() const { return _obj_name; }

    ConditionType getType() const { return _type; }

    void setMesh(std::string const mesh_name);

    virtual void setGeoObject(std::string geo_name, std::string obj_name);

    static std::string convertTypeToString(ConditionType type);

    static ConditionType convertStringToType(std::string const& type_str);

private:
    BaseObjType _base_type;
    ConditionType _type;
    std::string _process_var;
    std::string _param_name;
    std::string _base_obj_name;
    std::string _obj_name;
};

} // namespace
