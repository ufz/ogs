/**
 * \file
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
struct ProcessVariable
{
    std::string name;
    std::size_t order;
    std::size_t components;
};

enum class BaseObjType
{
    MESH = 0,
    GEOMETRY = 1
};

enum class ParameterType
{
    NONE = 0,
    CONSTANT,
    FUNCTION
};

/// Base class for boundary conditions, initial conditions and source terms.
class FemCondition
{
public:
    FemCondition(ProcessVariable const& process_var,
                 std::string const& param_name);
    virtual ~FemCondition() = default;

    /// Returns the name of the associated process variable
    std::string const getProcessVarName() const { return _process_var.name; }

    /// Returns the numerical order of the process vairable
    ProcessVariable const& getProcessVar() const { return _process_var; }

    /// Returns the name of the parameter associated with the condition
    std::string const getParamName() const { return _param_name; }

    /// Specifies if the condition is set a geometry or on a mesh
    BaseObjType getBaseObjType() const { return _base_type; }

    /// Returns the name of the base object (i.e. geometry or mesh)
    std::string getBaseObjName() const { return _base_obj_name; }

    /// Returns the name of the geometric object
    std::string const getObjName() const { return _obj_name; }

    /// Sets a mesh as corresponding object for the condition
    void setMesh(std::string const mesh_name);

    /// Sets a geometric object as corresponding object for the condition
    virtual void setGeoObject(std::string geo_name, std::string obj_name);

    /// Returns the type of condition for displaying purposes
    virtual std::string const getConditionClassStr() const = 0;

private:
    BaseObjType _base_type;
    ProcessVariable _process_var;
    std::string _param_name;
    std::string _base_obj_name;
    std::string _obj_name;
};

}  // namespace DataHolderLib
