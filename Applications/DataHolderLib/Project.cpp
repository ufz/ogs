/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "Project.h"

#include <algorithm>

#include <logog/include/logog.hpp>

#include "BaseLib/FileTools.h"
#include "BaseLib/uniqueInsert.h"

#include "MeshLib/Mesh.h"

namespace DataHolderLib
{

void Project::addMesh(std::unique_ptr<MeshLib::Mesh> mesh)
{
    std::string name = mesh->getName();
    getUniqueName(name);
    mesh->setName(name);
    _mesh_vec.push_back(std::move(mesh));
}

std::vector<std::unique_ptr<MeshLib::Mesh>>::const_iterator
Project::findMeshByName(std::string const& name) const
{
    return const_cast<Project&>(*this).findMeshByName(name);
}

std::vector<std::unique_ptr<MeshLib::Mesh>>::iterator
Project::findMeshByName(std::string const& name)
{
    return std::find_if(_mesh_vec.begin(), _mesh_vec.end(),
        [&name](std::unique_ptr<MeshLib::Mesh> & mesh)
        { return mesh && (name == mesh->getName()); });
}

const MeshLib::Mesh* Project::getMesh(const std::string &name) const
{
    auto it = findMeshByName(name);
    return (it == _mesh_vec.end() ? nullptr : it->get());
}

bool Project::removeMesh(const std::string &name)
{
    auto it = findMeshByName(name);
    if (it != _mesh_vec.end())
    {
        delete it->release();
        _mesh_vec.erase(it);
        return true;
    }
    return false;
}

bool Project::meshExists(const std::string &name) const
{
    return findMeshByName(name) != _mesh_vec.end();
}

bool Project::getUniqueName(std::string &name) const
{
    int count(0);
    bool isUnique(false);
    std::string cpName;

    while (!isUnique)
    {
        isUnique = true;
        cpName = name;

        count++;
        // If the original name already exists we start to add numbers to name for
        // as long as it takes to make the name unique.
        if (count > 1)
            cpName = cpName + "-" + std::to_string(count);

        for (auto & mesh : _mesh_vec)
            if ( cpName == mesh->getName())
                isUnique = false;
    }

    // At this point cpName is a unique name and isUnique is true.
    // If cpName is not the original name, "name" is changed and isUnique is set to false,
    // indicating that a vector with the original name already exists.
    if (count > 1)
    {
        isUnique = false;
        name = cpName;
    }
    return isUnique;
}

void Project::removePrimaryVariable(std::string const primary_var_name)
{
    std::size_t n_bc(_boundary_conditions.size());
    for (std::size_t i = 0; i<n_bc; ++i)
        if (_boundary_conditions[i]->getProcessVarName() == primary_var_name)
            removeBoundaryCondition(primary_var_name, _boundary_conditions[i]->getParamName());

    std::size_t n_st(_source_terms.size());
    for (std::size_t i = 0; i<n_bc; ++i)
        if (_source_terms[i]->getProcessVarName() == primary_var_name)
            removeSourceTerm(primary_var_name, _source_terms[i]->getParamName());
}

void Project::removeBoundaryCondition(std::string const primary_var_name, std::string const& param_name)
{
    std::size_t n_bc (_boundary_conditions.size());
    for (std::size_t i=0; i<n_bc; ++i)
    {
        if (_boundary_conditions[i]->getProcessVarName() == primary_var_name &&
            _boundary_conditions[i]->getParamName() == param_name)
        {
            BoundaryCondition* bc = _boundary_conditions[i].release();
            delete bc;
            _boundary_conditions.erase(_boundary_conditions.begin()+i);
            return;
        }
    }
}

void Project::removeSourceTerm(std::string const primary_var_name, std::string const& param_name)
{
    std::size_t n_st(_source_terms.size());
    for (std::size_t i = 0; i<n_st; ++i)
    {
        if (_source_terms[i]->getProcessVarName() == primary_var_name &&
            _source_terms[i]->getParamName() == param_name)
        {
            SourceTerm* st = _source_terms[i].release();
            delete st;
            _source_terms.erase(_source_terms.begin() + i);
            return;
        }
    }
}

} //namespace
