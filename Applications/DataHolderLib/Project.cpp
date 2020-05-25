/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "Project.h"

#include <algorithm>

#include "BaseLib/Algorithm.h"
#include "BaseLib/FileTools.h"
#include "MeshLib/Mesh.h"

namespace DataHolderLib
{

void Project::addMesh(std::unique_ptr<MeshLib::Mesh> mesh)
{
    std::string name = mesh->getName();
    getUniqueName(name);
    mesh->setName(name);
    mesh_vec_.push_back(std::move(mesh));
}

std::vector<std::unique_ptr<MeshLib::Mesh>>::const_iterator
Project::findMeshByName(std::string const& name) const
{
    return const_cast<Project&>(*this).findMeshByName(name);
}

std::vector<std::unique_ptr<MeshLib::Mesh>>::iterator
Project::findMeshByName(std::string const& name)
{
    return std::find_if(mesh_vec_.begin(), mesh_vec_.end(),
        [&name](std::unique_ptr<MeshLib::Mesh> & mesh)
        { return mesh && (name == mesh->getName()); });
}

const MeshLib::Mesh* Project::getMesh(const std::string &name) const
{
    auto it = findMeshByName(name);
    return (it == mesh_vec_.end() ? nullptr : it->get());
}

bool Project::removeMesh(const std::string &name)
{
    auto it = findMeshByName(name);
    if (it != mesh_vec_.end())
    {
        delete it->release();
        mesh_vec_.erase(it);
        return true;
    }
    return false;
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
        {
            cpName = cpName + "-" + std::to_string(count);
        }

        for (auto& mesh : mesh_vec_)
        {
            if (cpName == mesh->getName())
            {
                isUnique = false;
            }
        }
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

void Project::removePrimaryVariable(std::string const& primary_var_name)
{
    std::size_t const n_bc(boundary_conditions_.size());
    for (int i = n_bc - 1; i >= 0; --i)
    {
        if (boundary_conditions_[i]->getProcessVarName() == primary_var_name)
        {
            removeBoundaryCondition(primary_var_name,
                                    boundary_conditions_[i]->getParamName());
        }
    }

    std::size_t const n_st(source_terms_.size());
    for (int i = n_st - 1; i >= 0; --i)
    {
        if (source_terms_[i]->getProcessVarName() == primary_var_name)
        {
            removeSourceTerm(primary_var_name,
                             source_terms_[i]->getParamName());
        }
    }
}

void Project::removeBoundaryCondition(std::string const& primary_var_name,
                                      std::string const& param_name)
{
    std::size_t const n_bc(boundary_conditions_.size());
    for (std::size_t i = 0; i < n_bc; ++i)
    {
        if (boundary_conditions_[i]->getProcessVarName() == primary_var_name &&
            boundary_conditions_[i]->getParamName() == param_name)
        {
            boundary_conditions_[i].reset();
            boundary_conditions_.erase(boundary_conditions_.begin() + i);
            return;
        }
    }
}

void Project::removeSourceTerm(std::string const& primary_var_name,
                               std::string const& param_name)
{
    std::size_t const n_st(source_terms_.size());
    for (std::size_t i = 0; i < n_st; ++i)
    {
        if (source_terms_[i]->getProcessVarName() == primary_var_name &&
            source_terms_[i]->getParamName() == param_name)
        {
            source_terms_[i].reset();
            source_terms_.erase(source_terms_.begin() + i);
            return;
        }
    }
}

}  // namespace DataHolderLib
