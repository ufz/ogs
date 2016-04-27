/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "Project.h"

#include <algorithm>

// ThirdParty/logog
#include <logog/include/logog.hpp>

#include "BaseLib/FileTools.h"
#include "BaseLib/uniqueInsert.h"

#include "MeshLib/Mesh.h"

namespace DataHolderLib
{

Project::~Project()
{
	for (MeshLib::Mesh* m : _mesh_vec)
		delete m;
}

void Project::addMesh(MeshLib::Mesh* mesh)
{
	std::string name = mesh->getName();
	getUniqueName(name);
	mesh->setName(name);
	_mesh_vec.push_back(mesh);
}

std::vector<MeshLib::Mesh*>::const_iterator Project::findMeshByName(
		std::string const& name) const
{
	return const_cast<Project&>(*this).findMeshByName(name);
}

std::vector<MeshLib::Mesh*>::iterator Project::findMeshByName(
		std::string const& name)
{
	return std::find_if(_mesh_vec.begin(), _mesh_vec.end(),
			[&name](MeshLib::Mesh* mesh)
			{
				return mesh && (name == mesh->getName());
			});
}

const MeshLib::Mesh* Project::getMesh(const std::string &name) const
{
	std::vector<MeshLib::Mesh*>::const_iterator it = findMeshByName(name);
	return (it == _mesh_vec.end() ? nullptr : *it);
}

bool Project::removeMesh(const std::string &name)
{
	std::vector<MeshLib::Mesh*>::iterator it = findMeshByName(name);
	if (it != _mesh_vec.end())
	{
		delete *it;
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

		for (MeshLib::Mesh* mesh :_mesh_vec)
			if ( cpName.compare(mesh->getName()) == 0 )
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

} //namespace
