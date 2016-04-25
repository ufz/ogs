/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "DataExplorerProject.h"

#include <algorithm>

// ThirdParty/logog
#include "logog/include/logog.hpp"

#include "BaseLib/FileTools.h"
#include "BaseLib/uniqueInsert.h"

#include "MeshLib/Mesh.h"

DataExplorerProject::DataExplorerProject() = default;

DataExplorerProject::~DataExplorerProject()
{
	delete _geoObjects;

	for (MeshLib::Mesh* m : _mesh_vec)
		delete m;
}

void DataExplorerProject::addMesh(MeshLib::Mesh* mesh)
{
	std::string name = mesh->getName();
	isMeshNameUniqueAndProvideUniqueName(name);
	mesh->setName(name);
	_mesh_vec.push_back(mesh);
}

std::vector<MeshLib::Mesh*>::const_iterator DataExplorerProject::findMeshByName(
		std::string const& name) const
{
	return const_cast<DataExplorerProject&>(*this).findMeshByName(name);
}

std::vector<MeshLib::Mesh*>::iterator DataExplorerProject::findMeshByName(
		std::string const& name)
{
	return std::find_if(_mesh_vec.begin(), _mesh_vec.end(),
			[&name](MeshLib::Mesh* mesh)
			{
				return mesh && (name == mesh->getName());
			});
}

const MeshLib::Mesh* DataExplorerProject::getMesh(const std::string &name) const
{
	std::vector<MeshLib::Mesh*>::const_iterator it = findMeshByName(name);
	return (it == _mesh_vec.end() ? nullptr : *it);
}

bool DataExplorerProject::removeMesh(const std::string &name)
{
	bool mesh_found = false;
	std::vector<MeshLib::Mesh*>::iterator it = findMeshByName(name);
	while (it != _mesh_vec.end())
	{
		delete *it;
		*it = nullptr;
		it = findMeshByName(name);
		mesh_found = true;
	}

	_mesh_vec.erase(std::remove(_mesh_vec.begin(), _mesh_vec.end(), nullptr),
			_mesh_vec.end());
	return mesh_found;
}

bool DataExplorerProject::meshExists(const std::string &name) const
{
	return findMeshByName(name) != _mesh_vec.end();
}

bool DataExplorerProject::isMeshNameUniqueAndProvideUniqueName(std::string &name) const
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

		for (std::vector<MeshLib::Mesh*>::const_iterator it = _mesh_vec.begin();
				it != _mesh_vec.end(); ++it)
			if ( cpName.compare((*it)->getName()) == 0 )
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
