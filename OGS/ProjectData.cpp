/**
 * \author Karsten Rink
 * \date   2010-08-25
 * \brief  Implementation of the project data class.
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ProjectData.h"

#include "MeshLib/Mesh.h"

ProjectData::ProjectData()
:
#ifdef OGS_BUILD_GUI
	_geoObjects(new GEOModels())
#else
	_geoObjects(new GeoLib::GEOObjects())
#endif
{
}

ProjectData::~ProjectData()
{
	delete _geoObjects;
	for (std::vector<MeshLib::Mesh*>::iterator it = _msh_vec.begin(); it != _msh_vec.end(); ++it)
		delete *it;
}

void ProjectData::addMesh(MeshLib::Mesh* mesh)
{
	std::string name = mesh->getName();
	isUniqueMeshName(name);
	mesh->setName(name);
	_msh_vec.push_back(mesh);
}

const MeshLib::Mesh* ProjectData::getMesh(const std::string &name) const
{
	for (std::vector<MeshLib::Mesh*>::const_iterator it = _msh_vec.begin(); it != _msh_vec.end(); ++it)
		if (name.compare((*it)->getName()) == 0)
			return *it;
	return NULL;
}

bool ProjectData::removeMesh(const std::string &name)
{
	for (std::vector<MeshLib::Mesh*>::iterator it = _msh_vec.begin(); it != _msh_vec.end(); ++it)
		if (name.compare((*it)->getName()) == 0)
		{
			delete *it;
			_msh_vec.erase(it);
			return true;
		}
	return false;
}

bool ProjectData::meshExists(const std::string &name)
{
	for (std::vector<MeshLib::Mesh*>::const_iterator it = _msh_vec.begin(); it != _msh_vec.end(); ++it)
		if (name.compare((*it)->getName()) == 0)
			return true;
	return false;
}

bool ProjectData::isUniqueMeshName(std::string &name)
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

		for (std::vector<MeshLib::Mesh*>::iterator it = _msh_vec.begin(); it != _msh_vec.end(); ++it)
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
