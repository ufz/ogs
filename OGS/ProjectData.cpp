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

// ThirdParty/logog
#include "logog/include/logog.hpp"

#include "BaseLib/FileTools.h"

#include "MeshLib/Mesh.h"
#include "ProcessLib/Process.h"

// FileIO
#include "FileIO/XmlIO/Boost/BoostXmlGmlInterface.h"
#include "FileIO/readMeshFromFile.h"

namespace detail
{
void readGeometry(std::string const& fname, GeoLib::GEOObjects & geo_objects)

{
	DBUG("Reading geometry file \'%s\'.", fname.c_str());
	FileIO::BoostXmlGmlInterface gml_reader(geo_objects);
	gml_reader.readFile(fname);
}
}

ProjectData::ProjectData()
:
#ifdef OGS_BUILD_GUI
	_geoObjects(new GEOModels())
#else
	_geoObjects(new GeoLib::GEOObjects())
#endif
{
}

ProjectData::ProjectData(ConfigTree const& project_config,
	std::string const& path)
:
#ifdef OGS_BUILD_GUI
	_geoObjects(new GEOModels())
#else
	_geoObjects(new GeoLib::GEOObjects())
#endif
{
	// read geometry
	std::string const geometry_file = BaseLib::copyPathToFileName(
			project_config.get<std::string>("geometry"), path
		);
	detail::readGeometry(geometry_file, *_geoObjects);

	// read mesh
	std::string const mesh_file = BaseLib::copyPathToFileName(
			project_config.get<std::string>("mesh"), path
		);
	_mesh_vec.push_back(FileIO::readMeshFromFile(mesh_file));
}

ProjectData::~ProjectData()
{
	delete _geoObjects;

	for(ProcessLib::Process* p : _processes)
		delete p;

	for (MeshLib::Mesh* m : _mesh_vec)
		delete m;
}

void ProjectData::addMesh(MeshLib::Mesh* mesh)
{
	std::string name = mesh->getName();
	isUniqueMeshName(name);
	mesh->setName(name);
	_mesh_vec.push_back(mesh);
}

std::vector<MeshLib::Mesh*>::const_iterator ProjectData::findMeshByName(
		std::string const& name) const
{
	return findMeshByName(name);
}

std::vector<MeshLib::Mesh*>::iterator ProjectData::findMeshByName(
		std::string const& name)
{
	return std::find_if(_mesh_vec.begin(), _mesh_vec.end(),
			[&name](MeshLib::Mesh* mesh)
			{
				return name == mesh->getName();
			});

}

const MeshLib::Mesh* ProjectData::getMesh(const std::string &name) const
{
	std::vector<MeshLib::Mesh*>::const_iterator it = findMeshByName(name);
	return (it == _mesh_vec.end() ? nullptr : *it);
}

bool ProjectData::removeMesh(const std::string &name)
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

bool ProjectData::meshExists(const std::string &name)
{
	return findMeshByName(name) != _mesh_vec.end();
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

		for (std::vector<MeshLib::Mesh*>::iterator it = _mesh_vec.begin();
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
