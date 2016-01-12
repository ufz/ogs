/**
 * \author Karsten Rink
 * \date   2010-08-25
 * \brief  Implementation of the project data class.
 *
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ProjectData.h"

#include <algorithm>

// ThirdParty/logog
#include "logog/include/logog.hpp"

#include "BaseLib/FileTools.h"

#include "MeshLib/Mesh.h"

#include "NumLib/TimeStepping/Algorithms/FixedTimeStepping.h"

// FileIO
#include "FileIO/XmlIO/Boost/BoostXmlGmlInterface.h"
#include "FileIO/readMeshFromFile.h"

#include "BaseLib/ConfigTreeNew.h"

namespace detail
{
static
void readGeometry(std::string const& fname, GeoLib::GEOObjects & geo_objects)
{
	DBUG("Reading geometry file \'%s\'.", fname.c_str());
	FileIO::BoostXmlGmlInterface gml_reader(geo_objects);
	gml_reader.readFile(fname);
}

}

ProjectData::ProjectData(BaseLib::ConfigTree const& project_config,
	std::string const& path)
{
	// geometry
	std::string const geometry_file = BaseLib::copyPathToFileName(
			project_config.get<std::string>("geometry"), path
		);
	detail::readGeometry(geometry_file, *_geoObjects);

	// mesh
	std::string const mesh_file = BaseLib::copyPathToFileName(
			project_config.get<std::string>("mesh"), path
		);

	MeshLib::Mesh* const mesh = FileIO::readMeshFromFile(mesh_file);
	if (!mesh)
		ERR("Could not read mesh from \'%s\' file. No mesh added.",
			mesh_file.c_str());
	_mesh_vec.push_back(mesh);

	// process variables

	BaseLib::ConfigTreeNew var_conf(project_config.get_child("process_variables"));
	parseProcessVariables(var_conf);

	// parameters
	parseParameters(project_config.get_child("parameters"));

	// processes
	parseProcesses(project_config.get_child("processes"));

	// output
	parseOutput(project_config.get_child("output"), path);

	// timestepping
	parseTimeStepping(project_config.get_child("time_stepping"));
}

ProjectData::~ProjectData()
{
	delete _geoObjects;

	for(ProcessLib::Process<GlobalSetupType>* p : _processes)
		delete p;

	for (MeshLib::Mesh* m : _mesh_vec)
		delete m;
}

void ProjectData::addMesh(MeshLib::Mesh* mesh)
{
	std::string name = mesh->getName();
	isMeshNameUniqueAndProvideUniqueName(name);
	mesh->setName(name);
	_mesh_vec.push_back(mesh);
}

std::vector<MeshLib::Mesh*>::const_iterator ProjectData::findMeshByName(
		std::string const& name) const
{
	return const_cast<ProjectData&>(*this).findMeshByName(name);
}

std::vector<MeshLib::Mesh*>::iterator ProjectData::findMeshByName(
		std::string const& name)
{
	return std::find_if(_mesh_vec.begin(), _mesh_vec.end(),
			[&name](MeshLib::Mesh* mesh)
			{
				return mesh && (name == mesh->getName());
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

bool ProjectData::meshExists(const std::string &name) const
{
	return findMeshByName(name) != _mesh_vec.end();
}

bool ProjectData::isMeshNameUniqueAndProvideUniqueName(std::string &name) const
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

void ProjectData::parseProcessVariables(
	BaseLib::ConfigTreeNew& process_variables_config)
{
	DBUG("Parse process variables:")
	if (_geoObjects == nullptr) {
		ERR("Geometric objects are required to define process variables.");
		ERR("No geometric objects present.");
		return;
	}

	// TODO at the moment we have only one mesh, later there
	// can be several meshes. Then we have to check for correct mesh here and
	// assign the referenced mesh below.
	if (_mesh_vec.empty() || _mesh_vec[0] == nullptr) {
		ERR("A mesh is required to define process variables.");
		return;
	}

	// _process_variables.reserve(process_variables_config.size());

	for (auto var_config
		 : process_variables_config.getConfSubtreeList("process_variable")) {
		// TODO Extend to referenced meshes.
		_process_variables.emplace_back(var_config, *_mesh_vec[0], *_geoObjects);
	}
}

void ProjectData::parseParameters(BaseLib::ConfigTree const& parameters_config)
{
	using namespace ProcessLib;

	DBUG("Reading parameters:");
	for (auto pc_it : parameters_config)
	{
		// Skip non-parameter section.
		if (pc_it.first != "parameter")
			continue;
		BaseLib::ConfigTree const& parameter_config = pc_it.second;

		auto name = parameter_config.get_optional<std::string>("name");
		if (!name)
		{
			ERR("The parameter config does not provide a name tag.");
			std::abort();
		}

		auto type = parameter_config.get_optional<std::string>("type");
		if (!type)
		{
			ERR("Could not find required parameter type.");
			std::abort();
		}

		// Create parameter based on the provided type.
		if (*type == "Constant")
		{
			INFO("ConstantParameter: %s.", name->c_str());
			_parameters.push_back(createConstParameter(parameter_config));
			_parameters.back()->name = *name;
		}
		else if (*type == "MeshProperty")
		{
			INFO("MeshPropertyParameter: %s", name->c_str());
			_parameters.push_back(
			    createMeshPropertyParameter(parameter_config, *_mesh_vec[0]));
			_parameters.back()->name = *name;
		}
		else
		{
			ERR("Cannot construct property of given type \'%s\'.",
			    type->c_str());
			std::abort();
		}
	}
}

void ProjectData::parseProcesses(BaseLib::ConfigTree const& processes_config)
{
	DBUG("Reading processes:");
	for (auto pc_it : processes_config) {
		BaseLib::ConfigTree const& process_config = pc_it.second;

		// Check if the process type is specified.
		if (!process_config.get_optional<std::string>("type")) {
			ERR("The process config does not provide a type tag.");
			ERR("   Ignoring this process config.");
			continue;
		}
		_process_configs.push_back(process_config);
	}
}

void ProjectData::parseOutput(BaseLib::ConfigTree const& output_config,
	std::string const& path)
{
	DBUG("Parse output configuration:");

	auto const file = output_config.get_optional<std::string>("file");
	if (!file) {
		ERR("The output config does not provide a file tag.");
		ERR("    Output file not set.");
		return;
	}

	_output_file_prefix = path + *file;
}

void ProjectData::parseTimeStepping(BaseLib::ConfigTree const& timestepping_config)
{
	using namespace ProcessLib;

	DBUG("Reading timestepping configuration.");

	auto const type = timestepping_config.get_optional<std::string>("type");
	if (!type)
	{
		ERR("Could not find required timestepper type.");
		std::abort();
	}

	if (*type == "FixedTimeStepping")
	{
		_time_stepper.reset(NumLib::FixedTimeStepping::newInstance(timestepping_config));
	}
	else if (*type == "SingleStep")
	{
		_time_stepper.reset(new NumLib::FixedTimeStepping(0.0, 1.0, 1.0));
	}
	else
	{
		ERR("Unknown timestepper type: `%s'.", type->c_str());
		std::abort();
	}

	if (!_time_stepper)
	{
		ERR("Initialization of timestepper failed.");
		std::abort();
	}
}



