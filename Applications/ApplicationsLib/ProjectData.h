/**
 * \author Karsten Rink
 * \date   2010-08-25
 *
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROJECTDATA_H_
#define PROJECTDATA_H_

#include <memory>

#include "BaseLib/ConfigTree.h"

#include "GeoLib/GEOObjects.h"

#include "NumLib/TimeStepping/Algorithms/ITimeStepAlgorithm.h"

#include "ProcessLib/ProcessVariable.h"
#include "ProcessLib/Process.h"
#include "ProcessLib/Parameter.h"
#include "ProcessLib/GroundwaterFlowProcess-fwd.h"

namespace MeshLib {
	class Mesh;
}

/**
 * The ProjectData Object contains all the data needed for a certain project, i.e. all
 * geometric data (stored in a GEOObjects-object), all the meshes, processes,
 * and process variables.
 */
class ProjectData
{
public:
	/// The empty constructor used in the gui, for example, when the project's
	/// configuration is not loaded yet.
	ProjectData() = default;

	/// Constructs project data by parsing provided configuration.
	/// The additional  path is used to find files referenced in the
	/// configuration.
	ProjectData(BaseLib::ConfigTree const& config_tree,
	            std::string const& path);

	ProjectData(ProjectData&) = delete;
	virtual ~ProjectData();

	/// Returns the GEOObjects containing all points, polylines and surfaces.
	GeoLib::GEOObjects* getGEOObjects()
	{
		return _geoObjects;
	}

	/// Adds a new mesh under a (possibly new) unique name.
	/// \attention This might change the given mesh's name.
	void addMesh(MeshLib::Mesh* mesh);

	/// Returns the mesh with the given name or a \c nullptr if the mesh was not
	/// found.
	const MeshLib::Mesh* getMesh(const std::string &name) const;

	/// Returns all the meshes with their respective names
	/// \attention This method should be used only by methods garanteeing
	/// read-only access to the meshes.
	/// \todo This method breaks encapsulation.
	const std::vector<MeshLib::Mesh*>& getMeshObjects() const
	{
		return _mesh_vec;
	}

	/// Deletes all meshes with the given name and removes them from the list of
	/// saved meshes. If any mesh was found for removal, true is returned and
	/// false otherwise.
	bool removeMesh(const std::string &name);

	//
	// Process interface
	//

	/// Builder for processes. The supplied template defines the types of global
	/// vectors and matrices, and the global executor. These types are passed to
	/// every of the constructed processes.
	template <typename GlobalSetupType>
	void buildProcesses()
	{
		for (auto pc : _process_configs)
		{
			if (pc.get<std::string>("type") == "GROUNDWATER_FLOW") {
				// The existence check of the in the configuration referenced
				// process variables is checked in the physical process.
				// TODO at the moment we have only one mesh, later there can be
				// several meshes. Then we have to assign the referenced mesh
				// here.
				_processes.push_back(
					new ProcessLib::GroundwaterFlowProcess<GlobalSetupType>(
						*_mesh_vec[0], _process_variables, _parameters, pc));
			}
			else
			{
				WARN("Unknown process type: %s\n", pc.get<std::string>("type").c_str());
			}
		}
	}

	/// Iterator access for processes.
	/// Provides read access to the process container.
	std::vector<ProcessLib::Process*>::const_iterator
	processesBegin() const
	{
		return _processes.begin();
	}
	std::vector<ProcessLib::Process*>::iterator
	processesBegin()
	{
		return _processes.begin();
	}

	/// Iterator access for processes as in processesBegin().
	std::vector<ProcessLib::Process*>::const_iterator
	processesEnd() const
	{
		return _processes.end();
	}
	std::vector<ProcessLib::Process*>::iterator
	processesEnd()
	{
		return _processes.end();
	}

	std::string const&
	getOutputFilePrefix() const
	{
		return _output_file_prefix;
	}

	NumLib::ITimeStepAlgorithm const& getTimeStepper() const
	{
		return *_time_stepper;
	}
	NumLib::ITimeStepAlgorithm& getTimeStepper()
	{
		return *_time_stepper;
	}

private:
	/// Checks if a mesh with the same name exists and provides a unique name in
	/// case of already existing mesh. Returns true if the mesh name is unique.
	/// Returns false and changes the provided name to a unique name otherwise.
	bool isMeshNameUniqueAndProvideUniqueName(std::string &name) const;

	/// Returns true if a mesh with the same name exists and false otherwise.
	bool meshExists(const std::string &name) const;


	/// Returns an iterator to the first found mesh with the given name.
	std::vector<MeshLib::Mesh*>::const_iterator findMeshByName(
		std::string const& name) const;
	std::vector<MeshLib::Mesh*>::iterator findMeshByName(
		std::string const& name);

	/// Parses the process variables configuration and creates new variables for
	/// each variable entry passing the corresponding subtree to the process
	/// variable constructor.
	void parseProcessVariables(BaseLib::ConfigTreeNew& process_variables_config);

	/// Parses the parameters configuration and saves them in a list.
	/// Checks if a parameter has name tag.
	void parseParameters(BaseLib::ConfigTree const& parameters_config);

	/// Parses the processes configuration and creates new processes for each
	/// process entry passing the corresponding subtree to the process
	/// constructor.
	void parseProcesses(BaseLib::ConfigTree const& process_config);

	/// Parses the output configuration.
	/// Parses the file tag and sets output file prefix.
	void parseOutput(BaseLib::ConfigTree const& output_config, std::string const& path);

	void parseTimeStepping(BaseLib::ConfigTree const& timestepping_config);

private:
	GeoLib::GEOObjects *_geoObjects = new GeoLib::GEOObjects();
	std::vector<MeshLib::Mesh*> _mesh_vec;
	std::vector<ProcessLib::Process*> _processes;
	std::vector<ProcessLib::ProcessVariable> _process_variables;

	/// Buffer for each process' config used in the process building function.
	std::vector<BaseLib::ConfigTree> _process_configs;

	/// Buffer for each parameter config passed to the process.
	std::vector<std::unique_ptr<ProcessLib::ParameterBase>> _parameters;

	/// Output file path with project prefix.
	std::string _output_file_prefix;

	/// Timestepper
	std::unique_ptr<NumLib::ITimeStepAlgorithm> _time_stepper;
};

#endif //PROJECTDATA_H_
