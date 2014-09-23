/**
 * \author Karsten Rink
 * \date   2010-08-25
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROJECTDATA_H_
#define PROJECTDATA_H_

#include <boost/property_tree/ptree.hpp>

#ifdef OGS_BUILD_GUI
#include "GEOModels.h"
#else
#include "GeoLib/GEOObjects.h"
#endif

#include "ProcessLib/ProcessVariable.h"
namespace MeshLib {
	class Mesh;
}

namespace ProcessLib {
	class Process;
}

/**
 * The ProjectData Object contains all the data needed for a certain project, i.e. all
 * geometric data (stored in a GEOObjects-object), all the meshes, processes,
 * and process variables.
 */
class ProjectData
{
using ConfigTree = boost::property_tree::ptree;
public:
	/// The empty constructor used in the gui, for example, when the project's
	/// configuration is not loaded yet.
	ProjectData() = default;

	/// Constructs project data by parsing provided configuration.
	/// The additional  path is used to find files referenced in the
	/// configuration.
	ProjectData(ConfigTree const& config_tree, std::string const& path);

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
	void parseProcessVariables(ConfigTree const& process_variables_config);

	/// Parses the processes configuration and creates new processes for each
	/// process entry passing the corresponding subtree to the process
	/// constructor.
	void parseProcesses(ConfigTree const& process_config);

private:
#ifdef OGS_BUILD_GUI
	GEOModels *_geoObjects = new GEOModels();
#else
	GeoLib::GEOObjects *_geoObjects = new GeoLib::GEOObjects();
#endif
	std::vector<MeshLib::Mesh*> _mesh_vec;
	std::vector<ProcessLib::Process*> _processes;
	std::vector<ProcessLib::ProcessVariable> _process_variables;
};

#endif //PROJECTDATA_H_
