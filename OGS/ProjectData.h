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
#include "Gui/DataView/GEOModels.h"
#else
#include "GeoLib/GEOObjects.h"
#endif

#include "OGS/ProcessVariable.h"
namespace MeshLib {
	class Mesh;
}

namespace ProcessLib {
	class Process;
}

/**
 * The ProjectData Object contains all the data needed for a certain project, i.e. all
 * geometric data (stored in a GEOObjects-object), all the meshes, FEM Conditions (i.e.
 * Boundary Conditions, Source Terms and Initial Conditions), etc.
 * ProjectData does not administrate any of the objects, it is just a "container class"
 * to store them all in one place.
 * For each class of object stored in this container exists an add-, get- and remove-method.
 *
 * \sa GEOModels, FEMCondition
 */
class ProjectData
{
using ConfigTree = boost::property_tree::ptree;
public:
	ProjectData();
	ProjectData(ConfigTree const& config_tree, std::string const& path);
	virtual ~ProjectData();

	//** Geometry functionality **//

	// Returns the GEOObjects containing all points, polylines and surfaces
	GeoLib::GEOObjects* getGEOObjects() { return _geoObjects; }

	//** Mesh functionality **//

	/// Adds a new mesh
	virtual void addMesh(MeshLib::Mesh* mesh);

	/// Returns the mesh with the given name.
	const MeshLib::Mesh* getMesh(const std::string &name) const;

	/// Returns all the meshes with their respective names
	const std::vector<MeshLib::Mesh*>& getMeshObjects() const { return _mesh_vec; }

	/// Deletes all meshes with the given name and removes them from the list of
	//  saved meshes.
	virtual bool removeMesh(const std::string &name);

	/// Checks if the name of the mesh is already exists, if so it generates a unique name.
	bool isUniqueMeshName(std::string &name);

	bool meshExists(const std::string &name);

private:
	std::vector<MeshLib::Mesh*>::const_iterator findMeshByName(
		std::string const& name) const;
	std::vector<MeshLib::Mesh*>::iterator findMeshByName(
		std::string const& name);

	/// read the process variables from configuration
	void readProcessVariables(ConfigTree const& process_variables_config);

	// read the processes from configuration
	void readProcesses(ConfigTree const& process_config);

private:
#ifdef OGS_BUILD_GUI
	GEOModels *_geoObjects;
#else
	GeoLib::GEOObjects *_geoObjects;
#endif
	std::vector<MeshLib::Mesh*> _mesh_vec;
	std::vector<ProcessLib::Process*> _processes;
	std::vector<OGS::ProcessVariable> _process_variables;
};

#endif //PROJECTDATA_H_
