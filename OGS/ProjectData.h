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

#ifdef OGS_BUILD_GUI
#include "Gui/DataView/GEOModels.h"
#else
#include "GeoLib/GEOObjects.h"
#endif


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
public:
	ProjectData();
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
	const std::vector<MeshLib::Mesh*>& getMeshObjects() const { return _msh_vec; }

	/// Removes the mesh with the given name.
	virtual bool removeMesh(const std::string &name);

	/// Checks if the name of the mesh is already exists, if so it generates a unique name.
	bool isUniqueMeshName(std::string &name);

	bool meshExists(const std::string &name);

private:
#ifdef OGS_BUILD_GUI
	GEOModels *_geoObjects;
#else
	GeoLib::GEOObjects *_geoObjects;
#endif
	std::vector<MeshLib::Mesh*> _msh_vec;
	std::vector<ProcessLib::Process*> _processes;
};

#endif //PROJECTDATA_H_
