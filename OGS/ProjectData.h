/**
 * \file
 * \author Karsten Rink
 * \date   2010-08-25
 * \brief
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROJECTDATA_H_
#define PROJECTDATA_H_

#include "FEMCondition.h"
#include "FEMEnums.h"
#include "GEOObjects.h"

namespace MeshLib {
	class Mesh;
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

	// Returns the GEOObjects containing all points, polylines and surfaces
	void setGEOObjects(GeoLib::GEOObjects* geo_objects) { _geoObjects = geo_objects; }

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

	//** Process functionality **//

	/// Adds a new process
	virtual void addProcess(ProcessInfo* pcs);

	/// Returns a process of the given type
	const ProcessInfo* getProcess(FiniteElement::ProcessType type) const;

	/// Removes a process of the given type
	virtual bool removeProcess(FiniteElement::ProcessType type);

	//** FEM Condition functionality **//

	/// Adds a new FEM Condition
	virtual void addCondition(FEMCondition* cond);

	/// Adds new FEM Conditions
	virtual void addConditions(std::vector<FEMCondition*> conds);

	/// Returns the FEM Condition set on a GeoObject with the given name and type from a certain geometry.
	const FEMCondition* getCondition(const std::string &geo_name,
	                                 GeoLib::GEOTYPE type,
	                                 const std::string &cond_name) const;

	/// Returns all FEM Conditions with the given type from a certain geometry.
	const std::vector<FEMCondition*> getConditions(FiniteElement::ProcessType pcs_type = FiniteElement::INVALID_PROCESS,
												   std::string geo_name = "",
												   FEMCondition::CondType type = FEMCondition::UNSPECIFIED) const;

	/// Removes the FEM Condition set on a GeoObject with the given name and type from a certain geometry.
	virtual bool removeCondition(const std::string &geo_name,
	                             GeoLib::GEOTYPE type,
	                             const std::string &cond_name);

	/// Removes all FEM Conditions with the given type from the given process
	virtual void removeConditions(FiniteElement::ProcessType pcs_type = FiniteElement::INVALID_PROCESS,
								  std::string geo_name = "",
								  FEMCondition::CondType cond_type = FEMCondition::UNSPECIFIED);

private:
	GeoLib::GEOObjects* _geoObjects;
	std::vector<MeshLib::Mesh*> _msh_vec;
	std::vector<ProcessInfo*> _pcs_vec;
	std::vector<FEMCondition*> _cond_vec;
};

#endif //PROJECTDATA_H_
