/**
 * \file
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>
#include <vector>

#include "GeoLib/GEOObjects.h"

#include "MeshLib/Mesh.h"

#include "Applications/DataHolderLib/BoundaryCondition.h"
#include "Applications/DataHolderLib/SourceTerm.h"

//namespace MeshLib {
//    class Mesh;
//}

namespace DataHolderLib
{

/**
 * The ProjectData Object contains all the data needed for a certain project, i.e. all
 * geometric data (stored in a GEOObjects-object), all the meshes, processes,
 * and process variables.
 */
class Project final
{
public:
    /// Constructor
    Project() = default;

    Project(Project&) = delete;

    ~Project() = default;

    /// Returns the GEOObjects containing all points, polylines and surfaces.
    GeoLib::GEOObjects& getGEOObjects() { return _geoObjects; }

    /// Adds a new mesh under a (possibly new) unique name.
    /// \attention This might change the given mesh's name.
    void addMesh(std::unique_ptr<MeshLib::Mesh> mesh);

    /// Returns the mesh with the given name or a \c nullptr if the mesh was not
    /// found.
    const MeshLib::Mesh* getMesh(const std::string &name) const;

    /// Returns all the meshes with their respective names
    const std::vector<std::unique_ptr<MeshLib::Mesh>>& getMeshObjects() const
    {
        return _mesh_vec;
    }

    /// Deletes all meshes with the given name and removes them from the list of
    /// saved meshes. If any mesh was found for removal, true is returned and
    /// false otherwise.
    bool removeMesh(const std::string &name);

    /// Adds a boundary condition to the project
    void addBoundaryCondition(std::unique_ptr<BoundaryCondition> bc)
    {
        _boundary_conditions.push_back(std::move(bc));
    }

    /// Adds a source term to the project
    void addSourceTerm(std::unique_ptr<SourceTerm> st)
    {
        _source_terms.push_back(std::move(st));
    }

    /// Returns the vector of boundary conditions
    std::vector<std::unique_ptr<BoundaryCondition>> const&
    getBoundaryConditions() const
    {
        return _boundary_conditions;
    }

    /// Returns the vector of source terms
    std::vector<std::unique_ptr<SourceTerm>> const& getSourceTerms() const
    {
        return _source_terms;
    }

    /// Removes a primary variable incl. all associated conditions
    void removePrimaryVariable(std::string const primary_var_name);

    /// Removes one boundary condition
    void removeBoundaryCondition(std::string const& primary_var_name,
                                 std::string const& param_name);

    /// Remove one source term
    void removeSourceTerm(std::string const& primary_var_name,
                          std::string const& param_name);

private:
    /// Checks if a mesh with the same name exists and provides a unique name in
    /// case of already existing mesh. Returns true if the mesh name is unique.
    /// Returns false and changes the provided name to a unique name otherwise.
    bool getUniqueName(std::string &name) const;

    /// Returns true if a mesh with the same name exists and false otherwise.
    bool meshExists(const std::string &name) const;


    /// Returns an iterator to the first found mesh with the given name.
    std::vector<std::unique_ptr<MeshLib::Mesh>>::const_iterator findMeshByName(
        std::string const& name) const;
    std::vector<std::unique_ptr<MeshLib::Mesh>>::iterator findMeshByName(
        std::string const& name);

    GeoLib::GEOObjects _geoObjects;
    std::vector<std::unique_ptr<MeshLib::Mesh>> _mesh_vec;
    std::vector<std::unique_ptr<BoundaryCondition>> _boundary_conditions;
    std::vector<std::unique_ptr<SourceTerm>> _source_terms;
};

} // namespace
