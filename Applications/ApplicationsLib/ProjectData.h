/**
 * \author Karsten Rink
 * \date   2010-08-25
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <boost/optional/optional.hpp>
#include <map>
#include <memory>
#include <string>

#include "BaseLib/ConfigTree.h"

#include "GeoLib/GEOObjects.h"

#include "ProcessLib/Parameter/Parameter.h"
#include "ProcessLib/Process.h"
#include "ProcessLib/ProcessVariable.h"

namespace MathLib
{
class PiecewiseLinearInterpolation;
}
namespace MeshLib
{
class Mesh;
}

namespace NumLib
{
class NonlinearSolverBase;
}

namespace ProcessLib
{
class UncoupledProcessesTimeLoop;
}

/**
 * The ProjectData Object contains all the data needed for a certain project,
 * i.e. all geometric data (stored in a GEOObjects-object), all the meshes,
 * processes, and process variables.
 */
class ProjectData final
{
public:
    /// The time loop type used to solve this project's processes.
    using TimeLoop = ProcessLib::UncoupledProcessesTimeLoop;

    /// The empty constructor used in the gui, for example, when the project's
    /// configuration is not loaded yet.
    ProjectData();

    /// Constructs project data by parsing provided configuration.
    ///
    /// \param config_tree Configuration as read from the prj file.
    /// \param project_directory Where to look for files referenced in the
    ///                          \c config_tree.
    /// \param output_directory  Where to write simulation output files to.
    ProjectData(BaseLib::ConfigTree const& config_tree,
                std::string const& project_directory,
                std::string const& output_directory);

    ProjectData(ProjectData&) = delete;

    ~ProjectData();

    /// Returns the GEOObjects containing all points, polylines and surfaces.
    GeoLib::GEOObjects* getGEOObjects() { return _geoObjects; }
    /// Adds a new mesh under a (possibly new) unique name.
    /// \attention This might change the given mesh's name.
    void addMesh(MeshLib::Mesh* mesh);

    /// Returns the mesh with the given name or a \c nullptr if the mesh was not
    /// found.
    const MeshLib::Mesh* getMesh(const std::string& name) const;

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
    bool removeMesh(const std::string& name);

    //
    // Process interface
    //

    /// Provides read access to the process container.
    std::map<std::string, std::unique_ptr<ProcessLib::Process>> const&
    getProcesses() const
    {
        return _processes;
    }

    TimeLoop& getTimeLoop() { return *_time_loop; }
private:
    /// Checks if a mesh with the same name exists and provides a unique name in
    /// case of already existing mesh. Returns true if the mesh name is unique.
    /// Returns false and changes the provided name to a unique name otherwise.
    bool isMeshNameUniqueAndProvideUniqueName(std::string& name) const;

    /// Returns true if a mesh with the same name exists and false otherwise.
    bool meshExists(const std::string& name) const;

    /// Returns an iterator to the first found mesh with the given name.
    std::vector<MeshLib::Mesh*>::const_iterator findMeshByName(
        std::string const& name) const;
    std::vector<MeshLib::Mesh*>::iterator findMeshByName(
        std::string const& name);

    /// Parses the process variables configuration and creates new variables for
    /// each variable entry passing the corresponding subtree to the process
    /// variable constructor.
    void parseProcessVariables(
        BaseLib::ConfigTree const& process_variables_config);

    /// Parses the parameters configuration and saves them in a list.
    /// Checks if a parameter has name tag.
    void parseParameters(BaseLib::ConfigTree const& parameters_config);

    /// Parses the processes configuration and creates new processes for each
    /// process entry passing the corresponding subtree to the process
    /// constructor.
    void parseProcesses(BaseLib::ConfigTree const& process_config,
                        std::string const& project_directory,
                        std::string const& output_directory);

    /// Parses the time loop configuration.
    void parseTimeLoop(BaseLib::ConfigTree const& config,
                       const std::string& output_directory);

    void parseLinearSolvers(BaseLib::ConfigTree const& config);

    void parseNonlinearSolvers(BaseLib::ConfigTree const& config);

    void parseCurves(boost::optional<BaseLib::ConfigTree> const& config);

    GeoLib::GEOObjects* _geoObjects = new GeoLib::GEOObjects();
    std::vector<MeshLib::Mesh*> _mesh_vec;
    std::map<std::string, std::unique_ptr<ProcessLib::Process>> _processes;
    std::vector<ProcessLib::ProcessVariable> _process_variables;

    /// Buffer for each parameter config passed to the process.
    std::vector<std::unique_ptr<ProcessLib::ParameterBase>> _parameters;

    /// The time loop used to solve this project's processes.
    std::unique_ptr<TimeLoop> _time_loop;

    std::map<std::string, std::unique_ptr<GlobalLinearSolver>> _linear_solvers;

    std::map<std::string, std::unique_ptr<NumLib::NonlinearSolverBase>>
        _nonlinear_solvers;
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>>
        _curves;
};
