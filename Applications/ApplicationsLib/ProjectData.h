/**
 * \author Karsten Rink
 * \date   2010-08-25
 *
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <map>
#include <memory>
#include <optional>
#include <string>

#include "BaseLib/ConfigTree.h"
#include "ChemistryLib/ChemicalSolverInterface.h"
#include "MaterialLib/MPL/Medium.h"
#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
#include "ParameterLib/CoordinateSystem.h"
#include "ParameterLib/Parameter.h"
#include "ProcessLib/Process.h"
#include "ProcessLib/ProcessVariable.h"

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
class TimeLoop;
}

/**
 * The ProjectData Object contains all the data needed for a certain project,
 * i.e. all geometric data (stored in a GEOObjects-object), all the meshes,
 * processes, and process variables.
 */
class ProjectData final
{
public:
    /// The empty constructor used in the gui, for example, when the project's
    /// configuration is not loaded yet.
    ProjectData();

    /// Constructs project data by parsing provided configuration.
    ///
    /// \param project_config Configuration as read from the prj file.
    /// \param project_directory Where to look for files referenced in the
    ///                          \c config_tree.
    /// \param output_directory  Where to write simulation output files to.
    ProjectData(BaseLib::ConfigTree const& project_config,
                std::string const& project_directory,
                std::string const& output_directory);

    ProjectData(ProjectData&) = delete;


    //
    // Process interface
    //

    /// Provides read access to the process container.
    std::vector<std::unique_ptr<ProcessLib::Process>> const& getProcesses()
        const
    {
        return _processes;
    }

    ProcessLib::TimeLoop& getTimeLoop() { return *_time_loop; }

private:
    /// Parses the process variables configuration and creates new variables for
    /// each variable entry passing the corresponding subtree to the process
    /// variable constructor.
    void parseProcessVariables(
        BaseLib::ConfigTree const& process_variables_config);

    /// Parses the parameters configuration and saves them.
    /// Checks for double parameters' names. Returns names of vectors which are
    /// to be transformed using local coordinate system.
    std::vector<std::string> parseParameters(
        BaseLib::ConfigTree const& parameters_config);

    /// Parses media configuration and saves them in an object.
    void parseMedia(std::optional<BaseLib::ConfigTree> const& media_config);

    /// Parses the processes configuration and creates new processes for each
    /// process entry passing the corresponding subtree to the process
    /// constructor.
    void parseProcesses(BaseLib::ConfigTree const& processes_config,
                        std::string const& project_directory,
                        std::string const& output_directory,
                        std::unique_ptr<ChemistryLib::ChemicalSolverInterface>&&
                            chemical_solver_interface);

    /// Parses the time loop configuration.
    void parseTimeLoop(BaseLib::ConfigTree const& config,
                       const std::string& output_directory);

    void parseLinearSolvers(BaseLib::ConfigTree const& config);

    void parseNonlinearSolvers(BaseLib::ConfigTree const& config);

    void parseCurves(std::optional<BaseLib::ConfigTree> const& config);

    std::unique_ptr<ChemistryLib::ChemicalSolverInterface>
    parseChemicalSolverInterface(
        std::optional<BaseLib::ConfigTree> const& config,
        const std::string& output_directory);

    std::vector<std::unique_ptr<MeshLib::Mesh>> _mesh_vec;
    std::vector<std::unique_ptr<ProcessLib::Process>> _processes;
    std::vector<ProcessLib::ProcessVariable> _process_variables;

    /// Buffer for each parameter config passed to the process.
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> _parameters;

    std::optional<ParameterLib::CoordinateSystem> _local_coordinate_system;

    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> _media;

    /// The time loop used to solve this project's processes.
    std::unique_ptr<ProcessLib::TimeLoop> _time_loop;

    std::map<std::string, std::unique_ptr<GlobalLinearSolver>> _linear_solvers;

    std::map<std::string, std::unique_ptr<NumLib::NonlinearSolverBase>>
        _nonlinear_solvers;
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>>
        _curves;
};

/// Parses a comma separated list of integers.
/// Such lists occur in the medium definition in the OGS prj file.
/// Error messages in this function refer to this specific purpose.
std::vector<int> splitMaterialIdString(std::string const& material_id_string);
