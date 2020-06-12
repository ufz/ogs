/**
 * \author Karsten Rink
 * \date   2010-08-25
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
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
#include "MaterialLib/MPL/Medium.h"
#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"

#include "ChemistryLib/ChemicalSolverInterface.h"
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
        return processes_;
    }

    ProcessLib::TimeLoop& getTimeLoop() { return *time_loop_; }

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
    void parseMedia(boost::optional<BaseLib::ConfigTree> const& media_config);

    /// Parses the processes configuration and creates new processes for each
    /// process entry passing the corresponding subtree to the process
    /// constructor.
    void parseProcesses(BaseLib::ConfigTree const& processes_config,
                        std::string const& project_directory,
                        std::string const& output_directory);

    /// Parses the time loop configuration.
    void parseTimeLoop(BaseLib::ConfigTree const& config,
                       const std::string& output_directory);

    void parseLinearSolvers(BaseLib::ConfigTree const& config);

    void parseNonlinearSolvers(BaseLib::ConfigTree const& config);

    void parseCurves(boost::optional<BaseLib::ConfigTree> const& config);

    void parseChemicalSolverInterface(
        boost::optional<BaseLib::ConfigTree> const& config,
        const std::string& output_directory);

    std::vector<std::unique_ptr<MeshLib::Mesh>> mesh_vec_;
    std::vector<std::unique_ptr<ProcessLib::Process>> processes_;
    std::vector<ProcessLib::ProcessVariable> process_variables_;

    /// Buffer for each parameter config passed to the process.
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> parameters_;

    boost::optional<ParameterLib::CoordinateSystem> local_coordinate_system_;

    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> media_;

    /// The time loop used to solve this project's processes.
    std::unique_ptr<ProcessLib::TimeLoop> time_loop_;

    std::map<std::string, std::unique_ptr<GlobalLinearSolver>> linear_solvers_;

    std::map<std::string, std::unique_ptr<NumLib::NonlinearSolverBase>>
        nonlinear_solvers_;
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>>
        curves_;

    std::unique_ptr<ChemistryLib::ChemicalSolverInterface>
        chemical_solver_interface_;
};
