/**
 * \file
 * \brief  Declaration of class Simulation
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <optional>
#ifdef USE_PETSC
#include <vtkMPIController.h>
#include <vtkSmartPointer.h>
#endif

#include "Applications/ApplicationsLib/LinearSolverLibrarySetup.h"
#include "Applications/ApplicationsLib/TestDefinition.h"
#include "BaseLib/ExportSymbol.h"

class ProjectData;
namespace MeshLib
{
class Mesh;
}

class Simulation final
{
public:
    OGS_EXPORT_SYMBOL Simulation(int argc, char* argv[]);

    OGS_EXPORT_SYMBOL void initializeDataStructures(
        std::string const& project,
        std::vector<std::string> const& xml_patch_file_names,
        bool reference_path_is_set, std::string const& reference_path,
        bool nonfatal, std::string const& outdir, std::string const& mesh_dir,
        std::string const& script_dir, bool write_prj);

    double currentTime() const;
    double endTime() const;
    OGS_EXPORT_SYMBOL bool executeTimeStep();
    OGS_EXPORT_SYMBOL bool executeSimulation();
    OGS_EXPORT_SYMBOL void outputLastTimeStep() const;
    OGS_EXPORT_SYMBOL MeshLib::Mesh& getMesh(std::string const& name);

    OGS_EXPORT_SYMBOL std::optional<ApplicationsLib::TestDefinition>
    getTestDefinition() const;

    OGS_EXPORT_SYMBOL ~Simulation();

private:
    ApplicationsLib::LinearSolverLibrarySetup linear_solver_library_setup;
#if defined(USE_PETSC)
    vtkSmartPointer<vtkMPIController> controller;
#endif
    std::unique_ptr<ProjectData> project_data;
    std::optional<ApplicationsLib::TestDefinition> test_definition;
#if defined(OGS_USE_INSITU)
    bool isInsituConfigured = false;
#endif
};
