/**
 * \author Karsten Rink
 * \date   2010-08-25
 * \brief  Implementation of the project data class.
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ProjectData.h"

#include <algorithm>
#include <set>

#include <logog/include/logog.hpp>

#include "BaseLib/FileTools.h"
#include "BaseLib/uniqueInsert.h"

#include "MathLib/Curve/CreatePiecewiseLinearCurve.h"
#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
#include "MeshLib/Mesh.h"

#include "NumLib/ODESolver/ConvergenceCriterion.h"
#include "ProcessLib/CreateJacobianAssembler.h"

// FileIO
#include "GeoLib/IO/XmlIO/Boost/BoostXmlGmlInterface.h"
#include "MeshLib/IO/readMeshFromFile.h"

#include "ProcessLib/UncoupledProcessesTimeLoop.h"

#ifdef OGS_BUILD_PROCESS_COMPONENTTRANSPORT
#include "ProcessLib/ComponentTransport/CreateComponentTransportProcess.h"
#endif
#ifdef OGS_BUILD_PROCESS_GROUNDWATERFLOW
#include "ProcessLib/GroundwaterFlow/CreateGroundwaterFlowProcess.h"
#endif
#ifdef OGS_BUILD_PROCESS_HT
#include "ProcessLib/HT/CreateHTProcess.h"
#endif
#ifdef OGS_BUILD_PROCESS_HEATCONDUCTION
#include "ProcessLib/HeatConduction/CreateHeatConductionProcess.h"
#endif
#ifdef OGS_BUILD_PROCESS_HYDROMECHANICS
#include "ProcessLib/HydroMechanics/CreateHydroMechanicsProcess.h"
#endif
#ifdef OGS_BUILD_PROCESS_LIE
#include "ProcessLib/LIE/HydroMechanics/CreateHydroMechanicsProcess.h"
#include "ProcessLib/LIE/SmallDeformation/CreateSmallDeformationProcess.h"
#endif
#ifdef OGS_BUILD_PROCESS_LIQUIDFLOW
#include "ProcessLib/LiquidFlow/CreateLiquidFlowProcess.h"
#endif
#ifdef OGS_BUILD_PROCESS_PHASEFIELD
#include "ProcessLib/PhaseField/CreatePhaseFieldProcess.h"
#endif
#ifdef OGS_BUILD_PROCESS_RICHARDSCOMPONENTTRANSPORT
#include "ProcessLib/RichardsComponentTransport/CreateRichardsComponentTransportProcess.h"
#endif
#ifdef OGS_BUILD_PROCESS_RICHARDSFLOW
#include "ProcessLib/RichardsFlow/CreateRichardsFlowProcess.h"
#endif
#ifdef OGS_BUILD_PROCESS_SMALLDEFORMATION
#include "ProcessLib/SmallDeformation/CreateSmallDeformationProcess.h"
#endif
#ifdef OGS_BUILD_PROCESS_TES
#include "ProcessLib/TES/CreateTESProcess.h"
#endif
#ifdef OGS_BUILD_PROCESS_THERMALTWOPHASEFLOWWITHPP
#include "ProcessLib/ThermalTwoPhaseFlowWithPP/CreateThermalTwoPhaseFlowWithPPProcess.h"
#endif
#ifdef OGS_BUILD_PROCESS_THERMOMECHANICALPHASEFIELD
#include "ProcessLib/ThermoMechanicalPhaseField/CreateThermoMechanicalPhaseFieldProcess.h"
#endif
#ifdef OGS_BUILD_PROCESS_THERMOMECHANICS
#include "ProcessLib/ThermoMechanics/CreateThermoMechanicsProcess.h"
#endif
#ifdef OGS_BUILD_PROCESS_TWOPHASEFLOWWITHPP
#include "ProcessLib/TwoPhaseFlowWithPP/CreateTwoPhaseFlowWithPPProcess.h"
#endif
#ifdef OGS_BUILD_PROCESS_TWOPHASEFLOWWITHPRHO
#include "ProcessLib/TwoPhaseFlowWithPrho/CreateTwoPhaseFlowWithPrhoProcess.h"
#endif

namespace detail
{
static void readGeometry(std::string const& fname,
                         GeoLib::GEOObjects& geo_objects)
{
    DBUG("Reading geometry file \'%s\'.", fname.c_str());
    GeoLib::IO::BoostXmlGmlInterface gml_reader(geo_objects);
    gml_reader.readFile(fname);
}
}

ProjectData::ProjectData() = default;

ProjectData::ProjectData(BaseLib::ConfigTree const& project_config,
                         std::string const& project_directory,
                         std::string const& output_directory)
{
    std::string const geometry_file = BaseLib::copyPathToFileName(
        //! \ogs_file_param{prj__geometry}
        project_config.getConfigParameter<std::string>("geometry"),
        project_directory);
    detail::readGeometry(geometry_file, *_geoObjects);

    {
        //! \ogs_file_param{prj__mesh}
        auto const mesh_param = project_config.getConfigParameter("mesh");

        std::string const mesh_file = BaseLib::copyPathToFileName(
            mesh_param.getValue<std::string>(), project_directory);

        MeshLib::Mesh* const mesh = MeshLib::IO::readMeshFromFile(mesh_file);
        if (!mesh)
        {
            OGS_FATAL("Could not read mesh from \'%s\' file. No mesh added.",
                      mesh_file.c_str());
        }

        if (auto const axially_symmetric =
                //! \ogs_file_attr{prj__mesh__axially_symmetric}
            mesh_param.getConfigAttributeOptional<bool>("axially_symmetric"))
        {
            mesh->setAxiallySymmetric(*axially_symmetric);
        }
        _mesh_vec.push_back(mesh);
    }

    //! \ogs_file_param{prj__curves}
    parseCurves(project_config.getConfigSubtreeOptional("curves"));

    //! \ogs_file_param{prj__parameters}
    parseParameters(project_config.getConfigSubtree("parameters"));

    //! \ogs_file_param{prj__process_variables}
    parseProcessVariables(project_config.getConfigSubtree("process_variables"));

    //! \ogs_file_param{prj__processes}
    parseProcesses(project_config.getConfigSubtree("processes"),
                   project_directory, output_directory);

    //! \ogs_file_param{prj__linear_solvers}
    parseLinearSolvers(project_config.getConfigSubtree("linear_solvers"));

    //! \ogs_file_param{prj__nonlinear_solvers}
    parseNonlinearSolvers(project_config.getConfigSubtree("nonlinear_solvers"));

    //! \ogs_file_param{prj__time_loop}
    parseTimeLoop(project_config.getConfigSubtree("time_loop"),
                  output_directory);
}

ProjectData::~ProjectData()
{
    delete _geoObjects;

    for (MeshLib::Mesh* m : _mesh_vec)
        delete m;
}

void ProjectData::addMesh(MeshLib::Mesh* mesh)
{
    std::string name = mesh->getName();
    isMeshNameUniqueAndProvideUniqueName(name);
    mesh->setName(name);
    _mesh_vec.push_back(mesh);
}

std::vector<MeshLib::Mesh*>::const_iterator ProjectData::findMeshByName(
    std::string const& name) const
{
    return const_cast<ProjectData&>(*this).findMeshByName(name);
}

std::vector<MeshLib::Mesh*>::iterator ProjectData::findMeshByName(
    std::string const& name)
{
    return std::find_if(_mesh_vec.begin(), _mesh_vec.end(),
                        [&name](MeshLib::Mesh* mesh) {
                            return mesh && (name == mesh->getName());
                        });
}

const MeshLib::Mesh* ProjectData::getMesh(const std::string& name) const
{
    auto it = findMeshByName(name);
    return (it == _mesh_vec.end() ? nullptr : *it);
}

bool ProjectData::removeMesh(const std::string& name)
{
    bool mesh_found = false;
    auto it = findMeshByName(name);
    while (it != _mesh_vec.end())
    {
        delete *it;
        *it = nullptr;
        it = findMeshByName(name);
        mesh_found = true;
    }

    _mesh_vec.erase(std::remove(_mesh_vec.begin(), _mesh_vec.end(), nullptr),
                    _mesh_vec.end());
    return mesh_found;
}

bool ProjectData::meshExists(const std::string& name) const
{
    return findMeshByName(name) != _mesh_vec.end();
}

bool ProjectData::isMeshNameUniqueAndProvideUniqueName(std::string& name) const
{
    int count(0);
    bool isUnique(false);
    std::string cpName;

    while (!isUnique)
    {
        isUnique = true;
        cpName = name;

        count++;
        // If the original name already exists we start to add numbers to name
        // for
        // as long as it takes to make the name unique.
        if (count > 1)
            cpName = cpName + "-" + std::to_string(count);

        for (auto mesh : _mesh_vec)
            if (cpName == mesh->getName())
                isUnique = false;
    }

    // At this point cpName is a unique name and isUnique is true.
    // If cpName is not the original name, "name" is changed and isUnique is set
    // to false,
    // indicating that a vector with the original name already exists.
    if (count > 1)
    {
        isUnique = false;
        name = cpName;
    }
    return isUnique;
}

void ProjectData::parseProcessVariables(
    BaseLib::ConfigTree const& process_variables_config)
{
    DBUG("Parse process variables:")
    if (_geoObjects == nullptr)
    {
        ERR("Geometric objects are required to define process variables.");
        ERR("No geometric objects present.");
        return;
    }

    // TODO at the moment we have only one mesh, later there
    // can be several meshes. Then we have to check for correct mesh here and
    // assign the referenced mesh below.
    if (_mesh_vec.empty() || _mesh_vec[0] == nullptr)
    {
        ERR("A mesh is required to define process variables.");
        return;
    }

    std::set<std::string> names;

    for (auto var_config
         //! \ogs_file_param{prj__process_variables__process_variable}
         : process_variables_config.getConfigSubtreeList("process_variable"))
    {
        // TODO Extend to referenced meshes.
        auto pv = ProcessLib::ProcessVariable{var_config, *_mesh_vec[0],
                                              *_geoObjects, _parameters};
        if (!names.insert(pv.getName()).second)
            OGS_FATAL("A process variable with name `%s' already exists.",
                      pv.getName().c_str());

        _process_variables.push_back(std::move(pv));
    }
}

void ProjectData::parseParameters(BaseLib::ConfigTree const& parameters_config)
{
    using namespace ProcessLib;

    std::set<std::string> names;

    DBUG("Reading parameters:");
    for (auto parameter_config :
         //! \ogs_file_param{prj__parameters__parameter}
         parameters_config.getConfigSubtreeList("parameter"))
    {
        auto p =
            ProcessLib::createParameter(parameter_config, _mesh_vec, _curves);
        if (!names.insert(p->name).second)
            OGS_FATAL("A parameter with name `%s' already exists.",
                      p->name.c_str());

        _parameters.push_back(std::move(p));
    }

    for (auto& parameter : _parameters)
        parameter->initialize(_parameters);
}

void ProjectData::parseProcesses(BaseLib::ConfigTree const& processes_config,
                                 std::string const& project_directory,
                                 std::string const& output_directory)
{
    DBUG("Reading processes:");
    //! \ogs_file_param{prj__processes__process}
    for (auto process_config : processes_config.getConfigSubtreeList("process"))
    {
        auto const type =
            //! \ogs_file_param{prj__processes__process__type}
            process_config.peekConfigParameter<std::string>("type");

        auto const name =
            //! \ogs_file_param{prj__processes__process__name}
            process_config.getConfigParameter<std::string>("name");

        auto const integration_order =
            //! \ogs_file_param{prj__processes__process__integration_order}
            process_config.getConfigParameter<int>("integration_order");

        std::unique_ptr<ProcessLib::Process> process;

        auto jacobian_assembler = ProcessLib::createJacobianAssembler(
            //! \ogs_file_param{prj__processes__process__jacobian_assembler}
            process_config.getConfigSubtreeOptional("jacobian_assembler"));

#ifdef OGS_BUILD_PROCESS_GROUNDWATERFLOW
        if (type == "GROUNDWATER_FLOW")
        {
            // The existence check of the in the configuration referenced
            // process variables is checked in the physical process.
            // TODO at the moment we have only one mesh, later there can be
            // several meshes. Then we have to assign the referenced mesh
            // here.
            process = ProcessLib::GroundwaterFlow::createGroundwaterFlowProcess(
                *_mesh_vec[0], std::move(jacobian_assembler),
                _process_variables, _parameters, integration_order,
                process_config, project_directory, output_directory);
        }
        else
#endif
#ifdef OGS_BUILD_PROCESS_LIQUIDFLOW
            if (type == "LIQUID_FLOW")
        {
            process = ProcessLib::LiquidFlow::createLiquidFlowProcess(
                *_mesh_vec[0], std::move(jacobian_assembler),
                _process_variables, _parameters, integration_order,
                process_config);
        }
        else
#endif
#ifdef OGS_BUILD_PROCESS_TES
            if (type == "TES")
        {
            process = ProcessLib::TES::createTESProcess(
                *_mesh_vec[0], std::move(jacobian_assembler),
                _process_variables, _parameters, integration_order,
                process_config);
        }
        else
#endif
#ifdef OGS_BUILD_PROCESS_HEATCONDUCTION
            if (type == "HEAT_CONDUCTION")
        {
            process = ProcessLib::HeatConduction::createHeatConductionProcess(
                *_mesh_vec[0], std::move(jacobian_assembler),
                _process_variables, _parameters, integration_order,
                process_config);
        }
        else
#endif
#ifdef OGS_BUILD_PROCESS_HYDROMECHANICS
            if (type == "HYDRO_MECHANICS")
        {
            //! \ogs_file_param{prj__processes__process__HYDRO_MECHANICS__dimension}
            switch (process_config.getConfigParameter<int>("dimension"))
            {
                case 2:
                    process =
                        ProcessLib::HydroMechanics::createHydroMechanicsProcess<
                            2>(*_mesh_vec[0], std::move(jacobian_assembler),
                               _process_variables, _parameters,
                               integration_order, process_config);
                    break;
                case 3:
                    process =
                        ProcessLib::HydroMechanics::createHydroMechanicsProcess<
                            3>(*_mesh_vec[0], std::move(jacobian_assembler),
                               _process_variables, _parameters,
                               integration_order, process_config);
                    break;
                default:
                    OGS_FATAL(
                        "HYDRO_MECHANICS process does not support given "
                        "dimension");
            }
        }
        else
#endif
#ifdef OGS_BUILD_PROCESS_LIE
            if (type == "HYDRO_MECHANICS_WITH_LIE")
        {
            //! \ogs_file_param{prj__processes__process__HYDRO_MECHANICS_WITH_LIE__dimension}
            switch (process_config.getConfigParameter<int>("dimension"))
            {
                case 2:
                    process = ProcessLib::LIE::HydroMechanics::
                        createHydroMechanicsProcess<2>(
                            *_mesh_vec[0], std::move(jacobian_assembler),
                            _process_variables, _parameters, integration_order,
                            process_config);
                    break;
                case 3:
                    process = ProcessLib::LIE::HydroMechanics::
                        createHydroMechanicsProcess<3>(
                            *_mesh_vec[0], std::move(jacobian_assembler),
                            _process_variables, _parameters, integration_order,
                            process_config);
                    break;
                default:
                    OGS_FATAL(
                        "HYDRO_MECHANICS_WITH_LIE process does not support "
                        "given dimension");
            }
        }
        else
#endif
#ifdef OGS_BUILD_PROCESS_HT
            if (type == "HT")
        {
            process = ProcessLib::HT::createHTProcess(
                *_mesh_vec[0], std::move(jacobian_assembler),
                _process_variables, _parameters, integration_order,
                process_config);
        }
        else
#endif
#ifdef OGS_BUILD_PROCESS_COMPONENTTRANSPORT
            if (type == "ComponentTransport")
        {
            process =
                ProcessLib::ComponentTransport::createComponentTransportProcess(
                    *_mesh_vec[0], std::move(jacobian_assembler),
                    _process_variables, _parameters, integration_order,
                    process_config);
        }
        else
#endif
#ifdef OGS_BUILD_PROCESS_PHASEFIELD
            if (type == "PHASE_FIELD")
        {
            switch (_mesh_vec[0]->getDimension())
            {
                case 2:
                    process =
                        ProcessLib::PhaseField::createPhaseFieldProcess<
                            2>(*_mesh_vec[0], std::move(jacobian_assembler),
                               _process_variables, _parameters,
                               integration_order, process_config);
                    break;
                case 3:
                    process =
                        ProcessLib::PhaseField::createPhaseFieldProcess<
                            3>(*_mesh_vec[0], std::move(jacobian_assembler),
                               _process_variables, _parameters,
                               integration_order, process_config);
                    break;
            }
        }
        else
#endif
#ifdef OGS_BUILD_PROCESS_RICHARDSCOMPONENTTRANSPORT
            if (type == "RichardsComponentTransport")
        {
            process = ProcessLib::RichardsComponentTransport::
                createRichardsComponentTransportProcess(
                    *_mesh_vec[0], std::move(jacobian_assembler),
                    _process_variables, _parameters, integration_order,
                    process_config);
        }
        else
#endif
#ifdef OGS_BUILD_PROCESS_SMALLDEFORMATION
            if (type == "SMALL_DEFORMATION")
        {
            switch (_mesh_vec[0]->getDimension())
            {
                case 2:
                    process = ProcessLib::SmallDeformation::
                        createSmallDeformationProcess<2>(
                            *_mesh_vec[0], std::move(jacobian_assembler),
                            _process_variables, _parameters, integration_order,
                            process_config);
                    break;
                case 3:
                    process = ProcessLib::SmallDeformation::
                        createSmallDeformationProcess<3>(
                            *_mesh_vec[0], std::move(jacobian_assembler),
                            _process_variables, _parameters, integration_order,
                            process_config);
                    break;
                default:
                    OGS_FATAL(
                        "SMALL_DEFORMATION process does not support "
                        "given dimension");
            }
        }
        else
#endif
#ifdef OGS_BUILD_PROCESS_LIE
            if (type == "SMALL_DEFORMATION_WITH_LIE")
        {
            //! \ogs_file_param{prj__processes__process__SMALL_DEFORMATION_WITH_LIE__dimension}
            switch (process_config.getConfigParameter<int>("dimension"))
            {
                case 2:
                    process = ProcessLib::LIE::SmallDeformation::
                        createSmallDeformationProcess<2>(
                            *_mesh_vec[0], std::move(jacobian_assembler),
                            _process_variables, _parameters, integration_order,
                            process_config);
                    break;
                case 3:
                    process = ProcessLib::LIE::SmallDeformation::
                        createSmallDeformationProcess<3>(
                            *_mesh_vec[0], std::move(jacobian_assembler),
                            _process_variables, _parameters, integration_order,
                            process_config);
                    break;
                default:
                    OGS_FATAL(
                        "SMALL_DEFORMATION_WITH_LIE process does not support "
                        "given dimension");
            }
        }
        else
#endif
#ifdef OGS_BUILD_PROCESS_THERMOMECHANICALPHASEFIELD
            if (type == "THERMO_MECHANICAL_PHASE_FIELD")
        {
            switch (_mesh_vec[0]->getDimension())
            {
                case 2:
                    process =
                        ProcessLib::ThermoMechanicalPhaseField::createThermoMechanicalPhaseFieldProcess<
                            2>(*_mesh_vec[0], std::move(jacobian_assembler),
                               _process_variables, _parameters,
                               integration_order, process_config);
                    break;
                case 3:
                    process =
                        ProcessLib::ThermoMechanicalPhaseField::createThermoMechanicalPhaseFieldProcess<
                            3>(*_mesh_vec[0], std::move(jacobian_assembler),
                               _process_variables, _parameters,
                               integration_order, process_config);
                    break;
            }
        }
        else
#endif
#ifdef OGS_BUILD_PROCESS_THERMOMECHANICS
            if (type == "THERMO_MECHANICS")
        {
            switch (_mesh_vec[0]->getDimension())
            {
                case 2:
                    process = ProcessLib::ThermoMechanics::
                        createThermoMechanicsProcess<2>(
                            *_mesh_vec[0], std::move(jacobian_assembler),
                            _process_variables, _parameters, integration_order,
                            process_config);
                    break;
                case 3:
                    process = ProcessLib::ThermoMechanics::
                        createThermoMechanicsProcess<3>(
                            *_mesh_vec[0], std::move(jacobian_assembler),
                            _process_variables, _parameters, integration_order,
                            process_config);
                    break;
            }
        }
        else
#endif
#ifdef OGS_BUILD_PROCESS_RICHARDSFLOW
            if (type == "RICHARDS_FLOW")
        {
            process = ProcessLib::RichardsFlow::createRichardsFlowProcess(
                *_mesh_vec[0], std::move(jacobian_assembler),
                _process_variables, _parameters, integration_order,
                process_config, _curves);
        }
        else
#endif
#ifdef OGS_BUILD_PROCESS_TWOPHASEFLOWWITHPP
            if (type == "TWOPHASE_FLOW_PP")
        {
            process =
                ProcessLib::TwoPhaseFlowWithPP::createTwoPhaseFlowWithPPProcess(
                    *_mesh_vec[0], std::move(jacobian_assembler),
                    _process_variables, _parameters, integration_order,
                    process_config, _curves);
        }
        else
#endif
#ifdef OGS_BUILD_PROCESS_TWOPHASEFLOWWITHPRHO
            if (type == "TWOPHASE_FLOW_PRHO")
        {
            process = ProcessLib::TwoPhaseFlowWithPrho::
                createTwoPhaseFlowWithPrhoProcess(
                    *_mesh_vec[0], std::move(jacobian_assembler),
                    _process_variables, _parameters, integration_order,
                    process_config, _curves);
        }
        else
#endif
#ifdef OGS_BUILD_PROCESS_THERMALTWOPHASEFLOWWITHPP
            if (type == "THERMAL_TWOPHASE_WITH_PP")
        {
            process = ProcessLib::ThermalTwoPhaseFlowWithPP::
                createThermalTwoPhaseFlowWithPPProcess(
                    *_mesh_vec[0], std::move(jacobian_assembler),
                    _process_variables, _parameters, integration_order,
                    process_config, _curves);
        }
        else
#endif
        {
            OGS_FATAL("Unknown process type: %s", type.c_str());
        }

        BaseLib::insertIfKeyUniqueElseError(_processes,
                                            name,
                                            std::move(process),
                                            "The process name is not unique");
    }
}

void ProjectData::parseTimeLoop(BaseLib::ConfigTree const& config,
                                std::string const& output_directory)
{
    DBUG("Reading time loop configuration.");

    _time_loop = ProcessLib::createUncoupledProcessesTimeLoop(
        config, output_directory, _processes, _nonlinear_solvers);

    if (!_time_loop)
    {
        OGS_FATAL("Initialization of time loop failed.");
    }
}

void ProjectData::parseLinearSolvers(BaseLib::ConfigTree const& config)
{
    DBUG("Reading linear solver configuration.");

    //! \ogs_file_param{prj__linear_solvers__linear_solver}
    for (auto conf : config.getConfigSubtreeList("linear_solver"))
    {
        //! \ogs_file_param{prj__linear_solvers__linear_solver__name}
        auto const name = conf.getConfigParameter<std::string>("name");
        BaseLib::insertIfKeyUniqueElseError(
            _linear_solvers,
            name,
            std::make_unique<GlobalLinearSolver>("", &conf),
            "The linear solver name is not unique");
    }
}

void ProjectData::parseNonlinearSolvers(BaseLib::ConfigTree const& config)
{
    DBUG("Reading linear solver configuration.");

    //! \ogs_file_param{prj__nonlinear_solvers__nonlinear_solver}
    for (auto conf : config.getConfigSubtreeList("nonlinear_solver"))
    {
        auto const ls_name =
            //! \ogs_file_param{prj__nonlinear_solvers__nonlinear_solver__linear_solver}
            conf.getConfigParameter<std::string>("linear_solver");
        auto& linear_solver = BaseLib::getOrError(
            _linear_solvers, ls_name,
            "A linear solver with the given name does not exist.");

        //! \ogs_file_param{prj__nonlinear_solvers__nonlinear_solver__name}
        auto const name = conf.getConfigParameter<std::string>("name");
        BaseLib::insertIfKeyUniqueElseError(
            _nonlinear_solvers,
            name,
            NumLib::createNonlinearSolver(*linear_solver, conf).first,
            "The nonlinear solver name is not unique");
    }
}

void ProjectData::parseCurves(
    boost::optional<BaseLib::ConfigTree> const& config)
{
    if (!config)
        return;

    DBUG("Reading curves configuration.");

    //! \ogs_file_param{prj__curves__curve}
    for (auto conf : config->getConfigSubtreeList("curve"))
    {
        //! \ogs_file_param{prj__curves__curve__name}
        auto const name = conf.getConfigParameter<std::string>("name");
        BaseLib::insertIfKeyUniqueElseError(
            _curves,
            name,
            MathLib::createPiecewiseLinearCurve<
                MathLib::PiecewiseLinearInterpolation>(conf),
            "The curve name is not unique.");
    }
}
