/**
 * \author Karsten Rink
 * \date   2010-08-25
 * \brief  Implementation of the project data class.
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ProjectData.h"

#include <algorithm>
#include <cctype>
#include <set>

#ifdef OGS_USE_PYTHON
#include <pybind11/eval.h>
#endif

#include "BaseLib/Algorithm.h"
#include "BaseLib/ConfigTree.h"
#include "BaseLib/FileTools.h"
#include "BaseLib/Logging.h"
#include "BaseLib/StringTools.h"
#include "GeoLib/GEOObjects.h"
#include "InfoLib/CMakeInfo.h"
#include "MaterialLib/MPL/CreateMedium.h"
#include "MathLib/Curve/CreatePiecewiseLinearCurve.h"
#include "MeshGeoToolsLib/ConstructMeshesFromGeometries.h"
#include "MeshGeoToolsLib/CreateSearchLength.h"
#include "MeshGeoToolsLib/SearchLength.h"
#include "MeshLib/Mesh.h"
#include "NumLib/ODESolver/ConvergenceCriterion.h"
#include "ProcessLib/CreateJacobianAssembler.h"
#include "ProcessLib/DeactivatedSubdomain.h"

// FileIO
#include "GeoLib/IO/XmlIO/Boost/BoostXmlGmlInterface.h"
#include "MeshLib/IO/readMeshFromFile.h"

#include "ParameterLib/ConstantParameter.h"
#include "ParameterLib/Utils.h"
#include "ProcessLib/CreateTimeLoop.h"
#include "ProcessLib/TimeLoop.h"

#ifdef OGS_BUILD_PROCESS_COMPONENTTRANSPORT
#include "ChemistryLib/CreateChemicalSolverInterface.h"
#include "ProcessLib/ComponentTransport/CreateComponentTransportProcess.h"
#endif
#ifdef OGS_BUILD_PROCESS_STEADYSTATEDIFFUSION
#include "ProcessLib/SteadyStateDiffusion/CreateSteadyStateDiffusion.h"
#endif
#ifdef OGS_BUILD_PROCESS_HT
#include "ProcessLib/HT/CreateHTProcess.h"
#endif
#ifdef OGS_BUILD_PROCESS_HEATCONDUCTION
#include "ProcessLib/HeatConduction/CreateHeatConductionProcess.h"
#endif
#ifdef OGS_BUILD_PROCESS_HEATTRANSPORTBHE
#include "ProcessLib/HeatTransportBHE/CreateHeatTransportBHEProcess.h"
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

#ifdef OGS_BUILD_PROCESS_THERMORICHARDSMECHANICS
#include "ProcessLib/ThermoRichardsMechanics/CreateThermoRichardsMechanicsProcess.h"
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
#ifdef OGS_BUILD_PROCESS_RICHARDSMECHANICS
#include "ProcessLib/RichardsMechanics/CreateRichardsMechanicsProcess.h"
#endif
#ifdef OGS_BUILD_PROCESS_SMALLDEFORMATION
#include "ProcessLib/SmallDeformation/CreateSmallDeformationProcess.h"
#endif
#ifdef OGS_BUILD_PROCESS_SMALLDEFORMATIONNONLOCAL
#include "ProcessLib/SmallDeformationNonlocal/CreateSmallDeformationNonlocalProcess.h"
#endif
#ifdef OGS_BUILD_PROCESS_TES
#include "ProcessLib/TES/CreateTESProcess.h"
#endif
#ifdef OGS_BUILD_PROCESS_THERMALTWOPHASEFLOWWITHPP
#include "ProcessLib/ThermalTwoPhaseFlowWithPP/CreateThermalTwoPhaseFlowWithPPProcess.h"
#endif
#ifdef OGS_BUILD_PROCESS_THERMOHYDROMECHANICS
#include "ProcessLib/ThermoHydroMechanics/CreateThermoHydroMechanicsProcess.h"
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

namespace
{
void readGeometry(std::string const& fname, GeoLib::GEOObjects& geo_objects)
{
    DBUG("Reading geometry file '{:s}'.", fname);
    GeoLib::IO::BoostXmlGmlInterface gml_reader(geo_objects);
    gml_reader.readFile(fname);
}

std::unique_ptr<MeshLib::Mesh> readSingleMesh(
    BaseLib::ConfigTree const& mesh_config_parameter,
    std::string const& project_directory)
{
    std::string const mesh_file = BaseLib::copyPathToFileName(
        mesh_config_parameter.getValue<std::string>(), project_directory);
    DBUG("Reading mesh file '{:s}'.", mesh_file);

    auto mesh = std::unique_ptr<MeshLib::Mesh>(
        MeshLib::IO::readMeshFromFile(mesh_file));
    if (!mesh)
    {
        OGS_FATAL("Could not read mesh from '{:s}' file. No mesh added.",
                  mesh_file);
    }

#ifdef DOXYGEN_DOCU_ONLY
    //! \ogs_file_attr{prj__meshes__mesh__axially_symmetric}
    mesh_config_parameter.getConfigAttributeOptional<bool>("axially_symmetric");
#endif  // DOXYGEN_DOCU_ONLY

    if (auto const axially_symmetric =
            //! \ogs_file_attr{prj__mesh__axially_symmetric}
        mesh_config_parameter.getConfigAttributeOptional<bool>(
            "axially_symmetric"))
    {
        mesh->setAxiallySymmetric(*axially_symmetric);
    }

    return mesh;
}

std::vector<std::unique_ptr<MeshLib::Mesh>> readMeshes(
    BaseLib::ConfigTree const& config, std::string const& project_directory)
{
    std::vector<std::unique_ptr<MeshLib::Mesh>> meshes;

    //! \ogs_file_param{prj__meshes}
    auto optional_meshes = config.getConfigSubtreeOptional("meshes");
    if (optional_meshes)
    {
        DBUG("Reading multiple meshes.");
        //! \ogs_file_param{prj__meshes__mesh}
        auto const configs = optional_meshes->getConfigParameterList("mesh");
        std::transform(
            configs.begin(), configs.end(), std::back_inserter(meshes),
            [&project_directory](auto const& mesh_config) {
                return readSingleMesh(mesh_config, project_directory);
            });
    }
    else
    {  // Read single mesh with geometry.
        WARN(
            "Consider switching from mesh and geometry input to multiple "
            "meshes input. See "
            "https://www.opengeosys.org/docs/tools/model-preparation/"
            "constructmeshesfromgeometry/ tool for conversion.");
        //! \ogs_file_param{prj__mesh}
        meshes.push_back(readSingleMesh(config.getConfigParameter("mesh"),
                                        project_directory));

        std::string const geometry_file = BaseLib::copyPathToFileName(
            //! \ogs_file_param{prj__geometry}
            config.getConfigParameter<std::string>("geometry"),
            project_directory);
        GeoLib::GEOObjects geoObjects;
        readGeometry(geometry_file, geoObjects);

        std::unique_ptr<MeshGeoToolsLib::SearchLength> search_length_algorithm =
            MeshGeoToolsLib::createSearchLengthAlgorithm(config, *meshes[0]);
        bool const multiple_nodes_allowed = false;
        auto additional_meshes =
            MeshGeoToolsLib::constructAdditionalMeshesFromGeoObjects(
                geoObjects, *meshes[0], std::move(search_length_algorithm),
                multiple_nodes_allowed);

        std::move(begin(additional_meshes), end(additional_meshes),
                  std::back_inserter(meshes));
    }
    return meshes;
}

std::optional<ParameterLib::CoordinateSystem> parseLocalCoordinateSystem(
    std::optional<BaseLib::ConfigTree> const& config,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters)
{
    if (!config)
    {
        return {};
    }

    DBUG("Reading coordinate system configuration.");

    //
    // Fetch the first basis vector; its length defines the dimension.
    //
    auto const& basis_vector_0 = ParameterLib::findParameter<double>(
        *config,
        //! \ogs_file_param_special{prj__local_coordinate_system__basis_vector_0}
        "basis_vector_0", parameters, 0 /* any dimension */);
    int const dimension = basis_vector_0.getNumberOfGlobalComponents();

    // check dimension
    if (dimension != 2 && dimension != 3)
    {
        OGS_FATAL(
            "Basis vector parameter '{:s}' must have two or three components, "
            "but it has {:d}.",
            basis_vector_0.name, dimension);
    }

    //
    // Fetch the second basis vector, which must be of the same dimension as the
    // first one.
    //
    auto const& basis_vector_1 = ParameterLib::findParameter<double>(
        *config,
        //! \ogs_file_param_special{prj__local_coordinate_system__basis_vector_1}
        "basis_vector_1", parameters, dimension);

    //
    // For two dimensions, we are done; construct coordinate system;
    //
    if (dimension == 2)
    {
        return ParameterLib::CoordinateSystem{basis_vector_0, basis_vector_1};
    }

    //
    // Parse the third vector, for three dimensions.
    //
    auto const& basis_vector_2 = ParameterLib::findParameter<double>(
        *config,
        //! \ogs_file_param_special{prj__local_coordinate_system__basis_vector_2}
        "basis_vector_2", parameters, dimension);
    return ParameterLib::CoordinateSystem{basis_vector_0, basis_vector_1,
                                          basis_vector_2};
}
}  // namespace

ProjectData::ProjectData() = default;

ProjectData::ProjectData(BaseLib::ConfigTree const& project_config,
                         std::string const& project_directory,
                         std::string const& output_directory)
    : _mesh_vec(readMeshes(project_config, project_directory))
{
    if (auto const python_script =
            //! \ogs_file_param{prj__python_script}
        project_config.getConfigParameterOptional<std::string>("python_script"))
    {
#ifdef OGS_USE_PYTHON
        namespace py = pybind11;

        // Append to python's module search path
        auto py_path = py::module::import("sys").attr("path");
        py_path.attr("append")(project_directory);  // .prj directory
        // virtualenv
        py_path.attr("append")(
            CMakeInfoLib::CMakeInfo::python_virtualenv_sitepackages);

        auto const script_path =
            BaseLib::copyPathToFileName(*python_script, project_directory);

        // Evaluate in scope of main module
        py::object scope = py::module::import("__main__").attr("__dict__");
        // add (global) variables
        auto globals = py::dict(scope);
        globals["ogs_prj_directory"] = project_directory;
        py::eval_file(script_path, scope);
#else
        OGS_FATAL("OpenGeoSys has not been built with Python support.");
#endif  // OGS_USE_PYTHON
    }

    //! \ogs_file_param{prj__curves}
    parseCurves(project_config.getConfigSubtreeOptional("curves"));

    auto parameter_names_for_transformation =
        //! \ogs_file_param{prj__parameters}
        parseParameters(project_config.getConfigSubtree("parameters"));

    _local_coordinate_system = parseLocalCoordinateSystem(
        //! \ogs_file_param{prj__local_coordinate_system}
        project_config.getConfigSubtreeOptional("local_coordinate_system"),
        _parameters);

    for (auto& parameter : _parameters)
    {
        if (std::find(begin(parameter_names_for_transformation),
                      end(parameter_names_for_transformation),
                      parameter->name) !=
            end(parameter_names_for_transformation))
        {
            if (!_local_coordinate_system)
            {
                OGS_FATAL(
                    "The parameter '{:s}' is using the local coordinate system "
                    "but no local coordinate system was provided.",
                    parameter->name);
            }
            parameter->setCoordinateSystem(*_local_coordinate_system);
        }

        parameter->initialize(_parameters);
    }

    //! \ogs_file_param{prj__process_variables}
    parseProcessVariables(project_config.getConfigSubtree("process_variables"));

    //! \ogs_file_param{prj__media}
    parseMedia(project_config.getConfigSubtreeOptional("media"));

    auto chemical_solver_interface = parseChemicalSolverInterface(
        //! \ogs_file_param{prj__chemical_system}
        project_config.getConfigSubtreeOptional("chemical_system"),
        output_directory);

    //! \ogs_file_param{prj__processes}
    parseProcesses(project_config.getConfigSubtree("processes"),
                   project_directory, output_directory,
                   std::move(chemical_solver_interface));

    //! \ogs_file_param{prj__linear_solvers}
    parseLinearSolvers(project_config.getConfigSubtree("linear_solvers"));

    //! \ogs_file_param{prj__nonlinear_solvers}
    parseNonlinearSolvers(project_config.getConfigSubtree("nonlinear_solvers"));

    //! \ogs_file_param{prj__time_loop}
    parseTimeLoop(project_config.getConfigSubtree("time_loop"),
                  output_directory);
}

void ProjectData::parseProcessVariables(
    BaseLib::ConfigTree const& process_variables_config)
{
    DBUG("Parse process variables:");

    std::set<std::string> names;

    for (auto var_config
         //! \ogs_file_param{prj__process_variables__process_variable}
         : process_variables_config.getConfigSubtreeList("process_variable"))
    {
        // Either the mesh name is given, or the first mesh's name will be
        // taken. Taking the first mesh's value is deprecated.
        auto const mesh_name =
            //! \ogs_file_param{prj__process_variables__process_variable__mesh}
            var_config.getConfigParameter<std::string>("mesh",
                                                       _mesh_vec[0]->getName());

        auto& mesh = *BaseLib::findElementOrError(
            begin(_mesh_vec), end(_mesh_vec),
            [&mesh_name](auto const& m) { return m->getName() == mesh_name; },
            "Expected to find a mesh named " + mesh_name + ".");

        auto pv = ProcessLib::ProcessVariable{var_config, mesh, _mesh_vec,
                                              _parameters, _curves};
        if (!names.insert(pv.getName()).second)
        {
            OGS_FATAL("A process variable with name `{:s}' already exists.",
                      pv.getName());
        }

        _process_variables.push_back(std::move(pv));
    }
}

std::vector<std::string> ProjectData::parseParameters(
    BaseLib::ConfigTree const& parameters_config)
{
    using namespace ProcessLib;

    std::set<std::string> names;
    std::vector<std::string> parameter_names_for_transformation;

    DBUG("Reading parameters:");
    for (auto parameter_config :
         //! \ogs_file_param{prj__parameters__parameter}
         parameters_config.getConfigSubtreeList("parameter"))
    {
        auto p =
            ParameterLib::createParameter(parameter_config, _mesh_vec, _curves);
        if (!names.insert(p->name).second)
        {
            OGS_FATAL("A parameter with name `{:s}' already exists.", p->name);
        }

        auto const use_local_coordinate_system =
            //! \ogs_file_param{prj__parameters__parameter__use_local_coordinate_system}
            parameter_config.getConfigParameterOptional<bool>(
                "use_local_coordinate_system");
        if (!!use_local_coordinate_system && *use_local_coordinate_system)
        {
            parameter_names_for_transformation.push_back(p->name);
        }

        _parameters.push_back(std::move(p));
    }

    _parameters.push_back(
        std::make_unique<ParameterLib::ConstantParameter<double>>(
            ProcessLib::DeactivatedSubdomain::zero_parameter_name, 0.0));

    return parameter_names_for_transformation;
}

void ProjectData::parseMedia(
    std::optional<BaseLib::ConfigTree> const& media_config)
{
    if (!media_config)
    {
        return;
    }

    DBUG("Reading media:");

    if (_mesh_vec.empty() || _mesh_vec[0] == nullptr)
    {
        ERR("A mesh is required to define medium materials.");
        return;
    }

    for (auto const& medium_config :
         //! \ogs_file_param{prj__media__medium}
         media_config->getConfigSubtreeList("medium"))
    {
        auto material_id_string =
            //! \ogs_file_attr{prj__media__medium__id}
            medium_config.getConfigAttribute<std::string>("id", "0");
        material_id_string.erase(
            remove_if(begin(material_id_string), end(material_id_string),
                      [](unsigned char const c) { return std::isspace(c); }),
            end(material_id_string));
        auto const material_ids_strings =
            BaseLib::splitString(material_id_string, ',');

        // Convert strings to ints;
        std::vector<int> material_ids;
        std::transform(
            begin(material_ids_strings), end(material_ids_strings),
            std::back_inserter(material_ids), [](std::string const& m_id) {
                if (auto const it = std::find_if_not(
                        begin(m_id), end(m_id),
                        [](unsigned char const c) { return std::isdigit(c); });
                    it != end(m_id))
                {
                    OGS_FATAL(
                        "Could not parse material ID's from '{:s}'. Please "
                        "separate multiple material ID's by comma only. "
                        "Invalid character: '%c'",
                        m_id, *it);
                }
                return std::stoi(m_id);
            });

        for (auto const& id : material_ids)
        {
            if (_media.find(id) != end(_media))
            {
                OGS_FATAL(
                    "Multiple media were specified for the same material id "
                    "'{:d}'. "
                    "Keep in mind, that if no material id is specified, it is "
                    "assumed to be 0 by default.",
                    id);
            }

            _media[id] =
                (id == material_ids[0])
                    ? MaterialPropertyLib::createMedium(
                          _mesh_vec[0]->getDimension(), medium_config,
                          _parameters,
                          _local_coordinate_system ? &*_local_coordinate_system
                                                   : nullptr,
                          _curves)
                    : _media[material_ids[0]];
        }
    }

    if (_media.empty())
    {
        OGS_FATAL("No entity is found inside <media>.");
    }
}

std::unique_ptr<ChemistryLib::ChemicalSolverInterface>
ProjectData::parseChemicalSolverInterface(
    std::optional<BaseLib::ConfigTree> const& config,
    std::string const& output_directory)
{
    if (!config)
    {
        return nullptr;
    }

    std::unique_ptr<ChemistryLib::ChemicalSolverInterface>
        chemical_solver_interface;
#ifdef OGS_BUILD_PROCESS_COMPONENTTRANSPORT
    INFO(
        "Ready for initializing interface to a chemical solver for water "
        "chemistry calculation.");

    auto const chemical_solver =
        //! \ogs_file_attr{prj__chemical_system__chemical_solver}
        config->getConfigAttribute<std::string>("chemical_solver");

    if (boost::iequals(chemical_solver, "Phreeqc"))
    {
        INFO(
            "Configuring phreeqc interface for water chemistry "
            "calculation using file-based approach.");

        chemical_solver_interface = ChemistryLib::createChemicalSolverInterface<
            ChemistryLib::ChemicalSolver::Phreeqc>(_mesh_vec, *config,
                                                   output_directory);
    }
    else if (boost::iequals(chemical_solver, "PhreeqcKernel"))
    {
        OGS_FATAL(
            "The chemical solver option of PhreeqcKernel is not accessible "
            "for the time being. Please set 'Phreeqc'' as the chemical "
            "solver for reactive transport modeling.");
    }
    else
    {
        OGS_FATAL(
            "Unknown chemical solver. Please specify either Phreeqc or "
            "PhreeqcKernel as the solver for water chemistry calculation "
            "instead.");
    }
#else
    (void)output_directory;

    OGS_FATAL(
        "Found the type of the process to be solved is not component transport "
        "process. Please specify the process type to ComponentTransport. At "
        "the present, water chemistry calculation is only available for "
        "component transport process.");
#endif
    return chemical_solver_interface;
}

void ProjectData::parseProcesses(
    BaseLib::ConfigTree const& processes_config,
    std::string const& project_directory,
    std::string const& output_directory,
    [[maybe_unused]] std::unique_ptr<ChemistryLib::ChemicalSolverInterface>&&
        chemical_solver_interface)
{
    (void)project_directory;  // to avoid compilation warning
    (void)output_directory;   // to avoid compilation warning

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

#ifdef OGS_BUILD_PROCESS_STEADYSTATEDIFFUSION
        if (type == "STEADY_STATE_DIFFUSION")
        {
            // The existence check of the in the configuration referenced
            // process variables is checked in the physical process.
            // TODO at the moment we have only one mesh, later there can be
            // several meshes. Then we have to assign the referenced mesh
            // here.
            process =
                ProcessLib::SteadyStateDiffusion::createSteadyStateDiffusion(
                    name, *_mesh_vec[0], std::move(jacobian_assembler),
                    _process_variables, _parameters, integration_order,
                    process_config, _mesh_vec, _media);
        }
        else
#endif
#ifdef OGS_BUILD_PROCESS_LIQUIDFLOW
            if (type == "LIQUID_FLOW")
        {
            process = ProcessLib::LiquidFlow::createLiquidFlowProcess(
                name, *_mesh_vec[0], std::move(jacobian_assembler),
                _process_variables, _parameters, integration_order,
                process_config, _mesh_vec, _media);
        }
        else
#endif
#ifdef OGS_BUILD_PROCESS_TES
            if (type == "TES")
        {
            process = ProcessLib::TES::createTESProcess(
                name, *_mesh_vec[0], std::move(jacobian_assembler),
                _process_variables, _parameters, integration_order,
                process_config);
        }
        else
#endif
#ifdef OGS_BUILD_PROCESS_HEATCONDUCTION
            if (type == "HEAT_CONDUCTION")
        {
            process = ProcessLib::HeatConduction::createHeatConductionProcess(
                name, *_mesh_vec[0], std::move(jacobian_assembler),
                _process_variables, _parameters, integration_order,
                process_config, _media);
        }
        else
#endif
#ifdef OGS_BUILD_PROCESS_HEATTRANSPORTBHE
            if (type == "HEAT_TRANSPORT_BHE")
        {
            if (_mesh_vec[0]->getDimension() != 3)
            {
                OGS_FATAL(
                    "HEAT_TRANSPORT_BHE can only work with a 3-dimensional "
                    "mesh! ");
            }

            process =
                ProcessLib::HeatTransportBHE::createHeatTransportBHEProcess(
                    name, *_mesh_vec[0], std::move(jacobian_assembler),
                    _process_variables, _parameters, integration_order,
                    process_config, _curves, _media);
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
                            2>(name, *_mesh_vec[0],
                               std::move(jacobian_assembler),
                               _process_variables, _parameters,
                               _local_coordinate_system, integration_order,
                               process_config, _media);
                    break;
                case 3:
                    process =
                        ProcessLib::HydroMechanics::createHydroMechanicsProcess<
                            3>(name, *_mesh_vec[0],
                               std::move(jacobian_assembler),
                               _process_variables, _parameters,
                               _local_coordinate_system, integration_order,
                               process_config, _media);
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
                            name, *_mesh_vec[0], std::move(jacobian_assembler),
                            _process_variables, _parameters,
                            _local_coordinate_system, integration_order,
                            process_config);
                    break;
                case 3:
                    process = ProcessLib::LIE::HydroMechanics::
                        createHydroMechanicsProcess<3>(
                            name, *_mesh_vec[0], std::move(jacobian_assembler),
                            _process_variables, _parameters,
                            _local_coordinate_system, integration_order,
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
                name, *_mesh_vec[0], std::move(jacobian_assembler),
                _process_variables, _parameters, integration_order,
                process_config, _mesh_vec, _media);
        }
        else
#endif
#ifdef OGS_BUILD_PROCESS_COMPONENTTRANSPORT
            if (type == "ComponentTransport")
        {
            process =
                ProcessLib::ComponentTransport::createComponentTransportProcess(
                    name, *_mesh_vec[0], std::move(jacobian_assembler),
                    _process_variables, _parameters, integration_order,
                    process_config, _mesh_vec, _media,
                    std::move(chemical_solver_interface));
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
                        ProcessLib::PhaseField::createPhaseFieldProcess<2>(
                            name, *_mesh_vec[0], std::move(jacobian_assembler),
                            _process_variables, _parameters,
                            _local_coordinate_system, integration_order,
                            process_config);
                    break;
                case 3:
                    process =
                        ProcessLib::PhaseField::createPhaseFieldProcess<3>(
                            name, *_mesh_vec[0], std::move(jacobian_assembler),
                            _process_variables, _parameters,
                            _local_coordinate_system, integration_order,
                            process_config);
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
                    name, *_mesh_vec[0], std::move(jacobian_assembler),
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
                            name, *_mesh_vec[0], std::move(jacobian_assembler),
                            _process_variables, _parameters,
                            _local_coordinate_system, integration_order,
                            process_config);
                    break;
                case 3:
                    process = ProcessLib::SmallDeformation::
                        createSmallDeformationProcess<3>(
                            name, *_mesh_vec[0], std::move(jacobian_assembler),
                            _process_variables, _parameters,
                            _local_coordinate_system, integration_order,
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
#ifdef OGS_BUILD_PROCESS_SMALLDEFORMATIONNONLOCAL
            if (type == "SMALL_DEFORMATION_NONLOCAL")
        {
            switch (_mesh_vec[0]->getDimension())
            {
                case 2:
                    process = ProcessLib::SmallDeformationNonlocal::
                        createSmallDeformationNonlocalProcess<2>(
                            name, *_mesh_vec[0], std::move(jacobian_assembler),
                            _process_variables, _parameters,
                            _local_coordinate_system, integration_order,
                            process_config);
                    break;
                case 3:
                    process = ProcessLib::SmallDeformationNonlocal::
                        createSmallDeformationNonlocalProcess<3>(
                            name, *_mesh_vec[0], std::move(jacobian_assembler),
                            _process_variables, _parameters,
                            _local_coordinate_system, integration_order,
                            process_config);
                    break;
                default:
                    OGS_FATAL(
                        "SMALL_DEFORMATION_NONLOCAL process does not support "
                        "given dimension {:d}",
                        _mesh_vec[0]->getDimension());
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
                            name, *_mesh_vec[0], std::move(jacobian_assembler),
                            _process_variables, _parameters,
                            _local_coordinate_system, integration_order,
                            process_config);
                    break;
                case 3:
                    process = ProcessLib::LIE::SmallDeformation::
                        createSmallDeformationProcess<3>(
                            name, *_mesh_vec[0], std::move(jacobian_assembler),
                            _process_variables, _parameters,
                            _local_coordinate_system, integration_order,
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
#ifdef OGS_BUILD_PROCESS_THERMOHYDROMECHANICS
            if (type == "THERMO_HYDRO_MECHANICS")
        {
            //! \ogs_file_param{prj__processes__process__THERMO_HYDRO_MECHANICS__dimension}
            switch (process_config.getConfigParameter<int>("dimension"))
            {
                case 2:
                    process = ProcessLib::ThermoHydroMechanics::
                        createThermoHydroMechanicsProcess<2>(
                            name, *_mesh_vec[0], std::move(jacobian_assembler),
                            _process_variables, _parameters,
                            _local_coordinate_system, integration_order,
                            process_config, _media);
                    break;
                case 3:
                    process = ProcessLib::ThermoHydroMechanics::
                        createThermoHydroMechanicsProcess<3>(
                            name, *_mesh_vec[0], std::move(jacobian_assembler),
                            _process_variables, _parameters,
                            _local_coordinate_system, integration_order,
                            process_config, _media);
                    break;
                default:
                    OGS_FATAL(
                        "THERMO_HYDRO_MECHANICS process does not support given "
                        "dimension");
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
                    process = ProcessLib::ThermoMechanicalPhaseField::
                        createThermoMechanicalPhaseFieldProcess<2>(
                            name, *_mesh_vec[0], std::move(jacobian_assembler),
                            _process_variables, _parameters,
                            _local_coordinate_system, integration_order,
                            process_config);
                    break;
                case 3:
                    process = ProcessLib::ThermoMechanicalPhaseField::
                        createThermoMechanicalPhaseFieldProcess<3>(
                            name, *_mesh_vec[0], std::move(jacobian_assembler),
                            _process_variables, _parameters,
                            _local_coordinate_system, integration_order,
                            process_config);
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
                            name, *_mesh_vec[0], std::move(jacobian_assembler),
                            _process_variables, _parameters,
                            _local_coordinate_system, integration_order,
                            process_config, _media);
                    break;
                case 3:
                    process = ProcessLib::ThermoMechanics::
                        createThermoMechanicsProcess<3>(
                            name, *_mesh_vec[0], std::move(jacobian_assembler),
                            _process_variables, _parameters,
                            _local_coordinate_system, integration_order,
                            process_config, _media);
                    break;
            }
        }
        else
#endif
#ifdef OGS_BUILD_PROCESS_RICHARDSFLOW
            if (type == "RICHARDS_FLOW")
        {
            process = ProcessLib::RichardsFlow::createRichardsFlowProcess(
                name, *_mesh_vec[0], std::move(jacobian_assembler),
                _process_variables, _parameters, integration_order,
                process_config, _curves, _media);
        }
        else
#endif
#ifdef OGS_BUILD_PROCESS_RICHARDSMECHANICS
            if (type == "RICHARDS_MECHANICS")
        {
            //! \ogs_file_param{prj__processes__process__RICHARDS_MECHANICS__dimension}
            switch (process_config.getConfigParameter<int>("dimension"))
            {
                case 2:
                    process = ProcessLib::RichardsMechanics::
                        createRichardsMechanicsProcess<2>(
                            name, *_mesh_vec[0], std::move(jacobian_assembler),
                            _process_variables, _parameters,
                            _local_coordinate_system, integration_order,
                            process_config, _media);
                    break;
                case 3:
                    process = ProcessLib::RichardsMechanics::
                        createRichardsMechanicsProcess<3>(
                            name, *_mesh_vec[0], std::move(jacobian_assembler),
                            _process_variables, _parameters,
                            _local_coordinate_system, integration_order,
                            process_config, _media);
                    break;
            }
        }
        else
#endif

#ifdef OGS_BUILD_PROCESS_THERMORICHARDSMECHANICS
            if (type == "THERMO_RICHARDS_MECHANICS")
        {
            switch (_mesh_vec[0]->getDimension())
            {
                case 2:
                    process = ProcessLib::ThermoRichardsMechanics::
                        createThermoRichardsMechanicsProcess<2>(
                            name, *_mesh_vec[0], std::move(jacobian_assembler),
                            _process_variables, _parameters,
                            _local_coordinate_system, integration_order,
                            process_config, _media);
                    break;
                case 3:
                    process = ProcessLib::ThermoRichardsMechanics::
                        createThermoRichardsMechanicsProcess<3>(
                            name, *_mesh_vec[0], std::move(jacobian_assembler),
                            _process_variables, _parameters,
                            _local_coordinate_system, integration_order,
                            process_config, _media);
                    break;
            }
        }
        else
#endif

#ifdef OGS_BUILD_PROCESS_TWOPHASEFLOWWITHPP
            if (type == "TWOPHASE_FLOW_PP")
        {
            process =
                ProcessLib::TwoPhaseFlowWithPP::createTwoPhaseFlowWithPPProcess(
                    name, *_mesh_vec[0], std::move(jacobian_assembler),
                    _process_variables, _parameters, integration_order,
                    process_config, _curves, _media);
        }
        else
#endif
#ifdef OGS_BUILD_PROCESS_TWOPHASEFLOWWITHPRHO
            if (type == "TWOPHASE_FLOW_PRHO")
        {
            process = ProcessLib::TwoPhaseFlowWithPrho::
                createTwoPhaseFlowWithPrhoProcess(
                    name, *_mesh_vec[0], std::move(jacobian_assembler),
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
                    name, *_mesh_vec[0], std::move(jacobian_assembler),
                    _process_variables, _parameters, integration_order,
                    process_config, _curves);
        }
        else
#endif
        {
            OGS_FATAL("Unknown process type: {:s}", type);
        }

        if (BaseLib::containsIf(
                _processes,
                [&name](std::unique_ptr<ProcessLib::Process> const& p) {
                    return p->name == name;
                }))
        {
            OGS_FATAL("The process name '{:s}' is not unique.", name);
        }
        _processes.push_back(std::move(process));
    }
}

void ProjectData::parseTimeLoop(BaseLib::ConfigTree const& config,
                                std::string const& output_directory)
{
    DBUG("Reading time loop configuration.");

    _time_loop = ProcessLib::createTimeLoop(
        config, output_directory, _processes, _nonlinear_solvers, _mesh_vec);

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

void ProjectData::parseCurves(std::optional<BaseLib::ConfigTree> const& config)
{
    if (!config)
    {
        return;
    }

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
