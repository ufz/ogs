/**
 * @copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <tclap/CmdLine.h>

#include "Applications/FileIO/SWMM/SWMMInterface.h"

#include "Applications/ApplicationsLib/LogogSetup.h"

#include "BaseLib/FileTools.h"
#include "BaseLib/StringTools.h"

#include "GeoLib/GEOObjects.h"
#include "GeoLib/IO/XmlIO/Boost/BoostXmlGmlInterface.h"

#include "MeshLib/Mesh.h"
#include "MeshLib/Properties.h"
#include "MeshLib/IO/VtkIO/VtuInterface.h"

#include "Applications/FileIO/CsvInterface.h"

int writeGeoOutput(std::string input_file, std::string output_file)
{
    GeoLib::GEOObjects geo_objects;
    if (!FileIO::SwmmInterface::convertSwmmInputToGeometry(input_file, geo_objects, true))
        return -1;

    GeoLib::IO::BoostXmlGmlInterface xml(geo_objects);
    xml.setNameForExport(BaseLib::extractBaseNameWithoutExtension(input_file));
    xml.writeToFile(output_file);
    return 0;
}

int addObjectsToMesh(
    FileIO::SwmmInterface &swmm,
    MeshLib::Mesh &mesh,
    FileIO::SwmmObject const type,
    std::size_t const timestep)
{
    std::size_t const n_node_params(swmm.getNumberOfParameters(type));
    for (std::size_t j = 0; j<n_node_params; ++j)
    {
        std::string const vec_name(swmm.getArrayName(type, j));
        if (vec_name.empty())
            return -2;
        std::vector<double> data_vec = swmm.getArrayAtTimeStep(type, timestep, j);
        if (!swmm.addResultsToMesh(mesh, type, vec_name, data_vec))
            return -3;
    }
    return 0;
}

int writeMeshOutput(
    std::string const& input_file,
    std::string const& output_file,
    bool const node_args,
    bool const link_args)
{
    std::unique_ptr<FileIO::SwmmInterface> swmm = FileIO::SwmmInterface::create(input_file);
    if (swmm == nullptr)
    return -1;

    MeshLib::Mesh& mesh = swmm->getMesh();

    bool const no_output_file = !swmm->existsSwmmOutputFile();
    if (!(node_args || link_args) || no_output_file)
    {
        if (no_output_file)
            INFO("No output file found.");
        MeshLib::IO::VtuInterface vtkIO(&mesh, 0, false);
        vtkIO.writeToFile(output_file);
        return 0;
    }

    std::string const basename = BaseLib::dropFileExtension(output_file);
    std::string const extension = std::string("." + BaseLib::getFileExtension(output_file));
    std::size_t const n_time_steps(swmm->getNumberOfTimeSteps());
    INFO("Number of simulation time steps: %d", n_time_steps);
    for (std::size_t i=0; i<n_time_steps; i++)
    {
        if (node_args)
            addObjectsToMesh(*swmm, mesh, FileIO::SwmmObject::NODE, i);

        if (link_args)
            addObjectsToMesh(*swmm, mesh, FileIO::SwmmObject::LINK, i);

        MeshLib::IO::VtuInterface vtkio(&mesh, 0, false);
        std::string name(basename + BaseLib::tostring(i) + extension);
        vtkio.writeToFile(name);
    }
    return 0;
}

void writeObjectsOfSwmmTypeToCsv(
    FileIO::SwmmInterface &swmm,
    FileIO::SwmmObject const type,
    std::string const& base,
    std::string const& ext)
{
    std::size_t n_objects = swmm.getNumberOfObjects(type);
    std::string const& type_str (swmm.swmmObjectTypeToString(type));
    for (std::size_t i = 0; i<n_objects; ++i)
    {
        std::string const obj_name = swmm.getName(type, i);
        std::string const obj_file_name = std::string(base + "_" + type_str + "_" + obj_name + ext);
        swmm.writeCsvForObject(obj_file_name, type, i);
    }
}

int writeCsvOutput(
    std::string input_file,
    std::string output_file,
    bool const node_args,
    bool const link_args,
    bool const catchment_args,
    bool const system_args)
{
    std::unique_ptr<FileIO::SwmmInterface> swmm = FileIO::SwmmInterface::create(input_file);
    if (swmm == nullptr)
        return -1;

    if (!swmm->existsSwmmOutputFile())
    {
        INFO("No output file found, skipping data conversion to CSV.");
        return -1;
    }

    if (!(node_args || link_args || catchment_args || system_args))
    {
        INFO("No data category selected. Nothing to write.");
        return 0;
    }

    std::string const basename = BaseLib::dropFileExtension(output_file);
    std::string const extension = std::string("." + BaseLib::getFileExtension(output_file));

    if (node_args)
        writeObjectsOfSwmmTypeToCsv(*swmm, FileIO::SwmmObject::NODE, basename, extension);

    if (link_args)
        writeObjectsOfSwmmTypeToCsv(*swmm, FileIO::SwmmObject::LINK, basename, extension);

    if (catchment_args)
        writeObjectsOfSwmmTypeToCsv(*swmm, FileIO::SwmmObject::SUBCATCHMENT, basename, extension);

    if (system_args)
    {
        std::string const obj_file_name = std::string(basename + "_system" + extension);
        swmm->writeCsvForObject(obj_file_name, FileIO::SwmmObject::SYSTEM, 0);
    }
    return 0;
}

int main(int argc, char *argv[])
{
    ApplicationsLib::LogogSetup setup;

    TCLAP::CmdLine cmd
        ("Read files for the Storm Water Management Model (SWMM) and converts them into OGS data structures.", ' ', "0.1");
    TCLAP::ValueArg<std::string> mesh_output_arg
        ("m","mesh", "mesh output file (*.vtu)", false, "", "mesh output file");
    cmd.add(mesh_output_arg);
    TCLAP::ValueArg<std::string> geo_output_arg
        ("g","geo", "geometry output file (*.gml)", false, "", "geometry output file");
    cmd.add(geo_output_arg);
    TCLAP::ValueArg<std::string> csv_output_arg
        ("c", "csv", "csv output file (*.csv)", false, "", "CSV output file");
    cmd.add(csv_output_arg);
    TCLAP::ValueArg<std::string> swmm_input_arg
        ("i","input", "SWMM input file (*.inp)", true, "", "input file");
    cmd.add(swmm_input_arg);
    TCLAP::SwitchArg add_nodes_arg ("", "node_vars", "Read node variables and add to output mesh");
    cmd.add(add_nodes_arg);
    TCLAP::SwitchArg add_links_arg ("", "link_vars", "Read link variables and add to output mesh");
    cmd.add(add_links_arg);
    TCLAP::SwitchArg add_subcatchments_arg ("", "subcatchment_vars", "Read subcatchment variables and write to CSV-file");
    cmd.add(add_subcatchments_arg);
    TCLAP::SwitchArg add_system_arg ("", "system_vars", "Read system variables and write to CSV-file");
    cmd.add(add_system_arg);
    cmd.parse( argc, argv );

    if (!(geo_output_arg.isSet() || mesh_output_arg.isSet() || csv_output_arg.isSet()))
    {
        ERR ("No output format given. Please specify OGS geometry or mesh output file.");
        return -1;
    }

    if ((add_subcatchments_arg.getValue() || add_system_arg.getValue()) && !csv_output_arg.isSet())
    {
        ERR("Please specify csv output file for exporting subcatchment or system parameters.");
        return -1;
    }

    if (geo_output_arg.isSet())
        writeGeoOutput(swmm_input_arg.getValue(), geo_output_arg.getValue());

    if (mesh_output_arg.isSet())
        writeMeshOutput(swmm_input_arg.getValue(), mesh_output_arg.getValue(),
            add_nodes_arg.getValue(), add_links_arg.getValue());

    if (csv_output_arg.isSet())
        writeCsvOutput(
            swmm_input_arg.getValue(),
            csv_output_arg.getValue(),
            add_nodes_arg.getValue(),
            add_links_arg.getValue(),
            add_subcatchments_arg.getValue(),
            add_system_arg.getValue());

    return 0;
}
