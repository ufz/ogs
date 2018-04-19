/**
 * @copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 *
 */

#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include <tclap/CmdLine.h>

#include <Applications/ApplicationsLib/LogogSetup.h>

#include "BaseLib/BuildInfo.h"
#include "BaseLib/StringTools.h"

#include "MeshLib/Mesh.h"
#include "MeshLib/MeshGenerators/RasterToMesh.h"
#include "MeshLib/IO/VtkIO/VtuInterface.h"


std::string getValue(std::string const& line, std::string const& val_name, bool is_string)
{
    std::string value;
    std::size_t start(line.find(val_name));
    std::size_t end(std::string::npos);
    if (start == end)
    {
        ERR("Value not found.");
        return "";
    }
    value = line.substr(start + 1, std::string::npos);
    if (is_string)
    {
        start = value.find("\"");
        value = value.substr(start + 1, std::string::npos);
        end = value.find("\"");
    }
    else
    {
        start = value.find_first_not_of(" ");
        value = value.substr(start + 1, std::string::npos);
        BaseLib::trim(value);
        end = value.find_first_of(" ");
    }
    value = value.substr(0, end);
    BaseLib::trim(value);
    return value;
}

std::string getName(std::string const& line)
{
    return getValue(line, "T=", true);
}

std::pair<std::size_t, std::size_t> getDimensions(std::string const& line)
{
    std::pair<std::size_t, std::size_t> dims;
    std::stringstream start(getValue(line, "I=", false));
    start >> dims.first;
    std::stringstream end(getValue(line, "J=", false));
    end >> dims.second;
    return dims;
}

std::string trimVariable(std::string& var)
{
    std::size_t const start = var.find_first_not_of("\"");
    var = var.substr(start, std::string::npos);
    std::size_t const end = var.find_first_of("\"");
    return var.substr(0, end);
}

std::vector<std::string> getVariables(std::string const& line)
{
    std::string const var_str ("VARIABLES");
    std::size_t start(line.find(var_str));
    std::string all_vars = line.substr(start+var_str.length(), std::string::npos);
    start = all_vars.find("=");
    all_vars = all_vars.substr(start + 1, std::string::npos);
    BaseLib::trim(all_vars);

    std::vector<std::string> variables;
    std::size_t end = 0;
    while (end != std::string::npos)
    {
        end = all_vars.find_first_of(" ");
        std::string var = all_vars.substr(0, end);
        variables.push_back(trimVariable(var));
        all_vars = all_vars.substr(end+1, std::string::npos);
        BaseLib::trim(all_vars);
    }

    return variables;
}

bool dataCountError(std::string const& name, ::size_t const& current, std::size_t const& total)
{
    if (current != total)
    {
        ERR("Data rows found do not fit specified dimensions for section \"%s\".", name.c_str());
        return true;
    }
    return false;
}

bool dataCountError(std::ofstream& out,
                    std::string const& name,
                    std::size_t const& current,
                    std::size_t const& total)
{
    if (dataCountError(name, current, total))
    {
        out.close();
        return true;
    }
    return false;
}

void resetDataStructures(std::size_t const& n_scalars,
                         std::vector< std::vector<double> >& scalars,
                         std::size_t& val_count)
{
    scalars.clear();
    scalars.reserve(n_scalars);
    for (std::size_t i = 0; i<n_scalars; ++i)
        scalars.push_back(std::vector<double>(0));
    val_count = 0;
}

void writeTecPlotSection(std::ofstream& out,
                         std::string const& file_name,
                         std::size_t& write_count,
                         std::size_t& val_count,
                         std::size_t& val_total)
{
    if (write_count == 0 || val_total != 0)
    {
        std::size_t const delim_pos(file_name.find_last_of("."));
        std::string const base_name(file_name.substr(0, delim_pos+1));
        std::string const extension(file_name.substr(delim_pos, std::string::npos));

        val_count = 0;
        val_total = 0;
        INFO("Writing section #%i", write_count);
        out.close();
        out = std::ofstream(base_name + std::to_string(write_count++) + extension);
    }
}

int writeDataToMesh(std::string const& file_name,
                    std::size_t& write_count,
                    std::vector<std::string> const& vec_names,
                    std::vector< std::vector<double> > const& scalars,
                    std::pair<std::size_t, std::size_t> const& dims)
{
    double cellsize = 0;
    for (std::size_t i=0; i<vec_names.size(); ++i)
    {
        if (vec_names[i] == "x" || vec_names[i] == "X")
        {
            cellsize = scalars[i][1] - scalars[i][0];
            break;
        }
    }

    if (cellsize == 0)
    {
        ERR ("Cell size not found. Aborting...");
        return -4;
    }

    GeoLib::Point origin (0, 0, 0);
    GeoLib::RasterHeader header { dims.first, dims.second, 1, origin, cellsize, -9999 };

    std::unique_ptr<MeshLib::Mesh> mesh (MeshLib::RasterToMesh::convert(
        scalars[0].data(),
        header,
        MeshLib::MeshElemType::QUAD,
        MeshLib::UseIntensityAs::DATAVECTOR,
        vec_names[0]));
    MeshLib::Properties& properties = mesh->getProperties();
    for (std::size_t i=1; i<vec_names.size(); ++i)
    {
        auto* const prop = properties.createNewPropertyVector<double>(
            vec_names[i], MeshLib::MeshItemType::Cell, 1);
        if (!prop)
        {
            ERR("Error creating array \"%s\".", vec_names[i].c_str());
            return -5;
        }
        prop->reserve(scalars[i].size());
        std::copy(scalars[i].cbegin(), scalars[i].cend(), std::back_inserter(*prop));
    }

    std::size_t const delim_pos(file_name.find_last_of("."));
    std::string const base_name(file_name.substr(0, delim_pos + 1));
    std::string const extension(file_name.substr(delim_pos, std::string::npos));

    INFO("Writing section #%i", write_count);
    MeshLib::IO::VtuInterface vtu(mesh.get());
    vtu.writeToFile(base_name + std::to_string(write_count++) + extension);
    return 0;
}

void skipGeometrySection(std::ifstream& in, std::string& line)
{
    while (getline(in, line))
    {
        if ((line.find("TITLE") != std::string::npos) ||
            (line.find("VARIABLES") != std::string::npos) ||
            (line.find("ZONE") != std::string::npos))
            return;
    }
}

int splitFile(std::ifstream& in, std::string file_name)
{
    std::ofstream out;
    std::string line("");
    std::string name;
    std::pair<std::size_t, std::size_t> dims(0,0);
    std::size_t val_count(0), val_total(0);
    std::size_t write_count(0);
    while (getline(in, line))
    {
        if (line.find("TITLE") != std::string::npos)
        {
            if (dataCountError(out, name, val_count, val_total)) return -3;
            writeTecPlotSection(out, file_name, write_count, val_count, val_total);
            out << line << "\n";
            continue;
        }
        else if (line.find("VARIABLES") != std::string::npos)
        {
            if (dataCountError(out, name, val_count, val_total)) return -3;
            writeTecPlotSection(out, file_name, write_count, val_count, val_total);
            out << line << "\n";
            continue;
        }
        else if (line.find("ZONE") != std::string::npos)
        {
            if (dataCountError(out, name, val_count, val_total)) return -3;
            writeTecPlotSection(out, file_name, write_count, val_count, val_total);
            out << line << "\n";
            name = getName(line);
            dims = getDimensions(line);
            val_total = dims.first * dims.second;
            val_count = 0;
            continue;
        }

        out << line << "\n";
        val_count++;
    }
    if (dataCountError(out, name, val_count, val_total))
        return -3;
    INFO("Writing time step #%i", write_count);
    out.close();
    INFO("Finished split.");
    return 0;
}

int convertFile(std::ifstream& in, std::string file_name)
{
    std::string line(""), name("");
    std::pair<std::size_t, std::size_t> dims(0,0);
    std::vector<std::string> var_names;
    std::vector< std::vector<double> > scalars;
    std::size_t val_count(0), val_total(0);
    std::size_t write_count(0);
    while (getline(in, line))
    {
        if (line.find("GEOMETRY") != std::string::npos)
            skipGeometrySection(in, line);

        if (line.empty())
            continue;
        else if (line.find("TITLE") != std::string::npos)
        {
            if (dataCountError(name, val_count, val_total)) return -3;
            if (val_count != 0)
            {
                writeDataToMesh(file_name, write_count, var_names, scalars, dims);
                resetDataStructures(var_names.size(), scalars, val_count);
            }
            continue;
        }
        if (line.find("VARIABLES") != std::string::npos)
        {
            if (val_count != 0)
            {
                if (dataCountError(name, val_count, val_total)) return -3;
                writeDataToMesh(file_name, write_count, var_names, scalars, dims);
            }
            var_names.clear();
            var_names = getVariables(line);
            resetDataStructures(var_names.size(), scalars, val_count);
            continue;
        }
        if (line.find("ZONE") != std::string::npos)
        {
            if (val_count != 0)
            {
                if (dataCountError(name, val_count, val_total)) return -3;
                writeDataToMesh(file_name, write_count, var_names, scalars, dims);
                resetDataStructures(var_names.size(), scalars, val_count);
            }
            name = getName(line);
            dims = getDimensions(line);
            val_total = dims.first * dims.second;
            val_count = 0;
            continue;
        }

        double x;
        std::stringstream iss(line);
        std::size_t i(0);
        std::size_t const n_scalars (scalars.size());
        while (iss >> x)
        {
            if (i > n_scalars-1)
            {
                ERR("Too much data for existing scalar arrays");
                return -3;
            }
            scalars[i++].push_back(x);
        }
        if (i < n_scalars)
        {
            ERR("Not enough data for existing scalar arrays");
            return -3;
        }
        val_count++;
    }
    if (dataCountError(name, val_count, val_total))
        return -3;
    writeDataToMesh(file_name, write_count, var_names, scalars, dims);
    INFO("Finished conversion.");
    return 0;
}

int main(int argc, char *argv[])
{
    ApplicationsLib::LogogSetup logog_setup;

    TCLAP::CmdLine cmd("TecPlot Parser", ' ', BaseLib::BuildInfo::git_describe);
    TCLAP::ValueArg<std::string> input_arg("i", "input-file","TecPlot input file",true,"","string");
    cmd.add( input_arg );
    TCLAP::ValueArg<std::string> output_arg("o", "output-file","output mesh file",false,"","string");
    cmd.add( output_arg );
    TCLAP::SwitchArg split_arg("s","split","split time steps into seperate files");
    cmd.add(split_arg);
    TCLAP::SwitchArg convert_arg("c", "convert","convert TecPlot data into OGS meshes");
    cmd.add(convert_arg);
    cmd.parse( argc, argv );

    if (!input_arg.isSet())
    {
        ERR("No input file given. Please specify TecPlot (*.plt) file");
        return -1;
    }

    if (convert_arg.getValue() && !output_arg.isSet())
    {
        ERR("No output file given. Please specify OGS mesh (*.vtu) file");
        return -1;
    }

    std::ifstream in(input_arg.getValue().c_str());
    if (!in.is_open())
    {
        ERR("Could not open file %s.", input_arg.getValue().c_str());
        return -2;
    }

    if (!convert_arg.isSet() && !split_arg.isSet())
    {
        INFO("Nothing to do. Use -s to split or -c to convert.");
        return 0;
    }

    std::string const filename = (output_arg.isSet()) ? output_arg.getValue() : input_arg.getValue();
    bool const convert (convert_arg.getValue());
    int return_val (0);
    if (split_arg.getValue())
        return_val = splitFile(in, filename);
    else if (convert_arg.getValue())
        return_val = convertFile(in, filename);

    in.close();
    return return_val;
}
