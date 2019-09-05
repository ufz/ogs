/**
 * @copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

// STL
#include <cctype>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <utility>

#include <tclap/CmdLine.h>

#include "Applications/ApplicationsLib/LogogSetup.h"

#include <netcdf>

// BaseLib
#include "BaseLib/FileTools.h"
#include "BaseLib/LogogSimpleFormatter.h"
#include "InfoLib/GitInfo.h"

// GeoLib
#include "GeoLib/Raster.h"

// MeshLib
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshGenerators/RasterToMesh.h"

#include "MeshLib/IO/VtkIO/VtuInterface.h"

using namespace netCDF;

void checkExit(std::string const& str)
{
    if (str == "x" || str == "exit")
        exit(0);
}

void showErrorMessage(std::size_t error_id, std::size_t max = 0)
{
    if (error_id == 0)
    {
        ERR ("Input not valid.");
    }
    else if (error_id == 1)
    {
        ERR ("Index not valid. Valid indices are in [0,%d].", max);
    }
    else if (error_id == 2)
    {
        showErrorMessage(0);
        std::cout << "Type \"info\" to display the available options again. \"exit\" will exit the programme.\n";
    }
}

std::size_t parseInput(std::string const& request_str, std::size_t max_val, bool has_info = false)
{
    while (true)
    {
        std::cout << request_str;
        std::string str;
        std::size_t val;
        std::getline(std::cin, str);
        checkExit(str);
        if (has_info && str == "info")
            return (max_val);
        std::stringstream str_stream(str);
        if (!(str_stream >> val))
        {
            std::size_t error_val = (has_info) ? 2 : 0;
            showErrorMessage(error_val);
            continue;
        }
        if (val > max_val - 1)
        {
            showErrorMessage(1, max_val - 1);
            continue;
        }
        return val;
    }
    return std::numeric_limits<std::size_t>::max();
}

NcVar getDimVar(NcFile const& dataset, NcVar const& var, std::size_t const dim)
{
    NcDim const& dim_obj = var.getDim(dim);
    return dataset.getVar(dim_obj.getName());
}

std::vector<std::string> showArrays(NcFile const& dataset)
{
    std::size_t const n_vars(dataset.getDimCount());
    std::cout << "The NetCDF file contains the following " << n_vars << " arrays:\n\n";
    std::cout << "\tIndex\tArray Name\t#Dimensions\n";
    std::cout << "-------------------------------------------\n";
    auto const& names = dataset.getVars();
    std::vector<std::string> var_names;
    for (auto [name, var] : names)
    {
        std::cout << "\t" << var_names.size() << "\t" << name << "\t(" << var.getDimCount() << "D array)\n";
        var_names.push_back(name);
    }
    std::cout << "\n";
    return var_names;
}

bool showArraysDims(NcVar const& var)
{
    std::cout << "Data array \"" << var.getName()
              << "\" contains the following dimensions:\n";
    std::size_t const n_dims (var.getDimCount());
    for (std::size_t i=0; i<n_dims; ++i)
        std::cout << "\t" << i << "\t" << var.getDim(i).getName() << "\t(" << var.getDim(i).getSize() << " values)\n";
    std::cout << "\n";
    return true;
}

std::pair<double, double> getBoundaries(NcVar const& var)
{
    if ((var.getDimCount()) == 1)
    {
        double start, end;
        std::size_t size = var.getDim(0).getSize();
        var.getVar({0}, {1}, &start);
        var.getVar({size - 1}, {1}, &end);
        return std::make_pair(start, end);
    }
    return std::make_pair(0, 0);
}

MathLib::Point3d getOrigin(NcFile const& dataset, NcVar const& var)
{
    MathLib::Point3d origin(std::array<double, 3>{ {0, 0, 0}});
    std::size_t const n_dims = var.getDimCount();
    for (std::size_t i=1; i<n_dims; ++i)
    {
        NcVar const& dim = getDimVar(dataset, var, var.getDimCount() - i);
        auto const bounds = getBoundaries(dim);
        origin[n_dims - i] = (bounds.first < bounds.second) ? bounds.first : bounds.second;
    }
    return origin;
}

void reverseNorthSouth(double* data, std::size_t width, std::size_t height)
{
    std::size_t const total_length (width * height);
    double* cp_array = new double[total_length];

    for (std::size_t i=0; i<height; i++)
    {
        for (std::size_t j=0; j<width; j++)
        {
            std::size_t const old_index (total_length - (width*(i+1)));
            std::size_t const new_index (width*i);
            cp_array[new_index+j] = data[old_index+j];
        }
    }

    for (std::size_t i=0; i<total_length; i++)
        data[i] = cp_array[i];

    delete[] cp_array;
}

bool canConvert(NcVar const& var)
{
    bool ret (true);
    if (ret = (var.getDimCount() < 2))
        ERR ("Only 2+ dimensional variables can be converted into OGS Meshes.\n");
    return !ret;
}

std::string arraySelectionLoop(NcFile const& dataset)
{
    std::vector<std::string> const& names = showArrays(dataset);
    std::size_t idx = parseInput("Enter data array index: ", dataset.getVarCount(), true);

    if (idx == dataset.getVarCount() || !canConvert(dataset.getVar(names[idx])))
        return arraySelectionLoop(dataset);

    return names[idx];
}

bool dimensionSelectionLoop(NcVar const& var, std::array<std::size_t,4> &dim_idx_map)
{
    showArraysDims(var);
    std::size_t const n_dims(var.getDimCount());
    dim_idx_map[0] = std::numeric_limits<std::size_t>::max();
    bool is_time_dep (true);

    // get temporal dimension
    if (n_dims > 1)
    {
        std::string temp_str("");
        cout << "Is the parameter time-dependent?\n";
        while (dim_idx_map[0] == std::numeric_limits<std::size_t>::max() && is_time_dep == true)
        {
            cout << "Enter ID for temporal dimension or \"c\" to continue: ";
            std::getline(std::cin, temp_str);
            std::stringstream str_stream(temp_str);
            if (str_stream.str() == "c" || str_stream.str() == "continue")
                is_time_dep = false;
            else
            {
                if (!(str_stream >> dim_idx_map[0]))
                {
                    showErrorMessage(0);
                    dim_idx_map[0] = std::numeric_limits<std::size_t>::max();
                    continue;
                }
                if (dim_idx_map[0] > n_dims-1)
                {
                    showErrorMessage(1, var.getDimCount() - 1);
                    dim_idx_map[0] = std::numeric_limits<std::size_t>::max();
                }
            }
        }
    }
    else
        is_time_dep = false;

    // get spatial dimension(s)
    std::size_t const start_idx = (is_time_dep) ? 1 : 0;
    for (std::size_t i=start_idx; i<n_dims; ++i)
    {
        dim_idx_map[i] = std::numeric_limits<std::size_t>::max();

        std::string request_str ("Enter ID for dimension " + std::to_string(i) + ": ");
        std::size_t idx = parseInput(request_str, var.getDimCount(), true);

        if (idx == var.getDimCount())
        {
            showArraysDims(var);
            i--;
            continue;
        }
        dim_idx_map[i] = idx;
    }

    return is_time_dep;
}

std::pair<std::size_t, std::size_t> timestepSelectionLoop(NcFile const& dataset, NcVar const& var, std::size_t dim_idx)
{
    NcVar const& dim_var = getDimVar(dataset, var, dim_idx);
    std::size_t const n_time_steps = var.getDim(dim_idx).getSize();
    std::pair<std::size_t, std::size_t> bounds(
        std::numeric_limits<std::size_t>::max(),
        std::numeric_limits<std::size_t>::max());
    std::cout << "\nThe dataset contains " << n_time_steps << " time steps.\n";
    bounds.first = parseInput("Specify first time step to export: ", n_time_steps, false);
    while (bounds.first > bounds.second || bounds.second > n_time_steps)
        bounds.second = parseInput("Specify last time step to export: ", n_time_steps, false);
    return bounds;
}

MeshLib::MeshElemType elemSelectionLoop(std::size_t const dim)
{
    if (dim ==1)
        return MeshLib::MeshElemType::LINE;

    MeshLib::MeshElemType t = MeshLib::MeshElemType::INVALID;
    while (t == MeshLib::MeshElemType::INVALID)
    {
        std::cout << "\nSelect element type for result, choose ";

        if (dim==2) std::cout << "(t)riangle or (q)uadliteral: ";
        if (dim==3) std::cout << "(p)rism or (h)exahedron: ";
        std::string type ("");
        std::getline(std::cin, type);
        checkExit(type);
        if (dim==2)
        {
            if (type != "t" && type != "q" &&
                type != "tri" && type != "quad" &&
                type != "triangle" && type != "quatliteral")
                continue;
            if (type == "t" || type == "tri" || type == "triangle")
                return MeshLib::MeshElemType::TRIANGLE;
            else
                return MeshLib::MeshElemType::QUAD;
        }

        if (dim==3)
        {
            if (type != "p" || type != "h" ||
                type != "prism" || type != "hex" ||
                type != "hexahedron")
                continue;
            if (type == "p" || type == "prism")
                return MeshLib::MeshElemType::PRISM;
            else
                return MeshLib::MeshElemType::HEXAHEDRON;
        }
    }
    return t;
}

bool multFilesSelectionLoop(std::pair<std::size_t, std::size_t> const& time_bounds)
{
    std::size_t n_time_steps(time_bounds.second - time_bounds.first + 1);
    std::cout << "\nThe selection includes " << n_time_steps << " time steps.\n";
    std::cout << "0. Save data in " << n_time_steps << " mesh files with one scalar array each.\n";
    std::cout << "1. Save data in one mesh file with " << n_time_steps << " scalar arrays.\n";
    std::size_t ret = parseInput("Select preferred method: ", 2, false);
    return (ret == 0) ? true : false;
}

std::string getIterationString(std::size_t i, std::size_t max)
{
    std::size_t const max_length (std::to_string(max).length());
    std::string const current_str (std::to_string(i));
    return std::string(max_length - current_str.length(), '0') + current_str;
}

double getResolution(NcFile const& dataset, NcVar const& var)
{
    std::size_t dim_idx = var.getDimCount() - 1;
    auto const bounds = getBoundaries(getDimVar(dataset, var, dim_idx));
    return fabs(bounds.second - bounds.first) / static_cast<double>(var.getDim(dim_idx).getSize());
}

GeoLib::RasterHeader createRasterHeader(
    NcFile const& dataset, NcVar const& var,
    std::array<std::size_t, 4> const& dim_idx_map,
    std::vector<std::size_t> const& length, bool is_time_dep)
{
    MathLib::Point3d origin = getOrigin(dataset, var);
    double const res = getResolution(dataset, var);
    std::size_t n_dims = var.getDimCount();
    std::size_t temp_offset = (is_time_dep) ? 1 : 0;
    std::size_t z_length = (n_dims - temp_offset == 3) ? length[dim_idx_map.back()] : 1;
    return {length[dim_idx_map[0 + temp_offset]], length[dim_idx_map[1 + temp_offset]], z_length, origin, res, -9999};
}

int main (int argc, char* argv[])
{
    ApplicationsLib::LogogSetup logog_setup;

    TCLAP::CmdLine cmd(
        "Converts NetCDF into mesh file(s).\n\n "
        "OpenGeoSys-6 software, version " +
            GitInfoLib::GitInfo::ogs_version +
            ".\n"
            "Copyright (c) 2012-2019, OpenGeoSys Community "
            "(http://www.opengeosys.org)",
        ' ', GitInfoLib::GitInfo::ogs_version);

    TCLAP::ValueArg<std::string> input("i", "input-file",
                                       "the NetCDF input file", true,
                                       "", "the NetCDF input file");
    cmd.add(input);
    TCLAP::ValueArg<std::string> output("o", "output-file",
                                        "the OGS mesh output file", true,
                                        "", "the OGS mesh output file");
    cmd.add(output);
    cmd.parse(argc, argv);

    NcFile dataset(input.getValue().c_str(), NcFile::read);

    if (dataset.isNull())
    {
        ERR("Error opening file.");
        return -1;
    }

    std::cout << "OpenGeoSys NetCDF Converter\n";
    std::cout << "File " << input.getValue() << " loaded. Press ENTER to display available data arrays.\n";
    std::cin.ignore();
    std::string output_name (output.getValue());

    std::string const& var_name = arraySelectionLoop(dataset);
    NcVar& var = dataset.getVar(var_name);

    std::array<std::size_t, 4> dim_idx_map;
    bool const is_time_dep = dimensionSelectionLoop(var, dim_idx_map);

    std::size_t const temp_offset = (is_time_dep) ? 1 : 0;
    std::size_t const n_dims = (var.getDimCount());
    std::vector<std::size_t> offset(n_dims, 0);
    std::vector<std::size_t> length(n_dims, 1);
    std::size_t array_length = 1;
    for (std::size_t i=temp_offset; i<n_dims; ++i)
    {
        length[i] = var.getDim(i).getSize();
        array_length *= length[i];
    }

    GeoLib::RasterHeader header =
        createRasterHeader(dataset, var, dim_idx_map, length, is_time_dep);

    std::pair<std::size_t, std::size_t> time_bounds(0,0);
    if (is_time_dep)
        time_bounds = timestepSelectionLoop(dataset, var, dim_idx_map[0]);

    bool mult_files (false);
    if (is_time_dep && time_bounds.first != time_bounds.second)
        mult_files = multFilesSelectionLoop(time_bounds);

    std::unique_ptr<MeshLib::Mesh> mesh;
    MeshLib::MeshElemType const meshElemType = elemSelectionLoop(n_dims - temp_offset);
    for (std::size_t i=time_bounds.first; i<=time_bounds.second; ++i)
    {
        std::cout << "Converting time step " << i << "...\n";
        offset[0] = i;
        std::vector<double> data_array(array_length,0);
        var.getVar(offset, length, data_array.data());
        std::replace_if(data_array.begin(), data_array.end(),
                        [](double& x) { return x <= -9999; }, -9999);

        // reverse lines in vertical direction if the original file has its origin in the northwest corner
        NcVar const& dim_var = getDimVar(dataset, var, n_dims - 1);
        auto const bounds = getBoundaries(dim_var);
        if (bounds.first > bounds.second)
            reverseNorthSouth(data_array.data(), length[n_dims-2], length[n_dims-1]);

        MeshLib::UseIntensityAs const useIntensity = MeshLib::UseIntensityAs::DATAVECTOR;
        if (mult_files)
        {
            mesh.reset(MeshLib::RasterToMesh::convert(data_array.data(), header, meshElemType, useIntensity, var.getName()));
            std::string const output_file_name (BaseLib::dropFileExtension(output_name) + getIterationString(i, time_bounds.second) + ".vtu");
            MeshLib::IO::VtuInterface vtu(mesh.get());
            vtu.writeToFile(output_file_name);
        }
        else
        {
            std::string array_name (var.getName());
            if (time_bounds.first != time_bounds.second)
                array_name.append(getIterationString(i, time_bounds.second));
            if (i==time_bounds.first) // create persistent mesh
                mesh.reset(MeshLib::RasterToMesh::convert(data_array.data(), header, meshElemType, useIntensity, array_name));
            else // copy array to mesh
            {
                std::unique_ptr<MeshLib::Mesh> temp (MeshLib::RasterToMesh::convert(data_array.data(), header, meshElemType, useIntensity, array_name));
                MeshLib::PropertyVector<double> const*const vec = temp->getProperties().getPropertyVector<double>(array_name);
                if (vec == nullptr)
                    return EXIT_FAILURE;
                MeshLib::addPropertyToMesh<double>(*mesh, array_name, MeshLib::MeshItemType::Cell, 1, *vec);
            }
            if (i==time_bounds.second)
            {
                MeshLib::IO::VtuInterface vtu(mesh.get());
                vtu.writeToFile(output_name);
            }
        }
    }

    std::cout << "Conversion finished successfully.\n";
    return EXIT_SUCCESS;
}
