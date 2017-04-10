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

#include <netcdfcpp.h>

// BaseLib
#include "BaseLib/BuildInfo.h"
#include "BaseLib/LogogSimpleFormatter.h"
#include "BaseLib/FileTools.h"

// GeoLib
#include "GeoLib/Raster.h"

// MeshLib
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshGenerators/RasterToMesh.h"

#include "MeshLib/IO/VtkIO/VtuInterface.h"


void checkExit(std::string const& str)
{
    if (str == "exit")
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

void showArrays(NcFile const& dataset)
{
    std::size_t const n_vars (dataset.num_vars());
    std::cout << "The NetCDF file contains the following " << n_vars << " arrays:\n\n";
    std::cout << "\tIndex\tArray Name\t#Dimensions\n";
    std::cout << "-------------------------------------------\n";
    for (std::size_t i=0; i<n_vars; ++i)
    {
        NcVar const& var = *dataset.get_var(i);
        std::cout << "\t" << i << "\t" << var.name() << "\t(" << var.num_dims() << "D array)\n";
    }
    std::cout << "\n";
}

bool showArraysDims(NcVar const& var)
{
    std::cout << "Data array \"" << var.name() << "\" contains the following dimensions:\n";
    std::size_t const n_dims (var.num_dims());
    for (std::size_t i=0; i<n_dims; ++i)
        std::cout << "\t" << i << "\t" << var.get_dim(i)->name() << "\t(" << var.get_dim(i)->size() << " values)\n";
    std::cout << "\n";
    return true;
}

std::size_t getDimensionBoundaries(NcFile const& dataset, NcVar const& var, std::size_t dim, double &start, double &end)
{
    NcVar& dim_var = *dataset.get_var(var.get_dim(dim)->name());
    if ((dim_var.num_dims()) == 1)
    {
        std::size_t size = dim_var.get_dim(0)->size();

        std::size_t length[1] = {1};
        long idx[1] = {0};
        dim_var.set_cur(idx);
        double val_at_idx[1] = {0};
        dim_var.get(val_at_idx, length);
        start = val_at_idx[0];

        idx[0] = size-1;
        dim_var.set_cur(idx);
        dim_var.get(val_at_idx, length);
        end = val_at_idx[0];
        return size;
    }
    return 0;
}

MathLib::Point3d getOrigin(std::array<std::pair<double, double>,3> bounds, std::size_t n_dims)
{
    MathLib::Point3d origin(std::array<double, 3>{ {0, 0, 0}});
    for (std::size_t i=0; i<n_dims; ++i)
        origin[i] = (bounds[i].first < bounds[i].second) ? bounds[i].first : bounds[i].second;
    if (n_dims > 1)
        std::swap(origin[0], origin[1]);
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

std::size_t arraySelectionLoop(NcFile const& dataset)
{
    showArrays(dataset);
    std::size_t idx = parseInput("Enter data array index: ", dataset.num_vars(), true);

    if (idx == dataset.num_vars())
        return arraySelectionLoop(dataset);

    return idx;
}

bool dimensionSelectionLoop(NcVar const& var, std::array<std::size_t,4> &dim_idx_map)
{
    showArraysDims(var);
    std::size_t const n_dims (var.num_dims());
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
                    continue;
                }
                if (dim_idx_map[0] > n_dims-1)
                {
                    showErrorMessage(1, var.num_dims()-1);
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
        std::size_t idx = parseInput(request_str, var.num_dims(), true);

        if (idx == var.num_dims())
        {
            showArraysDims(var);
            i--;
            continue;
        }
        dim_idx_map[i] = idx;
    }

    return is_time_dep;
}

std::pair<std::size_t, std::size_t> timestepSelectionLoop(NcFile const& dataset, NcVar const& var, std::size_t dim)
{
    double start, end;
    std::string str("");
    std::pair<std::size_t, std::size_t> bounds(0,std::numeric_limits<std::size_t>::max());
    std::size_t n_time_steps = getDimensionBoundaries(dataset, var, dim, start, end);
    std::cout << "\nThe dataset contains " << n_time_steps << " time steps.\n";
    bounds.first = parseInput("Specify first time step to export: ", n_time_steps, false);
    while (bounds.first > bounds.second || bounds.second == std::numeric_limits<std::size_t>::max())
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
    std::size_t n_time_steps(time_bounds.second - time_bounds.first);
    std::cout << "\nThe selection includes " << n_time_steps << " time steps.\n";
    std::cout << "0. Save data in " << n_time_steps << " mesh files with one scalar array each.\n";
    std::cout << "1. Save data in one mesh file with " << n_time_steps << " scalar arrays.\n";
    std::size_t ret = parseInput("Select preferred method: ", 2, false);
    return (ret == 0) ? true : false;
}

void setTimeStep(NcVar &var, std::size_t time_step)
{
    long* newOrigin = new long[var.num_dims()];
    for (int i=0; i < var.num_dims(); ++i) 
        newOrigin[i]=0;
    newOrigin[0] = time_step;
    var.set_cur(newOrigin);
    delete [] newOrigin;
}

std::string getIterationString(std::size_t i, std::size_t max)
{
    std::size_t const max_length (std::to_string(max).length());
    std::string const current_str (std::to_string(i));
    return std::string(max_length - current_str.length(), '0') + current_str;
}

int main (int argc, char* argv[])
{
    ApplicationsLib::LogogSetup logog_setup;

    TCLAP::CmdLine cmd("Converts NetCDF into mesh file.", ' ', BaseLib::BuildInfo::git_describe);
    TCLAP::ValueArg<std::string> input("i", "input-file",
                                       "the NetCDF input file", true,
                                       "", "the NetCDF input file");
    cmd.add(input);
    TCLAP::ValueArg<std::string> output("o", "output-file",
                                        "the OGS mesh output file", true,
                                        "", "the OGS mesh output file");
    cmd.add(output);
    cmd.parse(argc, argv);

    NcFile dataset(input.getValue().c_str(), NcFile::ReadOnly);

    if (!dataset.is_valid())
    {
        ERR("Error opening file.");
        return -1;
    }

    std::cout << "OpenGeoSys NetCDF Converter\n";
    std::cout << "File " << input.getValue() << " loaded. Press ENTER to display available data arrays.\n";
    std::cin.ignore();
    std::string output_name (output.getValue());

    std::size_t const var_idx = arraySelectionLoop(dataset);
    NcVar& var = *dataset.get_var(var_idx);

    std::array<std::size_t, 4> dim_idx_map;
    bool const is_time_dep = dimensionSelectionLoop(var, dim_idx_map);

    std::size_t temp_offset = (is_time_dep) ? 1 : 0;
    std::size_t const n_dims = (var.num_dims());
    std::size_t* length = new std::size_t[n_dims];
    for(int i=0; i < n_dims; i++) 
        length[i]=1;
    std::array<std::size_t, 3> spatial_size = {0,0,0};
    std::array<std::pair<double, double>,3> spatial_bounds;
    for (std::size_t i=0; i<3; ++i)
        spatial_bounds[i] = std::pair<double, double>(0, 0);
    std::size_t array_length = 1;
    for (std::size_t i=temp_offset; i<n_dims; ++i)
    {
        length[dim_idx_map[i]] = getDimensionBoundaries(dataset, var, dim_idx_map[i], spatial_bounds[i-temp_offset].first, spatial_bounds[i-temp_offset].second);
        spatial_size[i-temp_offset] = length[dim_idx_map[i]];
        array_length *= length[dim_idx_map[i]];
    }

    MathLib::Point3d origin = getOrigin(spatial_bounds, n_dims);
    double const resolution = fabs(spatial_bounds[0].second - spatial_bounds[0].first) / (spatial_size[0] - 1);
    GeoLib::RasterHeader header = { spatial_size[1], spatial_size[0], origin, resolution, -9999 };
    MeshLib::MeshElemType meshElemType = elemSelectionLoop(n_dims - temp_offset);
    MeshLib::UseIntensityAs useIntensity = MeshLib::UseIntensityAs::DATAVECTOR;

    std::pair<std::size_t, std::size_t> time_bounds(0,0);
    if (is_time_dep)
        time_bounds = timestepSelectionLoop(dataset, var, dim_idx_map[0]);

    bool mult_files (false);
    if (is_time_dep && time_bounds.first != time_bounds.second)
        mult_files = multFilesSelectionLoop(time_bounds);

    double* data_array = nullptr;
    std::unique_ptr<MeshLib::Mesh> mesh;
    for (std::size_t i=time_bounds.first; i<=time_bounds.second; ++i)
    {
        if (is_time_dep)
            setTimeStep(var, i);

        data_array = new double[array_length];
        for (std::size_t i = 0; i < array_length; i++)
            data_array[i] = 0;
        var.get(data_array, length);

        for (std::size_t i=0; i < array_length; i++)
            if (data_array[i] < -9999 ) data_array[i] = -9999; // all values < -10000, set to "no-value"

        // reverse lines in vertical direction if the original file has its origin in the northwest corner
        if (spatial_bounds[0].first > spatial_bounds[0].second)
            reverseNorthSouth(data_array, spatial_size[1], spatial_size[0]);

        if (mult_files)
        {
            mesh.reset(MeshLib::RasterToMesh::convert(data_array, header, meshElemType, useIntensity, var.name()));
            std::string const output_file_name (BaseLib::dropFileExtension(output_name) + getIterationString(i, time_bounds.second) + ".vtu");
            MeshLib::IO::VtuInterface vtu(mesh.get());
            vtu.writeToFile(output_file_name);
        }
        else
        {
            std::string array_name (var.name());
            if (time_bounds.first != time_bounds.second)
                array_name.append(getIterationString(i, time_bounds.second));
            if (i==time_bounds.first) // create persistent mesh
                mesh.reset(MeshLib::RasterToMesh::convert(data_array, header, meshElemType, useIntensity, array_name));
            else
            {
                //copy array to mesh
            }
            if (i==time_bounds.second)
            {
                MeshLib::IO::VtuInterface vtu(mesh.get());
                vtu.writeToFile(output_name);
            }
        }
    }


    /***********************
    for (std::size_t j = time_bounds.first; j <= time_bounds.second; ++j)
    {
            boost::optional< MeshLib::PropertyVector<double>& > tmp_vec = tmp_mesh->getProperties().getPropertyVector<double>(var.name());
            if (!tmp_vec)
            {
                ERR("Vector not found");
                continue;
            }

            if (tmp_vec->size() != mesh->getNumberOfElements())
            {
                ERR("Vector size doesn't fit.");
                continue;
            }

            MeshLib::Properties& props (mesh->getProperties());

            std::string iter_str (getIterationString(j, time_bounds.second);

            boost::optional< MeshLib::PropertyVector<double>& > vec = 
                props.createNewPropertyVector<double>(std::string(var.name() + iter_str), MeshLib::MeshItemType::Cell, 1);

            if (!vec)
            {
                ERR ("Error creating vector");
                continue;
            }
            vec->reserve(tmp_vec->size());
            for (std::size_t i=0; i<tmp_vec->size(); ++i)
                vec->push_back((*tmp_vec)[i]);

            std::cout << "time step " << j << " written.\n";
    }
    /***********************/

    delete[] length;
    delete[] data_array;
    return 0;
}
