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

// GeoLib
#include "GeoLib/Raster.h"

// MeshLib
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshGenerators/RasterToMesh.h"

#include "MeshLib/IO/VtkIO/VtuInterface.h"


void showErrorMessage(std::size_t error_id, std::size_t max = 0)
{
    if (error_id == 0)
    {
        ERR ("Index not recognised.");
        std::cout << "\"info\" will display the available options again. \"exit\" will exit the programme.\n";
    }
    else if (error_id == 1)
    {
        ERR ("Index not valid. Valid indices are in [0,%d].", max);
    }
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

bool showArraysDims(NcFile const& dataset, NcVar const& var)
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
    NcVar dim_var = *dataset.get_var(var.get_dim(dim)->name());
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

void reverseNorthSouth(double* data, std::size_t width, std::size_t height)
{
    double* cp_array = new double[width*height];

    for (std::size_t i=0; i<height; i++)
    {
        for (std::size_t j=0; j<width; j++)
        {
            std::size_t const old_index ((width*height)-(width*(i+1)));
            std::size_t const new_index (width*i);
            cp_array[new_index+j] = data[old_index+j];
        }
    }

    std::size_t const length(height*width);
    for (std::size_t i=0; i<length; i++)
        data[i] = cp_array[i];

    delete[] cp_array;
}

std::size_t arraySelectionLoop(NcFile const& dataset)
{
    showArrays(dataset);
    while (true)
    {
        cout << "Enter data array index: ";
        std::string var_idx_str ("");
        std::getline(std::cin, var_idx_str);
        if (var_idx_str == "info")
        {
            showArrays(dataset);
            continue;
        }
        if (var_idx_str == "exit")
            exit(0);
    
        std::stringstream str_stream(var_idx_str);
        std::size_t var_idx;
        if (!(str_stream >> var_idx))
        {
            showErrorMessage(0);
            continue;
        }
        if (var_idx > dataset.num_vars()-1)
        {
            showErrorMessage(1, dataset.num_vars()-1);
            continue;
        }
        return var_idx;
    }
    return std::numeric_limits<std::size_t>::max();
}

bool dimensionSelectionLoop(NcFile const& dataset, NcVar const& var, std::array<std::size_t,4> &dim_idx_map)
{
    showArraysDims(dataset, var);
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
            if (str_stream.str() == "c")
                is_time_dep = false;
            else 
            {
                if (str_stream >> dim_idx_map[0])
                {
                    if (dim_idx_map[0] > n_dims-1)
                    {
                        showErrorMessage(1, var.num_dims()-1);
                        dim_idx_map[0] = std::numeric_limits<std::size_t>::max();
                    }
                }
            }
        }
    }

    // get spatial dimension(s)
    std::size_t const start_idx = (is_time_dep) ? 1 : 0;
    for (std::size_t i=start_idx; i<n_dims; ++i)
    {
        dim_idx_map[i] = std::numeric_limits<std::size_t>::max();
        while (dim_idx_map[i] == std::numeric_limits<std::size_t>::max())
        {
            std::string dim_idx_str;
            cout << "Enter ID for dimension " << i << ": ";
            std::getline(std::cin, dim_idx_str);
            if (dim_idx_str == "info")
            {
                showArraysDims(dataset, var);
                continue;
            }
            if (dim_idx_str == "exit")
                exit(0);

            std::stringstream str_stream(dim_idx_str);
            std::size_t dim_idx;
            if (!(str_stream >> dim_idx))
            {
                showErrorMessage(0);
                continue;
            }
            if (dim_idx > var.num_dims()-1)
            {
                showErrorMessage(1, var.num_dims()-1);
                continue;
            }
            dim_idx_map[i] = dim_idx;
        }
    }

    return is_time_dep;
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

        if (dim==2)
        {
            if (type != "t" || type != "q" ||
                type != "tri" || type != "quad" ||
                type != "triangle" || type != "quatliteral")
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

void setOrigin(NcVar &var, std::size_t time_step)
{
    long* newOrigin = new long[var.num_dims()];
    for (int i=0; i < var.num_dims(); ++i) 
        newOrigin[i]=0;
    newOrigin[0] = time_step;
    var.set_cur(newOrigin);
    delete [] newOrigin;
}

/*
MeshLib::Mesh convert(double* data_array, std::array<std::size_t,3> length, GeoLib::Point &origin, double resolution, MeshLib::MeshElemType meshElemType, MeshLib::UseIntensityAs useIntensity, std::string &name)
{
    return 
}
*/
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

    std::cout << "OpenGeoSys NetCDF Converter\n";
    std::cout << "File " << input.getValue() << " loaded. Press ENTER to display available data arrays.\n";
    std::cin.ignore();

    std::size_t const var_idx = arraySelectionLoop(dataset);

    std::array<std::size_t, 4> dim_idx_map;
    NcVar var = *dataset.get_var(var_idx);
    bool const is_time_dep = dimensionSelectionLoop(dataset, var, dim_idx_map);

    std::cin.ignore();

    std::size_t temp_offset = (is_time_dep) ? 1 : 0;
    std::size_t const n_dims = (var.num_dims());
    std::size_t* length = new std::size_t[n_dims];
    for(int i=0; i < n_dims; i++) 
        length[i]=1;
    std::array<double, 3> spatial_size = {0,0,0};
    std::array<double, 3> start = {0,0,0};
    std::array<double, 3> end   = {0,0,0};
    std::size_t array_length = 1;
    for (std::size_t i=temp_offset; i<n_dims; ++i)
    {
        length[i] = getDimensionBoundaries(dataset, var, dim_idx_map[i], start[i-temp_offset], end[i-temp_offset]);
        spatial_size[i-temp_offset] = length[i];
        array_length *= length[i];
    }

    setOrigin(var, 0);

    double* data_array = new double[array_length];
    for (std::size_t i=0; i < array_length; i++) 
        data_array[i]=0;
    var.get(data_array, length); //create Array of Values

    for (std::size_t i=0; i < array_length; i++)
        if (data_array[i] < -9999 ) data_array[i] = -9999; // all values < -10000, set to "no-value"

    double origin_x = (start[1] < end[1]) ? start[1] : end[1];
    double origin_y = (start[0] < end[0]) ? start[0] : end[0];
    MathLib::Point3d origin (std::array<double,3>{{origin_x, origin_y, 0}});
    double const resolution = fabs(end[0]-start[0])/(spatial_size[0]-1);

    // reverse lines in vertical direction if the original file has its origin in the northwest corner
    if (start[0] > end[0])
        reverseNorthSouth(data_array, spatial_size[1], spatial_size[0]);

    GeoLib::RasterHeader header = {spatial_size[1], spatial_size[0], origin, resolution, -9999};

    //MeshLib::MeshElemType meshElemType = elemSelectionLoop(n_dims-temp_offset);

    MeshLib::MeshElemType meshElemType = MeshLib::MeshElemType::TRIANGLE;
    MeshLib::UseIntensityAs useIntensity = MeshLib::UseIntensityAs::DATAVECTOR;

    MeshLib::Mesh* mesh = MeshLib::RasterToMesh::convert(data_array, header, meshElemType, useIntensity, var.name());

    /***********************
    std::size_t const n_time_steps = getDimensionBoundaries(dataset, var, dim_idx_map[2], start_lat, end_lat);
    std::size_t const array_length (sizeLat * sizeLon);
    for (std::size_t j=1; j<n_time_steps; ++j)
    {
        setOrigin(var, j);
        var.get(data_array, length); //create Array of Values

        for (std::size_t i=0; i < (sizeLat*sizeLon); i++)
            if (data_array[i] < -9999 ) data_array[i] = -9999; // all values < -10000, set to "no-value"

        double const resolution = fabs(end_lat-start_lat)/(sizeLat-1);

        // reverse lines in vertical direction if the original file has its origin in the northwest corner
        if (start_lat > end_lat)
            reverseNorthSouth(data_array, sizeLon, sizeLat);
        GeoLib::RasterHeader header = {sizeLon, sizeLat, origin, resolution, -9999};

        MeshLib::MeshElemType meshElemType = MeshLib::MeshElemType::TRIANGLE;
        MeshLib::UseIntensityAs useIntensity = MeshLib::UseIntensityAs::DATAVECTOR;
        std::unique_ptr<MeshLib::Mesh> tmp_mesh (MeshLib::RasterToMesh::convert(data_array, header, meshElemType, useIntensity, var.name()));

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

        std::string prefix ("");
        if (j<100) prefix = "0";
        if (j<10) prefix = "00";
        std::string var_name_str (prefix + std::to_string(j));

        boost::optional< MeshLib::PropertyVector<double>& > vec = 
            props.createNewPropertyVector<double>(std::string(var.name() + var_name_str), MeshLib::MeshItemType::Cell, 1);

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

    MeshLib::IO::VtuInterface vtu(mesh);
    vtu.writeToFile("d:/nctest.vtu");

    delete[] length;
    delete[] data_array;


    return 0;
}
