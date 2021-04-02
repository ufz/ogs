/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

// STL
#include <tclap/CmdLine.h>

#include <cctype>
#include <iostream>
#include <limits>
#include <memory>
#include <netcdf>
#include <numeric>
#include <sstream>
#include <string>
#include <utility>

#include "Applications/FileIO/AsciiRasterInterface.h"
#include "BaseLib/FileTools.h"
#include "BaseLib/Logging.h"
#include "GeoLib/Raster.h"
#include "InfoLib/GitInfo.h"
#include "MeshLib/IO/VtkIO/VtuInterface.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshGenerators/RasterToMesh.h"

using namespace netCDF;

static const double no_data_output = -9999;
static double no_data_input = -9999;

enum class OutputType
{
    INVALID,
    IMAGES,
    SINGLEMESH,
    MULTIMESH
};

static void checkExit(std::string const& str)
{
    if (str == "x" || str == "exit")
        exit(0);
}

static void showErrorMessage(std::size_t const error_id,
                             std::size_t const max = 0)
{
    if (error_id == 0)
    {
        ERR("Input not valid.");
    }
    else if (error_id == 1)
    {
        ERR("Index not valid. Valid indices are in [0,{:d}].", max);
    }
    else if (error_id == 2)
    {
        ERR("Input not valid.");
        std::cout << "Type \"info\" to display the available options again. "
                     "\"exit\" will exit the programme.\n";
    }
}

static std::size_t parseInput(std::string const& request_str,
                              std::size_t const max_val,
                              bool const has_info = false)
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
            std::size_t const error_val = (has_info) ? 2 : 0;
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

static NcVar getDimVar(NcFile const& dataset, NcVar const& var,
                       std::size_t const dim)
{
    NcDim const& dim_obj = var.getDim(dim);
    return dataset.getVar(dim_obj.getName());
}

static std::pair<double, double> getDimLength(NcVar const& var,
                                              std::size_t const dim)
{
    return std::make_pair(0.0, static_cast<double>(var.getDim(dim).getSize()));
}

static std::vector<std::string> getArrays(NcFile const& dataset)
{
    auto const& names = dataset.getVars();
    std::vector<std::string> var_names;
    for (auto [name, var] : names)
    {
        (void)var;
        var_names.push_back(name);
    }
    return var_names;
}

static void showArrays(NcFile const& dataset)
{
    std::size_t const n_vars(dataset.getDimCount());
    std::cout << "The NetCDF file contains the following " << n_vars
              << " arrays:\n\n";
    std::cout << "\tIndex\tArray Name\t#Dimensions\n";
    std::cout << "-------------------------------------------\n";
    auto const& names = dataset.getVars();
    std::size_t count = 0;
    for (auto [name, var] : names)
    {
        std::cout << "\t" << count++ << "\t" << name << "\t("
                  << var.getDimCount() << "D array)\n";
    }
    std::cout << "\n";
}

static void showArraysDims(NcVar const& var)
{
    std::cout << "Data array \"" << var.getName()
              << "\" contains the following dimensions:\n";
    std::size_t const n_dims(var.getDimCount());
    for (std::size_t i = 0; i < n_dims; ++i)
        std::cout << "\t" << i << "\t" << var.getDim(i).getName() << "\t("
                  << var.getDim(i).getSize() << " values)\n";
    std::cout << "\n";
}

static std::pair<double, double> getBoundaries(NcVar const& var)
{
    if (var.getDimCount() == 1)
    {
        double start, end;
        std::size_t const size = var.getDim(0).getSize();
        var.getVar({0}, {1}, &start);
        var.getVar({size - 1}, {1}, &end);
        return std::make_pair(start, end);
    }
    return std::make_pair(0, 0);
}

static MathLib::Point3d getOrigin(NcFile const& dataset, NcVar const& var,
                                  std::vector<std::size_t> const& dim_idx_map,
                                  std::size_t const time_offset)
{
    MathLib::Point3d origin(MathLib::ORIGIN);
    std::size_t const n_dims = var.getDimCount();
    for (std::size_t i = time_offset; i < n_dims; ++i)
    {
        NcVar const& dim = getDimVar(dataset, var, dim_idx_map[i]);
        auto const bounds = (dim.isNull()) ? getDimLength(var, dim_idx_map[i])
                                           : getBoundaries(dim);
        origin[i - time_offset] =
            (bounds.first < bounds.second) ? bounds.first : bounds.second;
    }
    return origin;
}

static void flipRaster(std::vector<double>& data, std::size_t const layers,
                       std::size_t const width, std::size_t const height)
{
    std::size_t const length(data.size());
    std::vector<double> tmp_vec;
    tmp_vec.reserve(length);
    for (std::size_t k = 0; k < layers; k++)
    {
        std::size_t const layer_end = (k + 1) * height * width;
        for (std::size_t i = 0; i < height; i++)
        {
            std::size_t const line_idx(layer_end - (width * (i + 1)));
            for (std::size_t j = 0; j < width; j++)
            {
                tmp_vec.push_back(data[line_idx + j]);
            }
        }
    }
    std::copy(tmp_vec.cbegin(), tmp_vec.cend(), data.begin());
}

static bool canConvert(NcVar const& var)
{
    bool ret(var.getDimCount() < 2);
    if (ret)
        ERR("Only 2+ dimensional variables can be converted into OGS "
            "Meshes.\n");
    return !ret;
}

static std::string arraySelectionLoop(NcFile const& dataset)
{
    std::vector<std::string> const& names = getArrays(dataset);
    showArrays(dataset);
    std::size_t const idx =
        parseInput("Enter data array index: ", dataset.getVarCount(), true);

    if (static_cast<int>(idx) == dataset.getVarCount() ||
        !canConvert(dataset.getVar(names[idx])))
        return arraySelectionLoop(dataset);

    return names[idx];
}

static bool dimensionSelectionLoop(NcVar const& var,
                                   std::vector<std::size_t>& dim_idx_map)
{
    showArraysDims(var);
    std::size_t const n_dims(var.getDimCount());
    dim_idx_map[0] = std::numeric_limits<std::size_t>::max();
    bool is_time_dep(true);

    // get temporal dimension
    if (n_dims > 1)
    {
        std::string temp_str("");
        cout << "Is the parameter time-dependent?\n";
        while (dim_idx_map[0] == std::numeric_limits<std::size_t>::max() &&
               is_time_dep == true)
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
                if (dim_idx_map[0] > n_dims - 1)
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
    std::array<std::string, 4> const dim_comment{
        "(x / longitude)", "(y / latitude)", "(z / height / depth)",
        "[Error: 4-dimensional non-temporal arrays are not supported]"};
    for (std::size_t i = start_idx; i < n_dims; ++i)
    {
        dim_idx_map[i] = std::numeric_limits<std::size_t>::max();

        std::string const request_str("Enter ID for dimension " +
                                      std::to_string(i) + " " +
                                      dim_comment[i - start_idx] + ": ");
        std::size_t const idx =
            parseInput(request_str, var.getDimCount(), true);

        if (static_cast<int>(idx) == var.getDimCount())
        {
            showArraysDims(var);
            i--;
            continue;
        }
        dim_idx_map[i] = idx;
    }

    return is_time_dep;
}

static std::pair<std::size_t, std::size_t> timestepSelectionLoop(
    NcVar const& var, std::size_t const dim_idx)
{
    std::size_t const n_time_steps = var.getDim(dim_idx).getSize();
    std::size_t const max_val = std::numeric_limits<std::size_t>::max();
    std::pair<std::size_t, std::size_t> bounds(max_val, max_val);
    std::cout << "\nThe dataset contains " << n_time_steps << " time steps.\n";
    while (bounds.first == max_val)
    {
        bounds.first = parseInput(
            "Specify first time step to export: ", n_time_steps, false);
    }
    while (bounds.first > bounds.second || bounds.second > n_time_steps)
    {
        bounds.second = parseInput(
            "Specify last time step to export: ", n_time_steps, false);
    }
    return bounds;
}

static MeshLib::MeshElemType elemSelectionLoop(std::size_t const dim)
{
    if (dim == 1)
        return MeshLib::MeshElemType::LINE;

    MeshLib::MeshElemType t = MeshLib::MeshElemType::INVALID;
    while (t == MeshLib::MeshElemType::INVALID)
    {
        std::cout << "\nSelect element type for result, choose ";

        if (dim == 2)
            std::cout << "(t)riangle or (q)uadliteral: ";
        if (dim == 3)
            std::cout << "(p)rism or (h)exahedron: ";
        std::string type("");
        std::getline(std::cin, type);
        checkExit(type);
        if (dim == 2)
        {
            if (type != "t" && type != "q" && type != "tri" && type != "quad" &&
                type != "triangle" && type != "quatliteral")
                continue;
            if (type == "t" || type == "tri" || type == "triangle")
                return MeshLib::MeshElemType::TRIANGLE;
            return MeshLib::MeshElemType::QUAD;
        }

        if (dim == 3)
        {
            if (type != "p" && type != "h" && type != "prism" &&
                type != "hex" && type != "hexahedron")
                continue;
            if (type == "p" || type == "prism")
                return MeshLib::MeshElemType::PRISM;
            return MeshLib::MeshElemType::HEXAHEDRON;
        }
    }
    return t;
}

static OutputType multFilesSelectionLoop(
    std::pair<std::size_t, std::size_t> const& time_bounds)
{
    OutputType t = OutputType::INVALID;
    while (t == OutputType::INVALID)
    {
        std::size_t const n_time_steps(time_bounds.second - time_bounds.first +
                                       1);
        std::cout << "\nThe selection includes " << n_time_steps
                  << " time steps.\n";
        std::cout << "0. Save data in " << n_time_steps
                  << " mesh files with one scalar array each.\n";
        std::cout << "1. Save data in one mesh file with " << n_time_steps
                  << " scalar arrays.\n";
        std::cout << "2. Save data as " << n_time_steps << " ASC images.\n";

        std::size_t const ret =
            parseInput("Select preferred method: ", 3, false);

        if (ret == 0)
            t = OutputType::MULTIMESH;
        else if (ret == 1)
            t = OutputType::SINGLEMESH;
        else if (ret == 2)
            t = OutputType::IMAGES;
    }
    return t;
}

static std::string getIterationString(std::size_t i, std::size_t max)
{
    std::size_t const max_length(std::to_string(max).length());
    std::string const current_str(std::to_string(i));
    return std::string(max_length - current_str.length(), '0') + current_str;
}

static double getResolution(NcFile const& dataset, NcVar const& var)
{
    std::size_t const dim_idx = var.getDimCount() - 1;
    NcVar const dim_var(getDimVar(dataset, var, dim_idx));
    auto const bounds = (dim_var.isNull()) ? getDimLength(var, dim_idx)
                                           : getBoundaries(dim_var);
    std::size_t const dim_size = var.getDim(dim_idx).getSize();
    if (dim_size == 0)
    {
        OGS_FATAL("Dimension '{:s}' has size 0. Aborting...",
                  var.getDim(dim_idx).getName());
    }
    return std::fabs(bounds.second - bounds.first) /
           static_cast<double>(dim_size);
}

static GeoLib::RasterHeader createRasterHeader(
    NcFile const& dataset, NcVar const& var,
    std::vector<std::size_t> const& dim_idx_map,
    std::vector<std::size_t> const& length, std::size_t const time_offset)
{
    MathLib::Point3d const origin =
        getOrigin(dataset, var, dim_idx_map, time_offset);
    double const res = getResolution(dataset, var);
    std::size_t n_dims = var.getDimCount();
    std::size_t z_length =
        (n_dims - time_offset == 3) ? length[dim_idx_map.back()] : 1;
    return {length[dim_idx_map[0 + time_offset]],
            length[dim_idx_map[1 + time_offset]],
            z_length,
            origin,
            res,
            no_data_output};
}

static std::vector<std::size_t> getLength(NcVar const& var,
                                          std::size_t const time_offset)
{
    std::size_t const n_dims = (var.getDimCount());
    std::vector<std::size_t> length(n_dims, 1);
    for (std::size_t i = time_offset; i < n_dims; ++i)
    {
        length[i] = var.getDim(i).getSize();
    }
    return length;
}

static std::vector<double> getData(NcFile const& dataset, NcVar const& var,
                                   std::size_t const total_length,
                                   std::size_t const time_step,
                                   std::vector<std::size_t> const& length)
{
    std::size_t const n_dims(var.getDimCount());
    std::vector<std::size_t> offset(n_dims, 0);
    offset[0] = time_step;
    std::vector<double> data_vec(total_length, 0);
    var.getVar(offset, length, data_vec.data());

    std::replace_if(
        data_vec.begin(), data_vec.end(),
        [](double const& x) { return x == no_data_input; }, no_data_output);

    return data_vec;
}

static bool assignDimParams(NcVar const& var,
                            std::vector<std::size_t>& dim_idx_map,
                            TCLAP::ValueArg<std::size_t>& arg_dim_time,
                            TCLAP::ValueArg<std::size_t>& arg_dim1,
                            TCLAP::ValueArg<std::size_t>& arg_dim2,
                            TCLAP::ValueArg<std::size_t>& arg_dim3)
{
    std::size_t dim_param_count = 0;
    if (arg_dim_time.isSet())
        dim_param_count++;
    if (arg_dim1.isSet())
        dim_param_count++;
    if (arg_dim2.isSet())
        dim_param_count++;
    if (arg_dim3.isSet())
        dim_param_count++;

    std::size_t const n_dims = var.getDimCount();
    if (dim_param_count != n_dims)
    {
        ERR("Number of parameters set does not fit number of parameters for "
            "specified variable.");
        return false;
    }

    if (arg_dim_time.getValue() >= n_dims || arg_dim1.getValue() >= n_dims ||
        arg_dim2.getValue() >= n_dims || arg_dim3.getValue() >= n_dims)
    {
        ERR("Maximum allowed dimension for variable \"{:s}\" is {:d}.",
            var.getName(), n_dims - 1);
        return false;
    }

    if (arg_dim_time.isSet())
        dim_idx_map[0] = arg_dim_time.getValue();
    std::size_t const temp_offset = (arg_dim_time.isSet()) ? 1 : 0;
    dim_idx_map[0 + temp_offset] = arg_dim1.getValue();
    dim_idx_map[1 + temp_offset] = arg_dim2.getValue();
    if (n_dims == (3 + temp_offset))
        dim_idx_map[2 + temp_offset] = arg_dim3.getValue();

    return true;
}

static std::pair<std::size_t, std::size_t> assignTimeBounds(
    NcVar const& var,
    TCLAP::ValueArg<std::size_t>& arg_time_start,
    TCLAP::ValueArg<std::size_t>& arg_time_end)
{
    auto const bounds = getBoundaries(var);
    if (arg_time_start.getValue() > bounds.second)
    {
        ERR("Start time step larger than total number of time steps. Resetting "
            "to 0.");
        arg_time_start.reset();
    }

    if (!arg_time_end.isSet())
        return {arg_time_start.getValue(), arg_time_start.getValue()};

    if (arg_time_end.getValue() > bounds.second)
    {
        ERR("End time step larger than total number of time steps. Resetting "
            "to starting time step");
        return {arg_time_start.getValue(), arg_time_start.getValue()};
    }

    if (arg_time_end.getValue() < arg_time_start.getValue())
    {
        ERR("End time step larger than starting time step. Swapping values");
        return {arg_time_end.getValue(), arg_time_start.getValue()};
    }

    return {arg_time_start.getValue(), arg_time_end.getValue()};
}

static MeshLib::MeshElemType assignElemType(
    TCLAP::ValueArg<std::string>& arg_elem_type)
{
    if (arg_elem_type.getValue() == "tri")
        return MeshLib::MeshElemType::TRIANGLE;
    if (arg_elem_type.getValue() == "quad")
        return MeshLib::MeshElemType::QUAD;
    if (arg_elem_type.getValue() == "prism")
        return MeshLib::MeshElemType::PRISM;
    if (arg_elem_type.getValue() == "hex")
        return MeshLib::MeshElemType::HEXAHEDRON;
    // this should never happen
    return MeshLib::MeshElemType::INVALID;
}

static bool convert(NcFile const& dataset, NcVar const& var,
                    std::string const& output_name,
                    std::vector<std::size_t> const& dim_idx_map,
                    std::size_t const time_offset,
                    std::pair<std::size_t, std::size_t> const& time_bounds,
                    OutputType const output,
                    MeshLib::MeshElemType const elem_type)
{
    std::unique_ptr<MeshLib::Mesh> mesh;
    std::vector<std::size_t> const length = getLength(var, time_offset);
    std::size_t const array_length = std::accumulate(
        length.cbegin(), length.cend(), 1, std::multiplies<std::size_t>());
    for (std::size_t i = time_bounds.first; i <= time_bounds.second; ++i)
    {
        std::string const step_str =
            (time_bounds.first != time_bounds.second)
                ? std::string(" time step " + std::to_string(i))
                : "";
        std::cout << "Converting" << step_str << "...\n";
        std::vector<double> data_vec =
            getData(dataset, var, array_length, i, length);

        // reverse lines in vertical direction if file has its origin in the
        // northwest corner
        std::size_t const n_dims = length.size();
        NcVar const dim_var(getDimVar(dataset, var, n_dims - 2));
        auto const bounds = (dim_var.isNull()) ? getDimLength(var, n_dims - 2)
                                               : getBoundaries(dim_var);
        if (bounds.first > bounds.second)
        {
            std::size_t n_layers =
                (length.size() - time_offset == 3) ? length[n_dims - 3] : 1;
            flipRaster(data_vec, n_layers, length[n_dims - 1],
                       length[n_dims - 2]);
        }

        GeoLib::RasterHeader const header =
            createRasterHeader(dataset, var, dim_idx_map, length, time_offset);
        MeshLib::UseIntensityAs const useIntensity =
            MeshLib::UseIntensityAs::DATAVECTOR;
        if (output == OutputType::MULTIMESH)
        {
            mesh = MeshLib::RasterToMesh::convert(data_vec.data(), header,
                                                  elem_type, useIntensity,
                                                  var.getName());
            std::string const output_file_name(
                BaseLib::dropFileExtension(output_name) +
                getIterationString(i, time_bounds.second) + ".vtu");
            MeshLib::IO::VtuInterface vtu(mesh.get());
            vtu.writeToFile(output_file_name);
        }
        else if (output == OutputType::SINGLEMESH)
        {
            std::string array_name(var.getName());
            if (time_bounds.first != time_bounds.second)
                array_name.append(getIterationString(i, time_bounds.second));
            if (i == time_bounds.first)  // create persistent mesh
                mesh = MeshLib::RasterToMesh::convert(data_vec.data(), header,
                                                      elem_type, useIntensity,
                                                      array_name);
            else  // copy array to mesh
            {
                std::unique_ptr<MeshLib::Mesh> const temp(
                    MeshLib::RasterToMesh::convert(data_vec.data(), header,
                                                   elem_type, useIntensity,
                                                   array_name));
                MeshLib::PropertyVector<double> const* const vec =
                    temp->getProperties().getPropertyVector<double>(array_name);
                if (vec == nullptr)
                    return false;
                MeshLib::addPropertyToMesh<double>(
                    *mesh, array_name, MeshLib::MeshItemType::Cell, 1, *vec);
            }
            if (i == time_bounds.second)
            {
                MeshLib::IO::VtuInterface vtu(mesh.get());
                std::string const output_file_name =
                    (BaseLib::getFileExtension(output_name) == ".vtu")
                        ? output_name
                        : output_name + ".vtu";
                vtu.writeToFile(output_file_name);
            }
        }
        else  // OutputType::IMAGES
        {
            GeoLib::Raster const raster(header, data_vec.data(),
                                        data_vec.data() + header.n_cols *
                                                              header.n_rows *
                                                              header.n_depth);
            std::string const output_file_name(
                BaseLib::dropFileExtension(output_name) +
                getIterationString(i, time_bounds.second) + ".asc");
            FileIO::AsciiRasterInterface::writeRasterAsASC(raster,
                                                           output_file_name);
        }
    }
    return true;
}

int main(int argc, char* argv[])
{
    TCLAP::CmdLine cmd(
        "Converts NetCDF data into mesh file(s).\n\n "
        "OpenGeoSys-6 software, version " +
            GitInfoLib::GitInfo::ogs_version +
            ".\n"
            "Copyright (c) 2012-2021, OpenGeoSys Community "
            "(http://www.opengeosys.org)",
        ' ', GitInfoLib::GitInfo::ogs_version);

    TCLAP::ValueArg<int> arg_nodata(
        "n", "nodata",
        "explicitly specifies the no data value used in the dataset (usually "
        "it is not necessary to set this)",
        false, no_data_input, "integer specifying no data value");
    cmd.add(arg_nodata);

    std::vector<std::string> allowed_elems{"tri", "quad", "prism", "hex"};
    TCLAP::ValuesConstraint<std::string> allowed_elem_vals(allowed_elems);
    TCLAP::ValueArg<std::string> arg_elem_type(
        "e", "elem-type", "the element type used in the resulting OGS mesh",
        false, "", &allowed_elem_vals);
    cmd.add(arg_elem_type);

    TCLAP::SwitchArg arg_images(
        "", "images",
        "if set, all time steps will be written as ESRI image files (*.asc)");
    cmd.add(arg_images);

    TCLAP::SwitchArg arg_multi_files(
        "", "multi-file",
        "if set, each time step will be written to a separate mesh file");
    cmd.add(arg_multi_files);

    TCLAP::SwitchArg arg_single_file(
        "", "single-file",
        "if set, all time steps will be written to a single mesh file (with "
        "one scalar array per time step)");
    cmd.add(arg_single_file);

    TCLAP::ValueArg<std::size_t> arg_time_end(
        "", "timestep-last",
        "last time step to be extracted (only for time-dependent variables!)",
        false, 0, "integer specifying index of time step");
    cmd.add(arg_time_end);

    TCLAP::ValueArg<std::size_t> arg_time_start(
        "", "timestep-first",
        "first time step to be extracted (only for time-dependent variables!)",
        false, 0, "integer specifying index of time step");
    cmd.add(arg_time_start);

    std::vector<std::size_t> allowed_dims{0, 1, 2, 3};
    TCLAP::ValuesConstraint<std::size_t> allowed_dim_vals(allowed_dims);
    TCLAP::ValueArg<std::size_t> arg_dim3(
        "", "dim3",
        "index of third dimension (z/height/depth) for the selected variable",
        false, 0, &allowed_dim_vals);
    cmd.add(arg_dim3);

    TCLAP::ValueArg<std::size_t> arg_dim2(
        "", "dim2",
        "index of second dimension (y/latitude) for the selected "
        "variable",
        false, 0, &allowed_dim_vals);
    cmd.add(arg_dim2);

    TCLAP::ValueArg<std::size_t> arg_dim1(
        "", "dim1",
        "index of first dimension (x/longitude) for the selected variable",
        false, 0, &allowed_dim_vals);
    cmd.add(arg_dim1);

    TCLAP::ValueArg<std::size_t> arg_dim_time(
        "t", "time",
        "index of the time-dependent dimension for the selected variable",
        false, 0, &allowed_dim_vals);
    cmd.add(arg_dim_time);

    TCLAP::ValueArg<std::string> arg_varname(
        "v", "var", "variable included in the the netCDF file", false, "",
        "string containing the variable name");
    cmd.add(arg_varname);

    TCLAP::ValueArg<std::string> arg_output(
        "o", "output", "the OGS mesh output file", true, "",
        "string containing the path and file name");
    cmd.add(arg_output);
    TCLAP::ValueArg<std::string> arg_input(
        "i", "input", "the netCDF input file", true, "",
        "string containing the path and file name");
    cmd.add(arg_input);
    cmd.parse(argc, argv);

    NcFile dataset(arg_input.getValue().c_str(), NcFile::read);

    if (dataset.isNull())
    {
        ERR("Error opening file.");
        return -1;
    }

    std::size_t const mutex =
        arg_single_file.isSet() + arg_multi_files.isSet() + arg_images.isSet();
    if (mutex > 1)
    {
        ERR("Only one output format can be specified (single-file, multi-file, "
            "or images)");
        return EXIT_FAILURE;
    }

    std::cout << "OpenGeoSys NetCDF Converter\n";
    if (!arg_varname.isSet())
    {
        std::cout << "File " << arg_input.getValue()
                  << " loaded. Press ENTER to display available data arrays.\n";
        std::cin.ignore();
    }

    std::string const& output_name(arg_output.getValue());
    std::string const& var_name = (arg_varname.isSet())
                                      ? arg_varname.getValue()
                                      : arraySelectionLoop(dataset);
    NcVar const& var = dataset.getVar(var_name);
    if (var.isNull())
    {
        ERR("Variable \"{:s}\" not found in file.", arg_varname.getValue());
        return EXIT_FAILURE;
    }

    std::vector<std::size_t> dim_idx_map(var.getDimCount(), 0);
    bool is_time_dep(false);
    if (arg_dim1.isSet() && arg_dim2.isSet())
    {
        is_time_dep = arg_dim_time.isSet();
        if (!assignDimParams(var, dim_idx_map, arg_dim_time, arg_dim1, arg_dim2,
                             arg_dim3))
            return EXIT_FAILURE;
    }
    else
    {
        is_time_dep = dimensionSelectionLoop(var, dim_idx_map);
    }

    std::pair<std::size_t, std::size_t> time_bounds(0, 0);
    if (is_time_dep)
        time_bounds =
            (arg_time_start.isSet())
                ? assignTimeBounds(getDimVar(dataset, var, dim_idx_map[0]),
                                   arg_time_start, arg_time_end)
                : timestepSelectionLoop(var, dim_idx_map[0]);

    OutputType output = OutputType::INVALID;
    if (arg_images.isSet())
    {
        output = OutputType::IMAGES;
    }
    else if (arg_multi_files.isSet())
    {
        output = OutputType::MULTIMESH;
    }
    else if (arg_single_file.isSet() || !is_time_dep ||
             time_bounds.first == time_bounds.second)
    {
        output = OutputType::SINGLEMESH;
    }
    else
    {
        output = multFilesSelectionLoop(time_bounds);
    }

    std::size_t const temp_offset = (is_time_dep) ? 1 : 0;
    std::size_t const n_dims = (var.getDimCount());
    MeshLib::MeshElemType elem_type = MeshLib::MeshElemType::INVALID;

    if (output != OutputType::IMAGES)
    {
        elem_type = (arg_elem_type.isSet())
                        ? assignElemType(arg_elem_type)
                        : elemSelectionLoop(n_dims - temp_offset);
        if (elem_type == MeshLib::MeshElemType::INVALID)
            elemSelectionLoop(n_dims - temp_offset);
    }

    if (arg_nodata.isSet())
    {
        no_data_input = arg_nodata.getValue();
    }

    if (!convert(dataset, var, output_name, dim_idx_map, is_time_dep,
                 time_bounds, output, elem_type))
        return EXIT_FAILURE;

    std::cout << "Conversion finished successfully.\n";
    return EXIT_SUCCESS;
}
