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
        ERR("Input not valid.");
    }
    else if (error_id == 1)
    {
        ERR("Index not valid. Valid indices are in [0,%d].", max);
    }
    else if (error_id == 2)
    {
        showErrorMessage(0);
        std::cout << "Type \"info\" to display the available options again. "
                     "\"exit\" will exit the programme.\n";
    }
}

std::size_t parseInput(std::string const& request_str, std::size_t max_val,
                       bool has_info = false)
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

NcVar getDimVar(NcFile const& dataset, NcVar const& var, std::size_t const dim)
{
    NcDim const& dim_obj = var.getDim(dim);
    return dataset.getVar(dim_obj.getName());
}

std::vector<std::string> showArrays(NcFile const& dataset)
{
    std::size_t const n_vars(dataset.getDimCount());
    std::cout << "The NetCDF file contains the following " << n_vars
              << " arrays:\n\n";
    std::cout << "\tIndex\tArray Name\t#Dimensions\n";
    std::cout << "-------------------------------------------\n";
    auto const& names = dataset.getVars();
    std::vector<std::string> var_names;
    for (auto [name, var] : names)
    {
        std::cout << "\t" << var_names.size() << "\t" << name << "\t("
                  << var.getDimCount() << "D array)\n";
        var_names.push_back(name);
    }
    std::cout << "\n";
    return var_names;
}

bool showArraysDims(NcVar const& var)
{
    std::cout << "Data array \"" << var.getName()
              << "\" contains the following dimensions:\n";
    std::size_t const n_dims(var.getDimCount());
    for (std::size_t i = 0; i < n_dims; ++i)
        std::cout << "\t" << i << "\t" << var.getDim(i).getName() << "\t("
                  << var.getDim(i).getSize() << " values)\n";
    std::cout << "\n";
    return true;
}

std::pair<double, double> getBoundaries(NcVar const& var)
{
    if ((var.getDimCount()) == 1)
    {
        double start, end;
        std::size_t const size = var.getDim(0).getSize();
        var.getVar({0}, {1}, &start);
        var.getVar({size - 1}, {1}, &end);
        return std::make_pair(start, end);
    }
    return std::make_pair(0, 0);
}

MathLib::Point3d getOrigin(NcFile const& dataset, NcVar const& var,
                           std::array<std::size_t, 4> const& dim_idx_map,
                           bool is_time_dep)
{
    std::size_t const temp_offset = (is_time_dep) ? 1 : 0;
    MathLib::Point3d origin(std::array<double, 3>{{0, 0, 0}});
    std::size_t const n_dims = var.getDimCount();
    for (std::size_t i = temp_offset; i < n_dims; ++i)
    {
        NcVar const& dim = getDimVar(dataset, var, dim_idx_map[i]);
        auto const bounds = getBoundaries(dim);
        origin[i - temp_offset] =
            (bounds.first < bounds.second) ? bounds.first : bounds.second;
    }
    return origin;
}

void flipRaster(std::vector<double>& data, std::size_t width,
                std::size_t height)
{
    std::size_t const length(data.size());
    std::vector<double> tmp_vec(length);
    for (std::size_t i = 0; i < height; i++)
    {
        std::size_t const line_idx(length - (width * (i + 1)));
        for (std::size_t j = 0; j < width; j++)
        {
            tmp_vec.push_back(data[line_idx + j]);
        }
    }
    std::copy(tmp_vec.cbegin(), tmp_vec.cend(), data.begin());
}

bool canConvert(NcVar const& var)
{
    bool ret(var.getDimCount() < 2);
    if (ret)
        ERR("Only 2+ dimensional variables can be converted into OGS Meshes.\n");
    return !ret;
}

std::string arraySelectionLoop(NcFile const& dataset)
{
    std::vector<std::string> const& names = showArrays(dataset);
    std::size_t const idx =
        parseInput("Enter data array index: ", dataset.getVarCount(), true);

    if (static_cast<int>(idx) == dataset.getVarCount() ||
        !canConvert(dataset.getVar(names[idx])))
        return arraySelectionLoop(dataset);

    return names[idx];
}

bool dimensionSelectionLoop(NcVar const& var,
                            std::array<std::size_t, 4>& dim_idx_map)
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

        std::string const request_str("Enter ID for dimension " + std::to_string(i) +
                                " " + dim_comment[i - start_idx] + ": ");
        std::size_t const idx = parseInput(request_str, var.getDimCount(), true);

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

std::pair<std::size_t, std::size_t> timestepSelectionLoop(NcFile const& dataset,
                                                          NcVar const& var,
                                                          std::size_t dim_idx)
{
    std::size_t const n_time_steps = var.getDim(dim_idx).getSize();
    std::pair<std::size_t, std::size_t> bounds(
        std::numeric_limits<std::size_t>::max(),
        std::numeric_limits<std::size_t>::max());
    std::cout << "\nThe dataset contains " << n_time_steps << " time steps.\n";
    bounds.first =
        parseInput("Specify first time step to export: ", n_time_steps, false);
    while (bounds.first > bounds.second || bounds.second > n_time_steps)
        bounds.second = parseInput(
            "Specify last time step to export: ", n_time_steps, false);
    return bounds;
}

MeshLib::MeshElemType elemSelectionLoop(std::size_t const dim)
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
            else
                return MeshLib::MeshElemType::QUAD;
        }

        if (dim == 3)
        {
            if (type != "p" && type != "h" && type != "prism" &&
                type != "hex" && type != "hexahedron")
                continue;
            if (type == "p" || type == "prism")
                return MeshLib::MeshElemType::PRISM;
            else
                return MeshLib::MeshElemType::HEXAHEDRON;
        }
    }
    return t;
}

bool multFilesSelectionLoop(
    std::pair<std::size_t, std::size_t> const& time_bounds)
{
    std::size_t const n_time_steps(time_bounds.second - time_bounds.first + 1);
    std::cout << "\nThe selection includes " << n_time_steps
              << " time steps.\n";
    std::cout << "0. Save data in " << n_time_steps
              << " mesh files with one scalar array each.\n";
    std::cout << "1. Save data in one mesh file with " << n_time_steps
              << " scalar arrays.\n";
    std::size_t ret = parseInput("Select preferred method: ", 2, false);
    return (ret == 0) ? true : false;
}

std::string getIterationString(std::size_t i, std::size_t max)
{
    std::size_t const max_length(std::to_string(max).length());
    std::string const current_str(std::to_string(i));
    return std::string(max_length - current_str.length(), '0') + current_str;
}

double getResolution(NcFile const& dataset, NcVar const& var)
{
    std::size_t const dim_idx = var.getDimCount() - 1;
    auto const bounds = getBoundaries(getDimVar(dataset, var, dim_idx));
    return fabs(bounds.second - bounds.first) /
           static_cast<double>(var.getDim(dim_idx).getSize());
}

GeoLib::RasterHeader createRasterHeader(
    NcFile const& dataset, NcVar const& var,
    std::array<std::size_t, 4> const& dim_idx_map,
    std::vector<std::size_t> const& length, bool is_time_dep)
{
    MathLib::Point3d const origin = getOrigin(dataset, var, dim_idx_map, is_time_dep);
    double const res = getResolution(dataset, var);
    std::size_t n_dims = var.getDimCount();
    std::size_t temp_offset = (is_time_dep) ? 1 : 0;
    std::size_t z_length =
        (n_dims - temp_offset == 3) ? length[dim_idx_map.back()] : 1;
    return {length[dim_idx_map[0 + temp_offset]],
            length[dim_idx_map[1 + temp_offset]],
            z_length, origin, res, -9999};
}

std::size_t getLength(NcVar const& var, bool const is_time_dep,
                      std::vector<std::size_t>& length)
{
    std::size_t const temp_offset = (is_time_dep) ? 1 : 0;
    std::size_t const n_dims = (var.getDimCount());
    length.resize(4, 1);
    std::size_t array_length = 1;
    for (std::size_t i = temp_offset; i < n_dims; ++i)
    {
        length[i] = var.getDim(i).getSize();
        array_length *= length[i];
    }
    return array_length;
}

std::vector<double> getData(NcFile const& dataset, NcVar const& var,
                            std::size_t total_length,
                            std::size_t const time_step,
                            std::vector<std::size_t> const& length)
{
    std::size_t const n_dims(var.getDimCount());
    std::vector<std::size_t> offset(n_dims, 0);
    offset[0] = time_step;
    std::vector<double> data_vec(total_length, 0);
    var.getVar(offset, length, data_vec.data());
    std::replace_if(data_vec.begin(), data_vec.end(),
                    [](double& x) { return x <= -9999; }, -9999);

    // reverse lines in vertical direction if the original file has its origin
    // in the northwest corner
    auto const bounds = getBoundaries(getDimVar(dataset, var, n_dims - 1));
    if (bounds.first > bounds.second)
        flipRaster(data_vec, length[n_dims - 2], length[n_dims - 1]);
    return data_vec;
}

bool assignDimParams(NcVar const& var,
                     std::array<std::size_t, 4>& dim_idx_map,
                     TCLAP::ValueArg<std::size_t>& arg_dim_time,
                     TCLAP::ValueArg<std::size_t>& arg_dim1,
                     TCLAP::ValueArg<std::size_t>& arg_dim2,
                     TCLAP::ValueArg<std::size_t>& arg_dim3)
{
    if (arg_dim3.isSet() && !arg_dim2.isSet())
    {
        ERR("Parameter for dim2 is not set. Ignoring dimension parameters.");
        return false;
    }
    if (!arg_dim1.isSet() || !arg_dim2.isSet())
    {
        ERR("dim1 and dim2 need to be set for extracting mesh.");
        return false;
    }

    std::size_t const n_dims = var.getDimCount();
    if (arg_dim_time.getValue() >= n_dims || arg_dim1.getValue() >= n_dims ||
        arg_dim2.getValue() >= n_dims || arg_dim3.getValue() >= n_dims)
    {
        ERR("Maximum allowed dimension for variable \"%s\" is %d.", var.getName().c_str(), n_dims-1);
        return false;
    }

    bool const is_time_dep = arg_dim_time.isSet();
    if (is_time_dep)
        dim_idx_map[0] = arg_dim_time.getValue();
    std::size_t const temp_offset = (is_time_dep) ? 1 : 0;
    dim_idx_map[0 + temp_offset] = arg_dim1.getValue();
    dim_idx_map[1 + temp_offset] = arg_dim2.getValue();
    if (n_dims == (3 + temp_offset))
        dim_idx_map[2 + temp_offset] = arg_dim3.getValue();

    return true;
}

std::pair<std::size_t, std::size_t> assignTimeBounds(NcVar const& var,
                     TCLAP::ValueArg<std::size_t>& arg_time_start,
                     TCLAP::ValueArg<std::size_t>& arg_time_end)
{
    auto const bounds = getBoundaries(var);
    if (arg_time_start.getValue() > bounds.second)
    {
        ERR("Start time step larger than total number of time steps. Resetting to 0.");
        arg_time_start.reset();
    }

    if (!arg_time_end.isSet())
        return {arg_time_start.getValue(), arg_time_start.getValue()};

    if (arg_time_end.getValue() > bounds.second)
    {
        ERR("End time step larger than total number of time steps. Resetting to starting time step");
        return {arg_time_start.getValue(), arg_time_start.getValue()};
    }

    if (arg_time_end.getValue() < arg_time_start.getValue())
    {
        ERR("End time step larger than starting time step. Swapping values");
        return {arg_time_end.getValue(), arg_time_start.getValue()};
    }

    return {arg_time_start.getValue(), arg_time_end.getValue()};
}

MeshLib::MeshElemType assignElemType(TCLAP::ValueArg<std::string>& arg_elem_type)
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

bool convert(NcFile const& dataset, NcVar const& var,
             std::string const& output_name,
             std::array<std::size_t, 4> const& dim_idx_map,
             bool const is_time_dep,
             std::pair<std::size_t, std::size_t> const& time_bounds,
             bool const use_single_file, MeshLib::MeshElemType const elem_type)
{
    std::unique_ptr<MeshLib::Mesh> mesh;
    std::vector<std::size_t> length;
    std::size_t const array_length = getLength(var, is_time_dep, length);
    for (std::size_t i = time_bounds.first; i <= time_bounds.second; ++i)
    {
        std::cout << "Converting time step " << i << "...\n";
        std::vector<double> const data_vec =
            getData(dataset, var, array_length, i, length);

        GeoLib::RasterHeader const header =
            createRasterHeader(dataset, var, dim_idx_map, length, is_time_dep);
        MeshLib::UseIntensityAs const useIntensity =
            MeshLib::UseIntensityAs::DATAVECTOR;
        if (!use_single_file)
        {
            mesh.reset(MeshLib::RasterToMesh::convert(
                data_vec.data(), header, elem_type, useIntensity,
                var.getName()));
            std::string const output_file_name(
                BaseLib::dropFileExtension(output_name) +
                getIterationString(i, time_bounds.second) + ".vtu");
            MeshLib::IO::VtuInterface vtu(mesh.get());
            vtu.writeToFile(output_file_name);
        }
        else
        {
            std::string array_name(var.getName());
            if (time_bounds.first != time_bounds.second)
                array_name.append(getIterationString(i, time_bounds.second));
            if (i == time_bounds.first)  // create persistent mesh
                mesh.reset(MeshLib::RasterToMesh::convert(
                    data_vec.data(), header, elem_type, useIntensity,
                    array_name));
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
                vtu.writeToFile(output_name);
            }
        }
    }
    return true;
}

int main(int argc, char* argv[])
{
    ApplicationsLib::LogogSetup logog_setup;

    TCLAP::CmdLine cmd(
        "Converts NetCDF data into mesh file(s).\n\n "
        "OpenGeoSys-6 software, version " +
            GitInfoLib::GitInfo::ogs_version +
            ".\n"
            "Copyright (c) 2012-2019, OpenGeoSys Community "
            "(http://www.opengeosys.org)",
        ' ', GitInfoLib::GitInfo::ogs_version);

    std::vector<std::string> allowed_elems{"tri", "quad", "prism", "hex"};
    TCLAP::ValuesConstraint<std::string> allowed_elem_vals(allowed_elems);
    TCLAP::ValueArg<std::string> arg_elem_type(
        "e", "elem-type", "the element type used in the resulting OGS mesh",
        false, "", &allowed_elem_vals);
    cmd.add(arg_elem_type);

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
        ERR("Variable \"%s\" not found in file.", arg_varname.getValue().c_str());
        return EXIT_FAILURE;
    }

    std::array<std::size_t, 4> dim_idx_map;
    bool is_time_dep (false);
    if (arg_dim1.isSet() && arg_dim2.isSet())
    {
        is_time_dep = arg_dim_time.isSet();
        if (!assignDimParams(var, dim_idx_map, arg_dim_time, arg_dim1, arg_dim2, arg_dim3))
            return EXIT_FAILURE;
    }
    else
    {
        is_time_dep = dimensionSelectionLoop(var, dim_idx_map);
    }

    std::pair<std::size_t, std::size_t> time_bounds(0, 0);
    if (is_time_dep)
        time_bounds = (arg_time_start.isSet())
                ? assignTimeBounds(getDimVar(dataset, var, dim_idx_map[0]),
                                   arg_time_start, arg_time_end)
                : timestepSelectionLoop(dataset, var, dim_idx_map[0]);

    bool use_single_file(true);
    if (arg_time_start.isSet())
    {
        use_single_file = arg_single_file.isSet();
    }
    else
    {
        if (is_time_dep && time_bounds.first != time_bounds.second)
            use_single_file = multFilesSelectionLoop(time_bounds);
    }

    std::size_t const temp_offset = (is_time_dep) ? 1 : 0;
    std::size_t const n_dims = (var.getDimCount());
    MeshLib::MeshElemType const elem_type = (arg_elem_type.isSet())
                            ? assignElemType(arg_elem_type)
                            : elemSelectionLoop(n_dims - temp_offset);

    if (!convert(dataset, var, output_name, dim_idx_map, is_time_dep,
                 time_bounds, use_single_file, elem_type))
        return EXIT_FAILURE;

    std::cout << "Conversion finished successfully.\n";
    return EXIT_SUCCESS;
}
