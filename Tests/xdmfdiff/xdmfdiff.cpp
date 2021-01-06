/**
 * \file
 * \copyright
 * Copyright (c) 2015-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <tclap/CmdLine.h>

#include <Xdmf.hpp>
#include <XdmfDomain.hpp>
#include <XdmfGridCollection.hpp>
#include <XdmfReader.hpp>
#include <XdmfUnstructuredGrid.hpp>
#include <boost/iterator/zip_iterator.hpp>
#include <boost/range.hpp>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <ios>
#include <memory>
#include <sstream>
#include <tuple>
// See https://stackoverflow.com/a/8513803/2706707
template <typename... Containers>
auto zip(Containers&&... containers) -> boost::iterator_range<
    boost::zip_iterator<decltype(boost::make_tuple(std::begin(containers)...))>>
{
    auto zip_begin =
        boost::make_zip_iterator(boost::make_tuple(std::begin(containers)...));
    auto zip_end =
        boost::make_zip_iterator(boost::make_tuple(std::end(containers)...));
    return boost::make_iterator_range(zip_begin, zip_end);
}

struct Args
{
    bool const quiet;
    bool const verbose;
    bool const meshcheck;
    double const abs_err_thr;
    double const rel_err_thr;
    std::string const xdmf_input_a;
    std::string const xdmf_input_b;
    std::string const data_array_a;
    std::string const data_array_b;
    unsigned int timestep_a;
    unsigned int timestep_b;
};

auto parseCommandLine(int argc, char* argv[]) -> Args
{
    TCLAP::CmdLine cmd(
        "XdmfDiff software.\n"
        "Copyright (c) 2015-2020, OpenGeoSys Community "
        "(http://www.opengeosys.org) "
        "Distributed under a Modified BSD License. "
        "See accompanying file LICENSE.txt or "
        "http://www.opengeosys.org/project/license",
        ' ',
        "0.1");

    TCLAP::UnlabeledValueArg<std::string> xdmf_input_a_arg(
        "input-file-a", "Path to the Xdmf input file.", true, "", "XDMF FILE");
    cmd.add(xdmf_input_a_arg);

    TCLAP::UnlabeledValueArg<std::string> xdmf_input_b_arg(
        "input-file-b",
        "Path to the second XDMF input file.",
        false,
        "",
        "XDMF FILE");
    cmd.add(xdmf_input_b_arg);

    TCLAP::ValueArg<std::string> data_array_a_arg(
        "a", "first_data_array", "First data array name for comparison", true,
        "", "NAME");

    TCLAP::ValueArg<std::string> data_array_b_arg(
        "b",
        "second_data_array",
        "Second data array name for comparison",
        false,
        "",
        "NAME");
    cmd.add(data_array_b_arg);

    TCLAP::ValueArg<unsigned int> time_a_arg(
        "", "timestep-a", "First data time step index (positive integer)",
        false, 0, "TIMESTEP");
    cmd.add(time_a_arg);

    TCLAP::ValueArg<unsigned int> time_b_arg(
        "", "timestep-b", "Second data time step index (positive integer)",
        false, 0, "TIMESTEP");
    cmd.add(time_b_arg);

    TCLAP::SwitchArg meshcheck_arg(
        "m", "mesh_check", "Compare mesh geometries using absolute tolerance.");
    cmd.xorAdd(data_array_a_arg, meshcheck_arg);

    TCLAP::SwitchArg quiet_arg("q", "quiet", "Suppress all but error output.");
    cmd.add(quiet_arg);

    TCLAP::SwitchArg verbose_arg("v", "verbose",
                                 "Also print which values differ.");
    cmd.add(verbose_arg);

    auto const double_eps_string =
        std::to_string(std::numeric_limits<double>::epsilon());

    TCLAP::ValueArg<double> abs_err_thr_arg(
        "",
        "abs",
        "Tolerance for the absolute error in the maximum norm (" +
            double_eps_string + ")",
        false,
        std::numeric_limits<double>::epsilon(),
        "FLOAT");
    cmd.add(abs_err_thr_arg);

    TCLAP::ValueArg<double> rel_err_thr_arg(
        "",
        "rel",
        "Tolerance for the componentwise relative error (" + double_eps_string +
            ")",
        false,
        std::numeric_limits<double>::epsilon(),
        "FLOAT");
    cmd.add(rel_err_thr_arg);

    cmd.parse(argc, argv);

    return Args{quiet_arg.getValue(),        verbose_arg.getValue(),
                meshcheck_arg.getValue(),    abs_err_thr_arg.getValue(),
                rel_err_thr_arg.getValue(),  xdmf_input_a_arg.getValue(),
                xdmf_input_b_arg.getValue(), data_array_a_arg.getValue(),
                data_array_b_arg.getValue(), time_a_arg.getValue(),
                time_b_arg.getValue()};
}

struct Grid
{
    std::vector<double> point_values;
    std::vector<int> topology_values;
    std::vector<double> attribute_values;
};

std::unique_ptr<Grid> readMesh(std::string const& filename,
                               unsigned int timestep,
                               std::string const& attribute_name)
{
    if (filename.empty())
    {
        return nullptr;
    }

    if (std::filesystem::path(filename).extension().string() != ".xdmf")
    {
        std::cerr << "Error: Expected a file with .xdmf extension."
                  << "File '" << filename << "' not read.";
        return nullptr;
    }
    auto const xreader = XdmfReader::New();
    auto const domain =
        shared_dynamic_cast<XdmfDomain>(xreader->read(filename));
    auto const gridcollection = domain->getGridCollection("Collection");
    auto const ungrid = gridcollection->getUnstructuredGrid(timestep);
    auto const geometry = ungrid->getGeometry();
    geometry->read();
    int const size_points = geometry->getSize();
    std::vector<double> points(size_points);
    geometry->getValues(0, points.data(), size_points);

    auto const topology = ungrid->getTopology();
    topology->read();
    int const size_topology = topology->getSize();
    std::vector<int> topology_values(size_topology);
    topology->getValues(0, topology_values.data(), size_topology);

    auto const attribute = ungrid->getAttribute(attribute_name);
    attribute->read();
    int const attribute_size = attribute->getSize();
    std::vector<double> attribute_values(attribute_size);
    attribute->getValues(0, attribute_values.data(), attribute_size);

    return std::make_unique<Grid>(
        Grid{points, topology_values, attribute_values});
}

std::tuple<std::unique_ptr<Grid>, std::unique_ptr<Grid>> readMeshes(
    std::string const& file_a_name,
    std::string const& file_b_name,
    unsigned int const& timestep_a,
    unsigned int const& timestep_b,
    std::string const& attribute_a_name,
    std::string const& attribute_b_name)
{
    return {readMesh(file_a_name, timestep_a, attribute_a_name),
            readMesh(file_b_name, timestep_b, attribute_b_name)};
}

bool compareCellTopology(std::vector<int> const& cells_a,
                         std::vector<int> const& cells_b)
{
    if (cells_a.size() != cells_b.size())
    {
        std::cerr << "Number of cells in the first mesh is " << cells_a.size()
                  << " and differs from the number of cells in the second "
                     "mesh, which is "
                  << cells_b.size() << "\n";
        return false;
    }

    for (std::size_t i = 0; i < cells_a.size(); ++i)
    {
        if (cells_a[i] != cells_b[i])
        {
            std::cerr << "Cell on position " << i << "\n"
                      << "in first mesh has celltype: " << cells_a[i] << "\n"
                      << "in second mesh has celltype: " << cells_b[i] << ".\n";

            return false;
        }
    }
    return true;
}

bool compareValues(std::vector<double> const& points_a,
                   std::vector<double> const& points_b, double const eps)
{
    if (points_a.size() != points_b.size())
    {
        std::cerr << "Number of points in the first mesh is " << points_a.size()
                  << " and differst from the number of point in the second "
                     "mesh, which is "
                  << points_b.size() << "\n";
        return false;
    }

    for (std::size_t i = 0; i < points_a.size(); ++i)
    {
        if (double const distance = std::abs(points_a[i] - points_b[i]);
            distance >= eps)
        {
            std::cerr << "Point on position " << i << "\n"
                      << "in first mesh has value: " << points_a[i] << "\n"
                      << "in second mesh has value: " << points_b[i] << ".\n"
                      << "The distance is: " << distance << "\n";
            return false;
        }
        return true;
    }
    return true;
}

int main(int argc, char* argv[])
{
    auto const digits10 = std::numeric_limits<double>::digits10;
    auto const args = parseCommandLine(argc, argv);

    // Setup the standard output and error stream numerical formats.
    std::cout << std::scientific << std::setprecision(digits10);
    std::cerr << std::scientific << std::setprecision(digits10);

    auto const [mesh_a, mesh_b] =
        readMeshes(args.xdmf_input_a, args.xdmf_input_b, args.timestep_a,
                   args.timestep_b, args.data_array_a, args.data_array_b);

    if (args.xdmf_input_a == args.xdmf_input_b)
    {
        std::cout << "Will not compare meshes from same input file.\n";
        return EXIT_SUCCESS;
    }
    if (!compareValues(mesh_a->point_values, mesh_b->point_values,
                       args.abs_err_thr))
    {
        std::cerr << "Error in mesh points' comparison occured.\n";
        return EXIT_FAILURE;
    }

    if (!compareCellTopology(mesh_a->topology_values, mesh_b->topology_values))
    {
        std::cerr << "Error in cells' topology comparison occured.\n";
        return EXIT_FAILURE;
    }

    if (!compareValues(mesh_a->attribute_values,
                       mesh_b->attribute_values,
                       args.abs_err_thr))
    {
        std::cerr << "Error in mesh attribute ( " << args.data_array_a << ") "
                  << "of first mesh and ( " << args.data_array_b << ") "
                  << "of second mesh occured.\n";
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}