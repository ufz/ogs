/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <tclap/CmdLine.h>

#ifdef USE_PETSC
#include <mpi.h>
#endif

#include <array>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <string>

#include "BaseLib/FileTools.h"
#include "BaseLib/Logging.h"
#include "BaseLib/MemWatch.h"
#include "BaseLib/RunTime.h"
#include "BaseLib/StringTools.h"
#include "InfoLib/GitInfo.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/IO/XDMF/XdmfHdfWriter.h"
#include "MeshLib/IO/readMeshFromFile.h"
#include "MeshLib/Mesh.h"

// TODO (naumov) use std::filesystem::path
std::vector<std::pair<double, std::string>> readPvd(
    std::string const& pvd_filename)
{
    DBUG("Start reading the PVD file {:s}", pvd_filename);
    boost::property_tree::ptree pvd;
    read_xml(pvd_filename, pvd,
             boost::property_tree::xml_parser::trim_whitespace);

    std::vector<std::pair<double, std::string>> timeseries;
    for (auto const& dataset : pvd.get_child("VTKFile.Collection"))
    {
        if (dataset.first != "DataSet")
        {
            OGS_FATAL("Expected DataSet tag but got '{:s}'", dataset.first);
        }

        auto const time = dataset.second.get<double>("<xmlattr>.timestep");
        auto const file = dataset.second.get<std::string>("<xmlattr>.file");
        timeseries.emplace_back(time, file);
    }
    DBUG("Finished reading the PVD file {:s}", pvd_filename);

    return timeseries;
}

template <typename T>
bool copyPropertyVector(MeshLib::Properties const& properties,
                        MeshLib::PropertyVectorBase* destination_pv)
{
    if (!dynamic_cast<MeshLib::PropertyVector<T>*>(destination_pv))
    {
        return false;
    }

    auto const* pv = properties.getPropertyVector<T>(
        destination_pv->getPropertyName(), destination_pv->getMeshItemType(),
        destination_pv->getNumberOfGlobalComponents());

    assert(pv != nullptr);

    std::copy(
        std::begin(*pv), std::end(*pv),
        std::begin(dynamic_cast<MeshLib::PropertyVector<T>&>(*destination_pv)));

    return true;
}

void testIfMeshesAreEqual(MeshLib::Mesh const& mesh,
                          MeshLib::Mesh const& main_mesh,
                          std::string error_message)
{
    if (mesh.getDimension() != main_mesh.getDimension())
    {
        OGS_FATAL(
            "{:s}The mesh has dimension {} which is different from the "
            "dimension {} of in the main mesh.",
            error_message, mesh.getDimension(), main_mesh.getDimension());
    }
    if (mesh.getNumberOfElements() != main_mesh.getNumberOfElements())
    {
        OGS_FATAL(
            "{:s}The mesh has {} elements which is different from the number "
            "of elements ({}) in the main mesh.",
            error_message, mesh.getNumberOfElements(),
            main_mesh.getNumberOfElements());
    }
    if (mesh.getNumberOfNodes() != main_mesh.getNumberOfNodes())
    {
        OGS_FATAL(
            "{:s}The mesh has {} nodes which is different from the number of "
            "nodes ({}) in the main mesh.",
            error_message, mesh.getNumberOfNodes(),
            main_mesh.getNumberOfNodes());
    }
}

int main(int argc, char* argv[])
{
    TCLAP::CmdLine cmd(
        "Converts a time series from PVD to XDMF format.\n\n"
        "OpenGeoSys-6 software, version " +
            GitInfoLib::GitInfo::ogs_version +
            ".\n"
            "Copyright (c) 2012-2024, OpenGeoSys Community "
            "(http://www.opengeosys.org)",
        ' ', GitInfoLib::GitInfo::ogs_version);

    TCLAP::ValueArg<std::string> log_level_arg(
        "l", "log-level",
        "the verbosity of logging messages: none, error, warn, info, "
        "debug, "
        "all",
        false,
#ifdef NDEBUG
        "info",
#else
        "all",
#endif
        "LOG_LEVEL");
    cmd.add(log_level_arg);

    TCLAP::UnlabeledValueArg<std::string> pvd_file_arg("pvd-file", "pvd file",
                                                       true, "", "file");
    cmd.add(pvd_file_arg);

    TCLAP::ValueArg<std::string> outdir_arg("o", "output-directory",
                                            "the output directory to write to",
                                            false, ".", "PATH");
    cmd.add(outdir_arg);

    cmd.parse(argc, argv);
#ifdef USE_PETSC
    MPI_Init(&argc, &argv);
#endif
    BaseLib::setConsoleLogLevel(log_level_arg.getValue());

    auto const pvd_file_dir = BaseLib::extractPath(pvd_file_arg.getValue());

    auto const timeseries = readPvd(pvd_file_arg.getValue());

    if (timeseries.empty())
    {
        OGS_FATAL("Empty time series.");
    }

    // Initialized from the first timeseries entry and data is updated for each
    // subsequent time step.
    std::unique_ptr<MeshLib::Mesh> main_mesh;

    std::filesystem::path const output_file{
        BaseLib::extractBaseNameWithoutExtension(pvd_file_arg.getValue()) +
        ".xdmf"};
    std::filesystem::path output_file_path{outdir_arg.getValue()};
    if (outdir_arg.getValue() != "")
    {
        output_file_path =
            BaseLib::joinPaths(outdir_arg.getValue(), output_file.string());
    }
    std::set<std::string> variable_output_names;
    std::unique_ptr<MeshLib::IO::XdmfHdfWriter> mesh_xdmf_hdf_writer;
    // read first file in the time series; it is determining variables.
    {
        auto [time, filename] = timeseries[0];
        DBUG("{} - {}", time, filename);

        main_mesh.reset(MeshLib::IO::readMeshFromFile(
            BaseLib::joinPaths(pvd_file_dir, filename)));
        if (main_mesh == nullptr)
        {
            OGS_FATAL("Could not read mesh from '{:s}'.",
                      BaseLib::joinPaths(pvd_file_dir, filename));
        }

        for (auto const& p : main_mesh->getProperties())
        {
            variable_output_names.insert(std::string(p.first));
        }
        mesh_xdmf_hdf_writer = std::make_unique<MeshLib::IO::XdmfHdfWriter>(
            std::vector{std::cref(*main_mesh)}, output_file_path,
            0 /*timestep*/, time, variable_output_names,
            true /*output_file.compression*/, 1 /*output_file.n_files*/,
            1048576 /*chunk_size_bytes*/);
    }

    for (std::size_t timestep = 1; timestep < timeseries.size(); ++timestep)
    {
        auto [time, filename] = timeseries[timestep];
        DBUG("{} - {}", time, filename);

        std::unique_ptr<MeshLib::Mesh> mesh{MeshLib::IO::readMeshFromFile(
            BaseLib::joinPaths(pvd_file_dir, filename))};
        if (mesh == nullptr)
        {
            OGS_FATAL("Could not read mesh from '{:s}'.",
                      BaseLib::joinPaths(pvd_file_dir, filename));
        }

        testIfMeshesAreEqual(
            *mesh, *main_mesh,
            fmt::format("Error in comparison of mesh from file '{:s}' to the "
                        "main mesh:\n",
                        BaseLib::joinPaths(pvd_file_dir, filename)));

        // We have to copy the values because the xdmf writer remembers the data
        // pointers. Therefore replacing the properties will not work.
        for (auto& [name, pv] : main_mesh->getProperties())
        {
            if (copyPropertyVector<double>(mesh->getProperties(), pv) ||
                copyPropertyVector<float>(mesh->getProperties(), pv) ||
                copyPropertyVector<int>(mesh->getProperties(), pv) ||
                copyPropertyVector<long>(mesh->getProperties(), pv) ||
                copyPropertyVector<unsigned>(mesh->getProperties(), pv) ||
                copyPropertyVector<unsigned long>(mesh->getProperties(), pv) ||
                copyPropertyVector<std::size_t>(mesh->getProperties(), pv) ||
                copyPropertyVector<char>(mesh->getProperties(), pv) ||
                copyPropertyVector<unsigned char>(mesh->getProperties(), pv))
            {
                continue;
            }
            OGS_FATAL(
                "Could not copy property vector '{:s}' from '{:s}' mesh to the "
                "main_mesh.",
                name, BaseLib::joinPaths(pvd_file_dir, filename));
        }

        mesh_xdmf_hdf_writer->writeStep(time);
    }
#ifdef USE_PETSC
    MPI_Finalize();
#endif
    return EXIT_SUCCESS;
}
