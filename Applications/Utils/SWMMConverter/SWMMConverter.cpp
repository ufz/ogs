/**
 * @copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <tclap/CmdLine.h>

#include "Applications/FileIO/SWMM/SWMMInterface.h"

#include "Applications/ApplicationsLib/LogogSetup.h"

#include "BaseLib/StringTools.h"

#include "GeoLib/GEOObjects.h"
#include "GeoLib/IO/XmlIO/Boost/BoostXmlGmlInterface.h"

#include "MeshLib/Mesh.h"
#include "MeshLib/Properties.h"
#include "MeshLib/IO/VtkIO/VtuInterface.h"

#include "Applications/FileIO/CsvInterface.h"

int main(int argc, char *argv[])
{
    ApplicationsLib::LogogSetup setup;

    TCLAP::CmdLine cmd
        ("Read files for the Storm Water Management Model (SWMM) and converts them to OGS.", ' ', "0.1");
    TCLAP::ValueArg<std::string> mesh_output_arg
        ("m","mesh", "mesh output file (*.vtu)", false, "", "mesh output file");
    cmd.add(mesh_output_arg);
    TCLAP::ValueArg<std::string> geo_output_arg
        ("g","geo", "geometry output file (*.gml)", false, "", "geometry output file");
    cmd.add(geo_output_arg);
    TCLAP::ValueArg<std::string> swmm_input_arg
        ("i","input", "SWMM input file (*.inp)", true, "", "input file");
    cmd.add(swmm_input_arg);
    cmd.parse( argc, argv );

    if (!(geo_output_arg.isSet() || mesh_output_arg.isSet()))
    {
        ERR ("No output format given. Please specify OGS geometry or mesh output file.");
        return -1;
    }

    if (geo_output_arg.isSet())
    {
        GeoLib::GEOObjects geo_objects;
        if (!FileIO::SwmmInterface::convertSwmmInputToGeometry(swmm_input_arg.getValue(), geo_objects, true))
            return -1;

        GeoLib::IO::BoostXmlGmlInterface xml(geo_objects);
        xml.setNameForExport(BaseLib::extractBaseNameWithoutExtension(swmm_input_arg.getValue()));
        xml.writeToFile(geo_output_arg.getValue());
        return 0;
    }

    std::unique_ptr<FileIO::SwmmInterface> swmm = nullptr;
    std::unique_ptr<MeshLib::Mesh> mesh = nullptr;
    if (mesh_output_arg.isSet())
    {
        swmm.reset(FileIO::SwmmInterface::create(swmm_input_arg.getValue()));
        if (swmm == nullptr)
            return -1;

        std::size_t a = swmm->getNumberOfObjects(FileIO::SwmmObject::SUBCATCHMENT);
        mesh.reset(swmm->getMesh());
        MeshLib::IO::VtuInterface vtkIO(mesh.get(), 0, false);
        vtkIO.writeToFile(mesh_output_arg.getValue());
    }

    std::cout << "Simulation time steps: " << swmm->getNumberOfTimeSteps() << std::endl;


    // Example code: Writing node information to csv file
    FileIO::CsvInterface csv;
    csv.addIndexVectorForWriting(swmm->getNumberOfObjects(FileIO::SwmmObject::NODE));
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;
    swmm->getNodeCoordinateVectors(x,y,z);
    csv.addVectorForWriting("x", x);
    csv.addVectorForWriting("y", y);
    csv.addVectorForWriting("z", z);
    for (std::size_t i=0; i<9; ++i)
    {
        std::vector<double> data_vec = swmm->getArrayAtTimeStep(FileIO::SwmmObject::NODE, 10, i);
        csv.addVectorForWriting(swmm->getArrayName(FileIO::SwmmObject::NODE, i), data_vec);
    }
    csv.writeToFile("d:/csvtest.csv");
    std::cout << "csv written" << std::endl;
    std::cin.ignore();

    for (std::size_t i=0; i<swmm->getNumberOfParameters(FileIO::SwmmObject::NODE); ++i)
        std::cout << i << "\t" << swmm->getArrayName(FileIO::SwmmObject::NODE, i) << std::endl;
    std::cin.ignore();

    // Example code: Add simulated parameter to mesh for each timestep and write result
    std::size_t n_time_steps (swmm->getNumberOfTimeSteps());
    for (std::size_t i=0; i<n_time_steps; i++)
    {
        FileIO::SwmmObject type = FileIO::SwmmObject::NODE;
        for (std::size_t j=6; j<9; ++j)
        {
            std::string vec_name (swmm->getArrayName(type, j));
            if (vec_name.empty())
                return -10;
            std::vector<double> data_vec = swmm->getArrayAtTimeStep(type, i, j);
            if (data_vec.empty())
                return -20;
            bool done = swmm->addResultsToMesh(*mesh, type, vec_name, data_vec);
            if (!done)
                return -30;
        }

        MeshLib::IO::VtuInterface vtkio(mesh.get(), 0, false);
        std::string name ("d:/swmmresults" + BaseLib::tostring(i) + ".vtu");
        vtkio.writeToFile(name);
        mesh->getProperties().removePropertyVector("P");
        mesh->getProperties().removePropertyVector("NH4");
        mesh->getProperties().removePropertyVector("CSB");
    }

    return 0;
}
