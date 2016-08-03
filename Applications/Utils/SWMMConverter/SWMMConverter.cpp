
#include <tclap/CmdLine.h>

#include "SWMMInterface.h"
#include "ThirdParty/SWMMInterface/swmm5_iface.h"

#include "Applications/ApplicationsLib/LogogSetup.h"

#include "GeoLib/GEOObjects.h"
#include "GeoLib/IO/XmlIO/Boost/BoostXmlGmlInterface.h"

#include "MeshLib/Mesh.h"
#include "MeshLib/IO/VtkIO/VtuInterface.h"

int main(int argc, char *argv[])
{
    ApplicationsLib::LogogSetup setup;

    TCLAP::CmdLine cmd
        ("Read files for the Storm Water Management Model (SWMM) and converts them to OGS.", ' ', "0.1");
    TCLAP::ValueArg<std::string> mesh_output_arg
        ("m","mesh-output", "mesh output file (*.vtu)", false, "", "mesh output file");
    cmd.add(mesh_output_arg);
    TCLAP::ValueArg<std::string> geo_output_arg
        ("g","geo-output", "geometry output file (*.gml)", false, "", "geometry output file");
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

    SwmmInterface swmm(swmm_input_arg.getValue());

    if (geo_output_arg.isSet())
    {
        GeoLib::GEOObjects geo_objects;
        std::string const inp_file_name (swmm_input_arg.getValue() + ".inp");
        if (swmm.SwmmInputToGeometry(inp_file_name, geo_objects) == 0)
        {
            GeoLib::IO::BoostXmlGmlInterface xml(geo_objects);
            xml.setNameForExport(swmm_input_arg.getValue());
            xml.writeToFile(geo_output_arg.getValue());
        }
    }

    if (mesh_output_arg.isSet())
    {
        MeshLib::Mesh* mesh = swmm.SwmmInputToLineMesh();
        MeshLib::IO::VtuInterface vtkIO(mesh, 0, false);
        vtkIO.writeToFile(mesh_output_arg.getValue());
    }

    return 0;
}
