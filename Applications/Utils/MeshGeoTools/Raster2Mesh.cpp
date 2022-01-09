/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <memory>
#include <string>

// ThirdParty
#include <tclap/CmdLine.h>

#include "Applications/FileIO/AsciiRasterInterface.h"
#include "GeoLib/Raster.h"
#include "InfoLib/GitInfo.h"
#include "MeshLib/IO/VtkIO/VtuInterface.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshEnums.h"
#include "MeshLib/MeshGenerators/RasterToMesh.h"

int main(int argc, char* argv[])
{
    TCLAP::CmdLine cmd(
        "Converts an ASCII raster file (*.asc) into a 2D triangle- or "
        "quad-mesh. Pixel values can be interpreted as elevation of mesh nodes "
        "or as scalar values for mesh elements.\nIt is highly recommended to "
        "create triangle meshes when interpreting pixel values as elevation, "
        "as it is likely that the resulting mesh will otherwise contain "
        "deformed (i.e. non-planar) quad-elements.\n\n"
        "OpenGeoSys-6 software, version " +
            GitInfoLib::GitInfo::ogs_version +
            ".\n"
            "Copyright (c) 2012-2022, OpenGeoSys Community "
            "(http://www.opengeosys.org)",
        ' ', GitInfoLib::GitInfo::ogs_version);
    TCLAP::ValueArg<std::string> array_name_arg(
        "n", "arrayname",
        "Name of the scalar array. Only required if assigning pixel values to "
        "cell data has been selected (default name is 'Values').",
        false, "", "name of data array");
    cmd.add(array_name_arg);
    std::vector<std::string> pixel_vals{"elevation", "materials", "scalar"};
    TCLAP::ValuesConstraint<std::string> pixel_val_options(pixel_vals);
    TCLAP::ValueArg<std::string> arg_pixel_type(
        "p", "pixel-type",
        "The choice how pixel values should be interpreted by the software: "
        "'elevation' adjusts z-coordinates; 'materials' sets (integer) "
        "material IDs; 'scalar' creates a (floating-point) array associated "
        "with mesh elements.",
        true, "", &pixel_val_options);
    cmd.add(arg_pixel_type);
    std::vector<std::string> allowed_elems{"tri", "quad"};
    TCLAP::ValuesConstraint<std::string> allowed_elem_vals(allowed_elems);
    TCLAP::ValueArg<std::string> arg_elem_type(
        "e", "elem-type", "The element type used in the resulting OGS mesh.",
        true, "", &allowed_elem_vals);
    cmd.add(arg_elem_type);
    TCLAP::ValueArg<std::string> output_arg("o", "output",
                                            "Name of the output mesh (*.vtu)",
                                            true, "", "output file name");
    cmd.add(output_arg);
    TCLAP::ValueArg<std::string> input_arg("i", "input",
                                           "Name of the input raster (*.asc)",
                                           true, "", "input file name");
    cmd.add(input_arg);
    cmd.parse(argc, argv);

    std::string const input_name = input_arg.getValue().c_str();
    std::string const output_name = output_arg.getValue().c_str();

    std::unique_ptr<GeoLib::Raster> const raster(
        FileIO::AsciiRasterInterface::getRasterFromASCFile(input_name));

    MeshLib::MeshElemType const elem_type =
        (arg_elem_type.getValue() == "tri") ? MeshLib::MeshElemType::TRIANGLE
                                            : MeshLib::MeshElemType::QUAD;
    MeshLib::UseIntensityAs intensity_type;
    if (arg_pixel_type.getValue() == "elevation")
        intensity_type = MeshLib::UseIntensityAs::ELEVATION;
    else if (arg_pixel_type.getValue() == "materials")
        intensity_type = MeshLib::UseIntensityAs::MATERIALS;
    else
        intensity_type = MeshLib::UseIntensityAs::DATAVECTOR;

    std::string array_name = "Values";
    if (intensity_type == MeshLib::UseIntensityAs::DATAVECTOR &&
        array_name_arg.isSet())
        array_name = array_name_arg.getValue().c_str();
    else if (intensity_type == MeshLib::UseIntensityAs::MATERIALS)
        array_name = "MaterialIDs";

    std::unique_ptr<MeshLib::Mesh> const mesh(MeshLib::RasterToMesh::convert(
        *raster, elem_type, intensity_type, array_name));

    if (mesh == nullptr)
    {
        ERR("Conversion failed.");
        return EXIT_FAILURE;
    }

    MeshLib::IO::VtuInterface vtu(mesh.get());
    vtu.writeToFile(output_name);
    return EXIT_SUCCESS;
}
