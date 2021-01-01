/**
 *
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <algorithm>
#include <cmath>
#include <memory>
#include <numeric>

#include <tclap/CmdLine.h>

#include "InfoLib/GitInfo.h"

#include "MeshLib/IO/readMeshFromFile.h"
#include "MeshLib/IO/writeMeshToFile.h"

#include "MeshLib/Mesh.h"

template <typename T1, typename T2>
std::pair<bool, std::string> castPropertyVectorToPropertyVector(
    MeshLib::Properties& properties,
    std::string const& property_vector_name_in,
    std::string const& property_vector_name_out)
{
    auto const* const orig_pv = properties.getPropertyVector<T1>(
        property_vector_name_in, MeshLib::MeshItemType::Cell, 1);
    if (!orig_pv)
    {
        return std::make_pair(false,
                              "Original property vector '" +
                                  property_vector_name_in + "' not found.");
    }
    auto* new_pv = properties.createNewPropertyVector<T2>(
        property_vector_name_out, MeshLib::MeshItemType::Cell, 1);
    if (!new_pv)
    {
        return std::make_pair(false,
                              "Could not create new property vector '" +
                                  property_vector_name_in + "' not found.");
    }
    new_pv->resize(orig_pv->getNumberOfTuples());
    for (std::size_t i(0); i < new_pv->getNumberOfTuples(); ++i)
    {
        (*new_pv)[i] = static_cast<T2>((*orig_pv)[i]);
    }
    return std::make_pair(true, "");
}

int main(int argc, char* argv[])
{
    TCLAP::CmdLine cmd(
        "Converts a double or floating point cell data array of a vtk "
        "unstructured grid into a int or double cell data array.\n\n"
        "OpenGeoSys-6 software, version " +
            GitInfoLib::GitInfo::ogs_version +
            ".\n"
            "Copyright (c) 2012-2021, OpenGeoSys Community "
            "(http://www.opengeosys.org)",
        ' ', GitInfoLib::GitInfo::ogs_version);

    TCLAP::ValueArg<std::string> new_property_data_type_arg(
        "t",
        "new-property-data-type",
        "the name of the data type as string (int or double)",
        false,
        "int",
        "data type as string");
    cmd.add(new_property_data_type_arg);

    TCLAP::ValueArg<std::string> new_property_arg(
        "n",
        "new-property-name",
        "the name of the new cell data array (PropertyVector) the values are "
        "stored",
        false,
        "MaterialIDs",
        "name of the new cell data array (PropertyVector) as string");
    cmd.add(new_property_arg);

    TCLAP::ValueArg<std::string> out_mesh_arg(
        "o", "out-mesh", "output mesh file name", true, "", "file name");
    cmd.add(out_mesh_arg);

    TCLAP::ValueArg<std::string> property_arg(
        "e",
        "existing-property-name",
        "the name of the existing cell data array (PropertyVector)",
        true,
        "",
        "property name as string");
    cmd.add(property_arg);

    TCLAP::ValueArg<std::string> mesh_arg(
        "i", "in-mesh", "input mesh file name", true, "", "file name");
    cmd.add(mesh_arg);

    cmd.parse(argc, argv);

    std::unique_ptr<MeshLib::Mesh> mesh(
        MeshLib::IO::readMeshFromFile(mesh_arg.getValue()));

    if (!mesh)
    {
        return -1;
    }

    bool success = false;
    std::string err_msg = "Could not find cell data array '" +
                          property_arg.getValue() + "' in the mesh '" +
                          mesh_arg.getValue() + "'";

    if (new_property_data_type_arg.getValue() == "int")
    {
        if (mesh->getProperties().existsPropertyVector<double>(
                property_arg.getValue(), MeshLib::MeshItemType::Cell, 1))
        {
            std::tie(success, err_msg) =
                castPropertyVectorToPropertyVector<double, int>(
                    mesh->getProperties(),
                    property_arg.getValue(),
                    new_property_arg.getValue());
        }

        if (mesh->getProperties().existsPropertyVector<float>(
                property_arg.getValue(), MeshLib::MeshItemType::Cell, 1))
        {
            std::tie(success, err_msg) =
                castPropertyVectorToPropertyVector<float, int>(
                    mesh->getProperties(),
                    property_arg.getValue(),
                    new_property_arg.getValue());
        }
    }
    if (new_property_data_type_arg.getValue() == "double")
    {
        if (mesh->getProperties().existsPropertyVector<float>(
                property_arg.getValue(), MeshLib::MeshItemType::Cell, 1))
        {
            std::tie(success, err_msg) =
                castPropertyVectorToPropertyVector<float, double>(
                    mesh->getProperties(),
                    property_arg.getValue(),
                    new_property_arg.getValue());
        }
    }

    if (!success)
    {
        ERR("{:s}", err_msg);
        return -1;
    }

    MeshLib::IO::writeMeshToFile(*mesh, out_mesh_arg.getValue());

    return EXIT_SUCCESS;
}
