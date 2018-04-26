/**
 * @file   removeMeshElements.cpp
 * @author Norihiro Watanabe
 * @date   2013/10/15
 * @brief  Remove mesh elements
 *
 * @copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include <memory>

#include <tclap/CmdLine.h>

#include "Applications/ApplicationsLib/LogogSetup.h"

#include "MeshLib/IO/readMeshFromFile.h"
#include "MeshLib/IO/writeMeshToFile.h"

#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/MeshEnums.h"
#include "MeshLib/MeshSearch/ElementSearch.h"
#include "MeshLib/MeshEditing/RemoveMeshComponents.h"

template <typename PROPERTY_TYPE>
void searchByProperty(std::string const& property_name,
                      std::vector<PROPERTY_TYPE> const& property_values,
                      MeshLib::ElementSearch& searcher)
{
    for (auto const& property_value : property_values) {
        const std::size_t n_marked_elements =
            searcher.searchByPropertyValue(property_name, property_value);
        INFO("%d elements with property value %s found.", n_marked_elements,
             std::to_string(property_value).c_str());
    }
}

int main (int argc, char* argv[])
{
    ApplicationsLib::LogogSetup logog_setup;

    TCLAP::CmdLine cmd(
        "Remove mesh elements. The documentation is available at "
        "https://docs.opengeosys.org/docs/tools/meshing/remove-mesh-elements.",
        ' ', "0.1");

    // Bounding box params
    TCLAP::ValueArg<double> zLargeArg("", "z-max", "largest allowed extent in z-dimension",
                                      false, std::numeric_limits<double>::max(), "value");
    cmd.add(zLargeArg);
    TCLAP::ValueArg<double> zSmallArg("", "z-min", "smallest allowed extent in z-dimension",
                                      false,  -1 * std::numeric_limits<double>::max(), "value");
    cmd.add(zSmallArg);
    TCLAP::ValueArg<double> yLargeArg("", "y-max", "largest allowed extent in y-dimension",
                                      false, std::numeric_limits<double>::max(), "value");
    cmd.add(yLargeArg);
    TCLAP::ValueArg<double> ySmallArg("", "y-min", "smallest allowed extent in y-dimension",
                                       false,  -1 * std::numeric_limits<double>::max(), "value");
    cmd.add(ySmallArg);
    TCLAP::ValueArg<double> xLargeArg("", "x-max", "largest allowed extent in x-dimension",
                                       false, std::numeric_limits<double>::max(), "value");
    cmd.add(xLargeArg);
    TCLAP::ValueArg<double> xSmallArg("", "x-min", "smallest allowed extent in x-dimension",
                                      false, -1 * std::numeric_limits<double>::max(), "value");
    cmd.add(xSmallArg);

    // Non-bounding-box params
    TCLAP::SwitchArg zveArg("z", "zero-volume", "remove zero volume elements", false);
    cmd.add(zveArg);
    TCLAP::MultiArg<std::string> eleTypeArg("t", "element-type",
                                          "element type to be removed", false, "element type");
    cmd.add(eleTypeArg);

    TCLAP::MultiArg<int> int_property_arg("", "int-property-value",
                                          "new property value (data type int)",
                                          false, "number");
    cmd.add(int_property_arg);
    TCLAP::ValueArg<std::string> property_name_arg(
        "n", "property-name", "name of property in the mesh", false,
        "MaterialIDs", "string");
    cmd.add(property_name_arg);

    // I/O params
    TCLAP::ValueArg<std::string> mesh_out("o", "mesh-output-file",
                                          "the name of the file the mesh will be written to", true,
                                          "", "file name of output mesh");
    cmd.add(mesh_out);
    TCLAP::ValueArg<std::string> mesh_in("i", "mesh-input-file",
                                         "the name of the file containing the input mesh", true,
                                         "", "file name of input mesh");
    cmd.add(mesh_in);

    cmd.parse(argc, argv);

    std::unique_ptr<MeshLib::Mesh const> mesh(
        MeshLib::IO::readMeshFromFile(mesh_in.getValue()));
    INFO("Mesh read: %d nodes, %d elements.", mesh->getNumberOfNodes(), mesh->getNumberOfElements());
    MeshLib::ElementSearch searcher(*mesh);

    // search elements IDs to be removed
    if (zveArg.isSet()) {
        INFO("%d zero volume elements found.", searcher.searchByContent());
    }
    if (eleTypeArg.isSet()) {
        const std::vector<std::string> eleTypeNames = eleTypeArg.getValue();
        for (const auto& typeName : eleTypeNames) {
            const MeshLib::MeshElemType type = MeshLib::String2MeshElemType(typeName);
            if (type == MeshLib::MeshElemType::INVALID) continue;
            INFO("%d %s elements found.", searcher.searchByElementType(type),
                 typeName.c_str());
        }
    }

    if (int_property_arg.isSet()) {
        searchByProperty(property_name_arg.getValue(),
                         int_property_arg.getValue(), searcher);
    }

    if (xSmallArg.isSet() || xLargeArg.isSet() ||
        ySmallArg.isSet() || yLargeArg.isSet() ||
        zSmallArg.isSet() || zLargeArg.isSet())
    {
        bool aabb_error (false);
        if (xSmallArg.getValue() >= xLargeArg.getValue())
        {
            ERR ("Minimum x-extent larger than maximum x-extent.");
            aabb_error = true;
        }
        if (ySmallArg.getValue() >= yLargeArg.getValue())
        {
            ERR ("Minimum y-extent larger than maximum y-extent.");
            aabb_error = true;
        }
        if (zSmallArg.getValue() >= zLargeArg.getValue())
        {
            ERR ("Minimum z-extent larger than maximum z-extent.");
            aabb_error = true;
        }
        if (aabb_error)
            return EXIT_FAILURE;

        std::array<MathLib::Point3d, 2> extent({{
            MathLib::Point3d(std::array<double,3>{{xSmallArg.getValue(),
                ySmallArg.getValue(), zSmallArg.getValue()}}),
            MathLib::Point3d(std::array<double,3>{{xLargeArg.getValue(),
                yLargeArg.getValue(), zLargeArg.getValue()}})}});
        INFO("%d elements found.",
             searcher.searchByBoundingBox(
                 GeoLib::AABB(extent.begin(), extent.end())));
    }

    // remove the elements and create a new mesh object.
    std::unique_ptr<MeshLib::Mesh const> new_mesh(MeshLib::removeElements(
        *mesh, searcher.getSearchedElementIDs(), mesh->getName()));

    if (new_mesh == nullptr)
        return EXIT_FAILURE;

    // write into a file
    MeshLib::IO::writeMeshToFile(*new_mesh, mesh_out.getValue());

    return EXIT_SUCCESS;
}



