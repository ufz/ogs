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
void searchByPropertyValue(std::string const& property_name,
                           std::vector<PROPERTY_TYPE> const& property_values,
                           MeshLib::ElementSearch& searcher)
{
    for (auto const& property_value : property_values)
    {
        std::size_t n_marked_elements = searcher.searchByPropertyValue<double>(
            property_name, property_value);
        if (n_marked_elements == 0)
            n_marked_elements = searcher.searchByPropertyValue<int>(
                property_name, property_value);

        INFO("%d elements with property value %s found.", n_marked_elements,
             std::to_string(property_value).c_str());
    }
}

void searchByPropertyRange(std::string const& property_name,
                           double const& min_value, double const& max_value,
                           bool const& outside,
                           MeshLib::ElementSearch& searcher)
{
    std::size_t n_marked_elements = searcher.searchByPropertyValueRange<double>(
        property_name, min_value, max_value, outside);

    if (n_marked_elements == 0)
        n_marked_elements = searcher.searchByPropertyValueRange<int>(
            property_name, static_cast<int>(min_value),
            static_cast<int>(max_value), outside);

    // add checks for other data types here (if n_marked_elements remains 0)

    INFO("%d elements in range [%s, %s] found.", n_marked_elements,
         std::to_string(min_value).c_str(), std::to_string(max_value).c_str());
}

int main (int argc, char* argv[])
{
    ApplicationsLib::LogogSetup logog_setup;

    TCLAP::CmdLine cmd("Removes mesh elements based on element type, element volume, scalar "
        "arrays, or bounding box . The documentation is available at "
        "https://docs.opengeosys.org/docs/tools/meshing/remove-mesh-elements.", ' ', "0.1");

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

    // scalar array params
    TCLAP::ValueArg<std::string> property_name_arg(
        "n", "property-name", "name of property in the mesh", false, "MaterialIDs", "string");
    cmd.add(property_name_arg);

    TCLAP::MultiArg<int> property_arg(
        "", "property-value", "value of selected property to be removed", false, "number");
    cmd.add(property_arg);

    TCLAP::ValueArg<double> min_property_arg(
        "", "min-value", "minimum value of range for selected property", false, 0, "number");
    cmd.add(min_property_arg);

    TCLAP::ValueArg<double> max_property_arg(
        "", "max-value", "maximum value of range for selected property", false, 0, "number");
    cmd.add(max_property_arg);

    TCLAP::SwitchArg outside_property_arg(
        "", "outside", "remove all elements outside the given property range");
    cmd.add(outside_property_arg);

    TCLAP::SwitchArg inside_property_arg(
        "", "inside", "remove all elements inside the given property range");
    cmd.add(inside_property_arg);

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
    if (mesh == nullptr)
        return EXIT_FAILURE;

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

    if (property_name_arg.isSet() || property_arg.isSet() ||
        min_property_arg.isSet() || max_property_arg.isSet())
    {
        if ((property_arg.isSet() || min_property_arg.isSet() || max_property_arg.isSet()) &&
            !property_name_arg.isSet())
        {
            ERR("Specify a property name for the value/range selected.");
            return EXIT_FAILURE;
        }

        if (property_name_arg.isSet() &&
            !((min_property_arg.isSet() && max_property_arg.isSet()) || property_arg.isSet()))
        {
            ERR("Specify a value or range (\"-min-value\" and \"-max_value\") "
                "for the property selected.");
            return EXIT_FAILURE;
        }

        // name + value
        if (property_arg.isSet() && property_name_arg.isSet())
        {
            searchByPropertyValue(property_name_arg.getValue(),
                                  property_arg.getValue(), searcher);
        }

        // name + range
        if (property_name_arg.isSet() && min_property_arg.isSet() &&
            max_property_arg.isSet())
        {
            if ((!outside_property_arg.isSet() &&
                 !inside_property_arg.isSet()) ||
                (outside_property_arg.isSet() && inside_property_arg.isSet()))
            {
                ERR("Specify if the inside or the outside of the selected "
                    "range should be removed.");
                return EXIT_FAILURE;
            }

            bool const outside = outside_property_arg.isSet();
            searchByPropertyRange(
                property_name_arg.getValue(), min_property_arg.getValue(),
                max_property_arg.getValue(), outside, searcher);
        }
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
             searcher.searchByBoundingBox(GeoLib::AABB(extent.begin(), extent.end())));
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
