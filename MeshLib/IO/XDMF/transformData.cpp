/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

// \TODO (tm) Extend xdmf lib with 64bit datatypes
#if !_MSC_VER
#pragma GCC diagnostic ignored "-Wnarrowing"
#endif

#include "transformData.h"

#include <XdmfArrayType.hpp>
#include <XdmfAttributeCenter.hpp>
#include <XdmfTopologyType.hpp>
#include <optional>
#include <variant>

#include "InfoLib/GitInfo.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshEnums.h"
#include "MeshLib/Node.h"

namespace
{
using XdmfDimType = unsigned int;
using Hdf5DimType = unsigned long long;
}  // namespace

namespace MeshLib::IO
{
// \TODO (tm) constexpr by other funtion signature can not be transformed to
// constexpr shared_ptr is not literal type / has non trivial destructor
boost::shared_ptr<const XdmfTopologyType> cellTypeOGS2XDMF(
    MeshLib::CellType const cell_type)
{
    static std::map<MeshLib::CellType const,
                    boost::shared_ptr<const XdmfTopologyType> const>
        elem_type_ogs2xdmf = {
            {MeshLib::CellType::POINT1, XdmfTopologyType::Polyvertex()},
            {MeshLib::CellType::LINE2, XdmfTopologyType::Polyline(2)},
            {MeshLib::CellType::LINE3, XdmfTopologyType::Polyline(3)},
            {MeshLib::CellType::QUAD4, XdmfTopologyType::Quadrilateral()},
            {MeshLib::CellType::QUAD8, XdmfTopologyType::Quadrilateral_8()},
            {MeshLib::CellType::QUAD9, XdmfTopologyType::Quadrilateral_9()},
            {MeshLib::CellType::TET4, XdmfTopologyType::Tetrahedron()},
            {MeshLib::CellType::TET10, XdmfTopologyType::Tetrahedron_10()},
            {MeshLib::CellType::TRI3, XdmfTopologyType::Triangle()},
            {MeshLib::CellType::TRI6, XdmfTopologyType::Triangle_6()},
            {MeshLib::CellType::PRISM6,
             XdmfTopologyType::Wedge()},  // parallel triangle wedge
            {MeshLib::CellType::PRISM15, XdmfTopologyType::Wedge_15()},
            {MeshLib::CellType::PRISM18, XdmfTopologyType::Wedge_18()},
            {MeshLib::CellType::PYRAMID13, XdmfTopologyType::Pyramid_13()},
            {MeshLib::CellType::PYRAMID5, XdmfTopologyType::Pyramid()},
            {MeshLib::CellType::HEX20, XdmfTopologyType::Hexahedron_20()},
            {MeshLib::CellType::HEX27, XdmfTopologyType::Hexahedron_27()},
            {MeshLib::CellType::HEX8, XdmfTopologyType::Hexahedron()}};

    return elem_type_ogs2xdmf.at(cell_type);
}

boost::shared_ptr<const XdmfAttributeCenter> elemTypeOGS2XDMF(
    MeshLib::MeshItemType const elem_type)
{
    std::map<MeshLib::MeshItemType,
             boost::shared_ptr<const XdmfAttributeCenter>>
        mesh_item_type_ogs2xdmf = {
            {MeshLib::MeshItemType::Cell, XdmfAttributeCenter::Cell()},
            {MeshLib::MeshItemType::Edge, XdmfAttributeCenter::Edge()},
            {MeshLib::MeshItemType::Face, XdmfAttributeCenter::Face()},
            {MeshLib::MeshItemType::Node, XdmfAttributeCenter::Node()},
            {MeshLib::MeshItemType::IntegrationPoint,
             XdmfAttributeCenter::Other()}};

    return mesh_item_type_ogs2xdmf.at(elem_type);
}

std::optional<AttributeMeta> transformAttribute(
    std::pair<std::string, PropertyVectorBase*> const& property)
{
    // lambda f : Transforms property to AttributeMeta, first parameter is
    // output parameter because the boolean return type is used to allow kind of
    // pipe using || operator
    auto f = [&](std::optional<AttributeMeta>& attribute_meta,
                 auto basic_type,
                 auto property_pair) -> bool {
        auto const property_base = property_pair.second;
        auto const name = property_pair.first;
        auto const typed_property =
            dynamic_cast<PropertyVector<decltype(basic_type)> const*>(
                property_base);
        if (typed_property == nullptr)
        {
            attribute_meta = std::nullopt;
            return false;
        }

        auto const mesh_item_type =
            elemTypeOGS2XDMF(property_base->getMeshItemType());
        auto const global_components =
            property_base->getNumberOfGlobalComponents();
        auto const size = typed_property->getNumberOfTuples();
        auto const vdims = [](XdmfDimType num_components,
                              XdmfDimType size) -> std::vector<XdmfDimType> {
            if (num_components > 1)
            {
                return {size, num_components};
            }
            return {size};
        }(global_components, size);

        // \TODO (tm) Remove code duplicationby eliminating the need for a
        // second vldim at all by modification of XdmfHdf5Controller
        // taking unsigned long long
        auto const vldims = [](XdmfDimType num_components,
                               XdmfDimType size) -> std::vector<Hdf5DimType> {
            if (num_components > 1)
            {
                return std::vector<Hdf5DimType>{size, num_components};
            }
            else
            {
                return std::vector<Hdf5DimType>{size};
            }
        }(global_components, size);

        auto data_type = XdmfArrayType::Float64();

        if constexpr (std::is_same_v<double, decltype(basic_type)>)
        {
            // The standard 64-bit IEEE 754 floating-point type (double-precision)
            // has a 53 bit fractional part (52 bits written, one implied)
            static_assert((std::numeric_limits<double>::digits == 53),
                          "Double has 52 bits fractional part");
            data_type = XdmfArrayType::Float64();
        }
        else if constexpr (std::is_same_v<float, decltype(basic_type)>)
        {
            // The standard 32-bit IEEE 754 floating-point type (single-precision)
            // has a 24 bit fractional part (23 bits written, one implied)
            static_assert((std::numeric_limits<float>::digits == 24),
                          "Float has 23 bits fractional part");
            data_type = XdmfArrayType::Float32();
        }
        else if constexpr (std::is_same_v<int, decltype(basic_type)>)
        {
            static_assert((std::numeric_limits<int>::digits == 31),
                          "Signed int has 32-1 bits");
            data_type = XdmfArrayType::Int32();
        }
        // \TODO (tm) Reimplement size checks
        // else if constexpr (std::is_same_v<long, decltype(basic_type)>)
        // {
        //     static_assert((std::numeric_limits<long>::digits == 63),
        //                   "Signed int has 64-1 bits");
        //     data_type = XdmfArrayType::Int64();
        // }
        else if constexpr (std::is_same_v<unsigned int, decltype(basic_type)>)
        {
            static_assert((std::numeric_limits<unsigned int>::digits == 32),
                          "Unsigned int has 32 bits");
            data_type = XdmfArrayType::UInt32();
        }
        // else if constexpr (std::is_same_v<unsigned long, decltype(basic_type)>)
        // {
        //     static_assert((std::numeric_limits<unsigned long>::digits == 64),
        //                   "Unsigned long has 64 bits");
        //     // \TODO (tm) Extend XdmfLibrary with 64bit datatypes
        //     data_type = XdmfArrayType::UInt32();
        // }
        else if constexpr (std::is_same_v<std::size_t, decltype(basic_type)>)
        {
            static_assert((std::numeric_limits<std::size_t>::digits == 64),
                          "size_t has 64 bits");
            // \TODO (tm) Extend XdmfLibrary with 64bit datatypes
            data_type = XdmfArrayType::UInt32();
        }
        else if constexpr (std::is_same_v<char, decltype(basic_type)>)
        {
            static_assert((std::numeric_limits<char>::digits == 7),
                          "char has 8-1 bits");
            data_type = XdmfArrayType::Int8();
        }
        else if constexpr (std::is_same_v<unsigned char, decltype(basic_type)>)
        {
            static_assert((std::numeric_limits<unsigned char>::digits == 8),
                          "Unsigned char has 8 bits");
            data_type = XdmfArrayType::UInt8();
        }
        else
        {
            return false;
        }

        std::vector<XdmfDimType> const starts = {0, 0, 0};
        std::vector<XdmfDimType> const strides = {1, 1, 1};

        attribute_meta = {typed_property->data(),
                          name,
                          mesh_item_type,
                          data_type,
                          starts,
                          strides,
                          vdims,
                          vldims};

        return true;
    };

    std::optional<AttributeMeta> attribute;
    f(attribute, double{}, property) || f(attribute, float{}, property) ||
        f(attribute, int{}, property) || f(attribute, long{}, property) ||
        f(attribute, unsigned{}, property) || f(attribute, long{}, property) ||
        f(attribute, static_cast<unsigned long>(0), property) ||
        f(attribute, std::size_t{}, property) ||
        f(attribute, char{}, property) ||
        f(attribute, static_cast<unsigned char>(0), property);

    if (!attribute)
    {
        OGS_FATAL("Could not apply function to PropertyVector '{:s}'.",
                  property.second->getPropertyName());
        return std::nullopt;
    }

    return attribute;
}

std::vector<AttributeMeta> transformAttributes(MeshLib::Mesh const& mesh)
{
    MeshLib::Properties const& properties = mesh.getProperties();

    // \TODO (tm) use c++20 ranges
    // a = p | filter (first!=OGS_VERSION) | filter null_opt | transformAttr |
    // optional_value
    std::vector<AttributeMeta> attributes;
    for (auto [name, property_base] : properties)
    {
        if (name == GitInfoLib::GitInfo::OGS_VERSION)
        {
            continue;
        }

        auto attribute = transformAttribute(std::pair(name, property_base));
        if (attribute)
        {
            attributes.push_back(attribute.value());
        }
    }
    return attributes;
}

Geometry transformGeometry(MeshLib::Mesh const& mesh)
{
    std::vector<MeshLib::Node*> const& nodes = mesh.getNodes();

    std::vector<double> values;
    values.reserve(nodes.size() * 3);
    for (auto n : nodes)
    {
        const double* x = n->getCoords();
        values.insert(values.cend(), x, x + 3);
    }

    std::vector<XdmfDimType> const vdims = {static_cast<XdmfDimType>(nodes.size()), 3};
    std::vector<Hdf5DimType> const vldims = {nodes.size(), 3};
    std::vector<XdmfDimType> const starts = {0, 0, 0};
    std::vector<XdmfDimType> const strides = {1, 1, 1};
    return Geometry{std::move(values), starts, strides, vdims, vldims};
}

Topology transformTopology(MeshLib::Mesh const& mesh)
{
    std::vector<MeshLib::Element*> elements = mesh.getElements();
    // \TODO (tm) Precalculate exact size
    std::vector<int> values;
    values.reserve(elements.size());

    for (auto const& cell : elements)
    {
        auto const celltype = cellTypeOGS2XDMF(cell->getCellType());
        auto const celltype_id = celltype->getID();
        values.push_back(celltype_id);
        if (celltype_id == 2 || celltype_id == 3)
        {
            values.push_back(celltype->getNodesPerElement());
        }

        for (std::size_t i = 0; i < cell->getNumberOfNodes(); ++i)
        {
            MeshLib::Node const* node = cell->getNode(i);
            values.push_back(node->getID());
        }
    }

    std::vector<XdmfDimType> const starts = {0};
    std::vector<XdmfDimType> const strides = {1};
    std::vector<XdmfDimType> const vdims = {static_cast<XdmfDimType>(values.size())};
    std::vector<Hdf5DimType> const vldims = {values.size()};

    return Topology{std::move(values), starts, strides, vdims, vldims};
}
}  // namespace MeshLib::IO
