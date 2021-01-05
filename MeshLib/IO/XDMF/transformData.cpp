/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

// \TODO (tm) Extend xdmf lib with 64bit data types

#include "transformData.h"

#include <XdmfTopologyType.hpp>
#include <optional>
#include <string>

#include "InfoLib/GitInfo.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshEnums.h"
#include "MeshLib/Node.h"
#include "MeshPropertyDataType.h"
#include "partition.h"

namespace MeshLib::IO
{
// \TODO (tm) constexpr by other function signature can not be transformed to
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

std::optional<AttributeMeta> transformAttribute(
    std::pair<std::string, PropertyVectorBase*> const& property_pair)
{
    // 3 data that will be captured and written by lambda f below
    MeshPropertyDataType data_type = MeshPropertyDataType::unknown;
    std::size_t num_of_tuples = 0;
    void const* data_ptr = 0;

    // lambda f : Collects properties from the propertyVectorBase. It captures
    // (and overwrites) data that can only be collected via the typed property.
    // It has boolean return type to allow kind of pipe using || operator.
    auto f = [&data_type, &num_of_tuples, &data_ptr,
              &property_pair](auto basic_type) -> bool {
        auto const property_base = property_pair.second;
        auto const typed_property =
            dynamic_cast<PropertyVector<decltype(basic_type)> const*>(
                property_base);
        if (typed_property == nullptr)
        {
            return false;
        }
        // overwrite captured data
        num_of_tuples = typed_property->getNumberOfTuples();
        data_ptr = typed_property->data();

        if constexpr (std::is_same_v<double, decltype(basic_type)>)
        {
            // The standard 64-bit IEEE 754 floating-point type
            // (double-precision) has a 53 bit fractional part (52 bits written,
            // one implied)
            static_assert((std::numeric_limits<double>::digits == 53),
                          "Double has 52 bits fractional part");
            data_type = MeshPropertyDataType::float64;
        }
        else if constexpr (std::is_same_v<float, decltype(basic_type)>)
        {
            // The standard 32-bit IEEE 754 floating-point type
            // (single-precision) has a 24 bit fractional part (23 bits written,
            // one implied)
            static_assert((std::numeric_limits<float>::digits == 24),
                          "Float has 23 bits fractional part");
            data_type = MeshPropertyDataType::float32;
        }
        else if constexpr (std::is_same_v<int, decltype(basic_type)>)
        {
            static_assert((std::numeric_limits<int>::digits == 31),
                          "Signed int has 32-1 bits");
            data_type = MeshPropertyDataType::int32;
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
            data_type = MeshPropertyDataType::uint32;
        }
        // else if constexpr (std::is_same_v<unsigned long,
        // decltype(basic_type)>)
        // {
        //     static_assert((std::numeric_limits<unsigned long>::digits == 64),
        //                   "Unsigned long has 64 bits");
        //     // \TODO (tm) Extend XdmfLibrary with 64bit data types
        //     data_type = XdmfArrayType::UInt32();
        // }
        else if constexpr (std::is_same_v<std::size_t, decltype(basic_type)>)
        {
            static_assert((std::numeric_limits<std::size_t>::digits == 64),
                          "size_t has 64 bits");
            // \TODO (tm) Extend XdmfLibrary with 64bit data types
            data_type = MeshPropertyDataType::uint32;
        }
        else if constexpr (std::is_same_v<char, decltype(basic_type)>)
        {
            static_assert((std::numeric_limits<char>::digits == 7),
                          "Signed char has 8-1 bits");
            data_type = MeshPropertyDataType::int8;
        }
        else if constexpr (std::is_same_v<unsigned char, decltype(basic_type)>)
        {
            static_assert((std::numeric_limits<unsigned char>::digits == 8),
                          "Unsigned char has 8 bits");
            data_type = MeshPropertyDataType::uint8;
        }
        else
        {
            return false;
        }
        return true;
    };

    f(double{}) || f(float{}) || f(int{}) || f(long{}) || f(unsigned{}) ||
        f(long{}) || f(static_cast<unsigned long>(0)) || f(std::size_t{}) ||
        f(char{}) || f(static_cast<unsigned char>(0));

    if (data_type == MeshPropertyDataType::unknown)
    {
        return std::nullopt;
    }

    auto const& property_base = property_pair.second;
    auto const& global_components =
        property_base->getNumberOfGlobalComponents();
    // TODO (tm) property_pair vector::getNumberOfGlobalComponents should return
    // unsigned value. Then explicit cast from signed to unsigned int and
    // assert can be removed here. Implicit cast to long long is fine and
    // can be kept
    assert(global_components >= 0);
    auto const ui_global_components =
        static_cast<unsigned int>(global_components);

    MeshLib::MeshItemType const mesh_item_type =
        property_base->getMeshItemType();

    std::string const& name = property_base->getPropertyName();

    HdfData hdf =
        HdfData(data_ptr, num_of_tuples, ui_global_components, name, data_type);

    XdmfData xdmf = XdmfData(num_of_tuples, ui_global_components, data_type,
                             name, mesh_item_type);

    return AttributeMeta{std::move(hdf), std::move(xdmf)};
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
    std::string const name = "geometry";
    std::vector<MeshLib::Node*> const& nodes = mesh.getNodes();

    int const point_size = 3;
    std::vector<double> values;
    values.reserve(nodes.size() * point_size);
    for (auto const& n : nodes)
    {
        const double* x = n->getCoords();
        values.insert(values.cend(), x, x + point_size);
    }

    auto const& partition_dim = nodes.size();

    HdfData hdf = HdfData(values.data(),
                          partition_dim,
                          point_size,
                          name,
                          MeshPropertyDataType::float64);
    XdmfData xdmf = XdmfData(partition_dim,
                             point_size,
                             MeshPropertyDataType::float64,
                             name,
                             std::nullopt);

    return Geometry{std::move(values), std::move(hdf), std::move(xdmf)};
}

Topology transformTopology(MeshLib::Mesh const& mesh, std::size_t const offset)
{
    std::string const name = "topology";
    std::vector<MeshLib::Element*> const& elements = mesh.getElements();
    // \TODO (tm) Precalculate exact size
    std::vector<int> values;
    values.reserve(elements.size());

    for (auto const& cell : elements)
    {
        auto const cell_type = cellTypeOGS2XDMF(cell->getCellType());
        auto const cell_type_id = cell_type->getID();
        values.push_back(cell_type_id);
        if (cell_type_id == 2 || cell_type_id == 3)
        {
            values.push_back(cell_type->getNodesPerElement());
        }

        for (std::size_t i = 0; i < cell->getNumberOfNodes(); ++i)
        {
            MeshLib::Node const* node = cell->getNode(i);
            values.push_back(node->getID() + offset);
        }
    }

    HdfData hdf = HdfData(values.data(), values.size(), 1, name,
                          MeshPropertyDataType::int32);
    XdmfData xdmf = XdmfData(values.size(), 1, MeshPropertyDataType::int32,
                             name, std::nullopt);

    return Topology{std::move(values), std::move(hdf), std::move(xdmf)};
}
}  // namespace MeshLib::IO
