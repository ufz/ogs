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

#include <algorithm>
#include <array>
#include <optional>
#include <string>

#include "BaseLib/cpp23.h"
#include "InfoLib/GitInfo.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshEnums.h"
#include "MeshLib/Node.h"
#include "MeshPropertyDataType.h"
#include "partition.h"

using namespace BaseLib;
namespace MeshLib::IO
{
struct XdmfTopology
{
    // https://www.xdmf.org/index.php/XDMF_Model_and_Format#Topology, Section
    // Arbitrary
    unsigned int id;
    unsigned int number_of_nodes;
};

static constexpr auto elemOGSTypeToXDMFType()
{
    std::array<XdmfTopology, to_underlying(MeshLib::CellType::enum_length)>
        elem_type{};
    elem_type[to_underlying(MeshLib::CellType::POINT1)] = {0x1, 1};
    elem_type[to_underlying(MeshLib::CellType::LINE2)] = {0x2, 2};
    elem_type[to_underlying(MeshLib::CellType::LINE3)] = {0x2, 3};
    elem_type[to_underlying(MeshLib::CellType::TRI3)] = {0x4, 3};
    elem_type[to_underlying(MeshLib::CellType::TRI6)] = {0x24, 6};
    elem_type[to_underlying(MeshLib::CellType::QUAD4)] = {0x5, 4};
    elem_type[to_underlying(MeshLib::CellType::QUAD8)] = {0x25, 8};
    elem_type[to_underlying(MeshLib::CellType::QUAD9)] = {0x23, 9};
    elem_type[to_underlying(MeshLib::CellType::TET4)] = {0x6, 4};
    elem_type[to_underlying(MeshLib::CellType::TET10)] = {0x26, 10};
    elem_type[to_underlying(MeshLib::CellType::HEX8)] = {0x9, 8};
    elem_type[to_underlying(MeshLib::CellType::HEX20)] = {0x30, 20};
    elem_type[to_underlying(MeshLib::CellType::HEX27)] = {0x32, 27};
    elem_type[to_underlying(MeshLib::CellType::PRISM6)] = {
        0x8, 6};  // parallel triangle wedge
    elem_type[to_underlying(MeshLib::CellType::PRISM15)] = {0x28, 15};
    elem_type[to_underlying(MeshLib::CellType::PRISM18)] = {0x29, 18};
    elem_type[to_underlying(MeshLib::CellType::PYRAMID5)] = {0x7, 5};
    elem_type[to_underlying(MeshLib::CellType::PYRAMID13)] = {0x27, 13};
    return elem_type;
}

constexpr auto elem_type_ogs2xdmf = elemOGSTypeToXDMFType();

constexpr auto cellTypeOGS2XDMF(MeshLib::CellType const& cell_type)
{
    return elem_type_ogs2xdmf[to_underlying(cell_type)];
}

std::optional<XdmfHdfData> transformAttribute(
    std::pair<std::string, PropertyVectorBase*> const& property_pair)
{
    // 3 data that will be captured and written by lambda f below
    MeshPropertyDataType data_type = MeshPropertyDataType::unknown;
    std::size_t num_of_tuples = 0;
    void const* data_ptr = 0;

    // lambda f : Collects properties from the propertyVectorBase. It captures
    // (and overwrites) data that can only be collected via the typed property.
    // It has boolean return type to allow kind of pipe using || operator.
    auto f = [&data_type, &num_of_tuples, &data_ptr, &property_pair](
                 auto basic_type) -> bool
    {
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
        // ToDo (tm) These tests are platform specific and currently fail on
        // Windows else if constexpr (std::is_same_v<long,
        // decltype(basic_type)>)
        //{
        //    static_assert((std::numeric_limits<long>::digits == 63),
        //                  "Signed int has 64-1 bits");
        //    data_type = MeshPropertyDataType::int64;
        //}
        // else if constexpr (std::is_same_v<unsigned long,
        // decltype(basic_type)>)
        //{
        //    static_assert((std::numeric_limits<unsigned long>::digits == 64),
        //                  "Unsigned long has 64 bits");
        //    data_type = MeshPropertyDataType::uint64;
        //}
        else if constexpr (std::is_same_v<unsigned int, decltype(basic_type)>)
        {
            static_assert((std::numeric_limits<unsigned int>::digits == 32),
                          "Unsigned int has 32 bits");
            data_type = MeshPropertyDataType::uint32;
        }
        else if constexpr (std::is_same_v<std::size_t, decltype(basic_type)>)
        {
            static_assert((std::numeric_limits<std::size_t>::digits == 64),
                          "size_t has 64 bits");
            data_type = MeshPropertyDataType::uint64;
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

    HdfData hdf = {data_ptr, num_of_tuples, ui_global_components, name,
                   data_type};

    XdmfData xdmf = {num_of_tuples, ui_global_components, data_type,
                     name,          mesh_item_type,       0};

    return XdmfHdfData{std::move(hdf), std::move(xdmf)};
}

std::vector<XdmfHdfData> transformAttributes(MeshLib::Mesh const& mesh)
{
    MeshLib::Properties const& properties = mesh.getProperties();

    // \TODO (tm) use c++20 ranges
    // a = p | filter (first!=OGS_VERSION) | filter null_opt | transformAttr |
    std::vector<XdmfHdfData> attributes;
    for (auto const& [name, property_base] : properties)
    {
        if (name == GitInfoLib::GitInfo::OGS_VERSION)
        {
            continue;
        }

        if (auto const attribute =
                transformAttribute(std::pair(name, property_base)))
        {
            attributes.push_back(attribute.value());
        }
        else
        {
            WARN("Could not create attribute meta of {:s}.", name);
        }
    }
    return attributes;
}

std::vector<double> transformToXDMFGeometry(MeshLib::Mesh const& mesh)
{
    std::vector<MeshLib::Node*> const& nodes = mesh.getNodes();

    int const point_size = 3;
    std::vector<double> values;
    values.reserve(nodes.size() * point_size);
    for (auto const& n : nodes)
    {
        const double* x = n->getCoords();
        values.insert(values.cend(), x, x + point_size);
    }

    return values;
}

XdmfHdfData transformGeometry(MeshLib::Mesh const& mesh, double const* data_ptr)
{
    std::string const name = "geometry";
    std::vector<MeshLib::Node*> const& nodes = mesh.getNodes();

    int const point_size = 3;
    auto const& partition_dim = nodes.size();

    HdfData const hdf = {data_ptr, partition_dim, point_size, name,
                         MeshPropertyDataType::float64};
    XdmfData const xdmf = {
        partition_dim, point_size,   MeshPropertyDataType::float64,
        name,          std::nullopt, 2};

    return XdmfHdfData{std::move(hdf), std::move(xdmf)};
}

std::vector<int> transformToXDMFTopology(MeshLib::Mesh const& mesh,
                                         std::size_t const offset)
{
    std::vector<MeshLib::Element*> const& elements = mesh.getElements();
    std::vector<int> values;
    values.reserve(elements.size());

    for (auto const& cell : elements)
    {
        auto const ogs_cell_type = cell->getCellType();
        auto const xdmf_cell_id = cellTypeOGS2XDMF(ogs_cell_type).id;
        values.push_back(xdmf_cell_id);
        if (ogs_cell_type == MeshLib::CellType::LINE2 ||
            ogs_cell_type == MeshLib::CellType::LINE3)
        {
            values.push_back(cellTypeOGS2XDMF(ogs_cell_type).number_of_nodes);
        }

        for (std::size_t i = 0; i < cell->getNumberOfNodes(); ++i)
        {
            MeshLib::Node const* node = cell->getNode(i);
            values.push_back(node->getID() + offset);
        }
    }
    return values;
}

XdmfHdfData transformTopology(std::vector<int> const& values)
{
    std::string const name = "topology";
    HdfData const hdf = {values.data(), values.size(), 1, name,
                         MeshPropertyDataType::int32};
    XdmfData const xdmf = {values.size(), 1, MeshPropertyDataType::int32, name,
                           std::nullopt,  3};

    return XdmfHdfData{std::move(hdf), std::move(xdmf)};
}
}  // namespace MeshLib::IO
