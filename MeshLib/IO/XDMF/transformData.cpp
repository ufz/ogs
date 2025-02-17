/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "transformData.h"

#include <algorithm>
#include <array>
#include <optional>
#include <range/v3/action/push_back.hpp>
#include <range/v3/view/transform.hpp>
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
    // https://xdmf.org/index.php/XDMF_Model_and_Format#Topology, Section
    // Arbitrary
    ParentDataType id;
    unsigned int number_of_nodes;
};

static constexpr auto elemOGSTypeToXDMFType()
{
    std::array<XdmfTopology, to_underlying(MeshLib::CellType::enum_length)>
        elem_type{};
    elem_type[to_underlying(MeshLib::CellType::POINT1)] = {
        ParentDataType::POLYVERTEX, 1};
    elem_type[to_underlying(MeshLib::CellType::LINE2)] = {
        ParentDataType::POLYLINE, 2};
    elem_type[to_underlying(MeshLib::CellType::LINE3)] = {
        ParentDataType::EDGE_3, 3};
    elem_type[to_underlying(MeshLib::CellType::TRI3)] = {
        ParentDataType::TRIANGLE, 3};
    elem_type[to_underlying(MeshLib::CellType::TRI6)] = {
        ParentDataType::TRIANGLE_6, 6};
    elem_type[to_underlying(MeshLib::CellType::QUAD4)] = {
        ParentDataType::QUADRILATERAL, 4};
    elem_type[to_underlying(MeshLib::CellType::QUAD8)] = {
        ParentDataType::QUADRILATERAL_8, 8};
    elem_type[to_underlying(MeshLib::CellType::QUAD9)] = {
        ParentDataType::QUADRILATERAL_9, 9};
    elem_type[to_underlying(MeshLib::CellType::TET4)] = {
        ParentDataType::TETRAHEDRON, 4};
    elem_type[to_underlying(MeshLib::CellType::TET10)] = {
        ParentDataType::TETRAHEDRON_10, 10};
    elem_type[to_underlying(MeshLib::CellType::HEX8)] = {
        ParentDataType::HEXAHEDRON, 8};
    elem_type[to_underlying(MeshLib::CellType::HEX20)] = {
        ParentDataType::HEXAHEDRON_20, 20};
    elem_type[to_underlying(MeshLib::CellType::HEX27)] = {
        ParentDataType::HEXAHEDRON_27, 27};
    elem_type[to_underlying(MeshLib::CellType::PRISM6)] = {
        ParentDataType::WEDGE, 6};  // parallel triangle wedge
    elem_type[to_underlying(MeshLib::CellType::PRISM15)] = {
        ParentDataType::WEDGE_15, 15};
    elem_type[to_underlying(MeshLib::CellType::PRISM18)] = {
        ParentDataType::WEDGE_18, 18};
    elem_type[to_underlying(MeshLib::CellType::PYRAMID5)] = {
        ParentDataType::PYRAMID, 5};
    elem_type[to_underlying(MeshLib::CellType::PYRAMID13)] = {
        ParentDataType::PYRAMID_13, 13};
    return elem_type;
}

constexpr auto elem_type_ogs2xdmf = elemOGSTypeToXDMFType();

constexpr auto cellTypeOGS2XDMF(MeshLib::CellType const& cell_type)
{
    return elem_type_ogs2xdmf[to_underlying(cell_type)];
}

std::optional<XdmfHdfData> transformAttribute(
    std::pair<std::string, PropertyVectorBase*> const& property_pair,
    unsigned int const n_files, unsigned int const chunk_size_bytes)
{
    // 3 data that will be captured and written by lambda f below
    MeshPropertyDataType data_type = MeshPropertyDataType::unknown;
    std::size_t num_of_tuples = 0;
    void const* data_ptr = 0;

    // lambda f : Collects properties from the PropertyVectorBase. It captures
    // (and overwrites) data that can only be collected via the typed property.
    // It has boolean return type to allow kind of pipe using || operator.
    auto f = [&data_type, &num_of_tuples, &data_ptr,
              &property_pair](auto basic_type) -> bool
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
        else if constexpr (std::is_same_v<int32_t, decltype(basic_type)>)
        {
            data_type = MeshPropertyDataType::int32;
        }
        else if constexpr (std::is_same_v<uint32_t, decltype(basic_type)>)
        {
            data_type = MeshPropertyDataType::uint32;
        }
        else if constexpr (std::is_same_v<int64_t, decltype(basic_type)>)
        {
            data_type = MeshPropertyDataType::int64;
        }
        else if constexpr (std::is_same_v<uint64_t, decltype(basic_type)>)
        {
            data_type = MeshPropertyDataType::uint64;
        }
        else if constexpr (std::is_same_v<int8_t, decltype(basic_type)>)
        {
            data_type = MeshPropertyDataType::int8;
        }
        else if constexpr (std::is_same_v<uint8_t, decltype(basic_type)>)
        {
            data_type = MeshPropertyDataType::uint8;
        }
        else if constexpr (std::is_same_v<char, decltype(basic_type)>)
        {
            data_type = MeshPropertyDataType::char_native;
        }
        else if constexpr (std::is_same_v<unsigned char, decltype(basic_type)>)
        {
            static_assert((std::numeric_limits<unsigned char>::digits == 8),
                          "Unsigned char has 8 bits");
            data_type = MeshPropertyDataType::uchar;
        }
        else if constexpr (std::is_same_v<unsigned long, decltype(basic_type)>)
        {
            if (sizeof(unsigned long) == 8 &&
                std::numeric_limits<unsigned long>::digits == 64)
            {
                data_type = MeshPropertyDataType::uint64;
            }
            else if (sizeof(unsigned long) == 4 &&
                     std::numeric_limits<unsigned long>::digits == 32)
            {
                data_type = MeshPropertyDataType::uint32;
            }
            else
            {
                return false;
            }
        }
        else if constexpr (std::is_same_v<std::size_t, decltype(basic_type)>)
        {
            if (sizeof(std::size_t) == 8 &&
                std::numeric_limits<std::size_t>::digits == 64)
            {
                data_type = MeshPropertyDataType::uint64;
            }
            else if (sizeof(std::size_t) == 4 &&
                     std::numeric_limits<std::size_t>::digits == 32)
            {
                data_type = MeshPropertyDataType::uint32;
            }
            else
            {
                return false;
            }
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

    HdfData hdf = {data_ptr,  num_of_tuples, ui_global_components, name,
                   data_type, n_files,       chunk_size_bytes};

    XdmfData xdmf{num_of_tuples, ui_global_components, data_type,
                  name,          mesh_item_type,       0,
                  n_files,       std::nullopt};

    return XdmfHdfData{std::move(hdf), std::move(xdmf)};
}

std::vector<XdmfHdfData> transformAttributes(
    MeshLib::Mesh const& mesh, unsigned int const n_files,
    unsigned int const chunk_size_bytes)
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

        if (!property_base->is_for_output)
        {
            continue;
        }

        if (auto const attribute = transformAttribute(

                std::pair(std::string(name), property_base), n_files,
                chunk_size_bytes))

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
        const double* x = n->data();
        values.insert(values.cend(), x, x + point_size);
    }

    return values;
}

XdmfHdfData transformGeometry(MeshLib::Mesh const& mesh, double const* data_ptr,
                              unsigned int const n_files,
                              unsigned int const chunk_size_bytes)
{
    std::string const name = "geometry";
    std::vector<MeshLib::Node*> const& nodes = mesh.getNodes();

    int const point_size = 3;
    auto const& partition_dim = nodes.size();

    HdfData const hdf = {data_ptr,
                         partition_dim,
                         point_size,
                         name,
                         MeshPropertyDataType::float64,
                         n_files,
                         chunk_size_bytes};
    XdmfData const xdmf = {
        partition_dim, point_size,   MeshPropertyDataType::float64,
        name,          std::nullopt, 2,
        n_files,       std::nullopt};

    return XdmfHdfData{std::move(hdf), std::move(xdmf)};
}

ParentDataType getTopologyType(MeshLib::Mesh const& mesh)
{
    auto const& elements = mesh.getElements();

    if (elements.empty())
    {
        return ParentDataType::MIXED;
    }
    auto const ogs_cell_type = elements[0]->getCellType();

    if (std::any_of(elements.begin(), elements.end(),
                    [&](auto const& cell)
                    { return cell->getCellType() != ogs_cell_type; }))
    {
        return ParentDataType::MIXED;
    }

    return cellTypeOGS2XDMF(ogs_cell_type).id;
}

std::pair<std::vector<std::size_t>, ParentDataType> transformToXDMFTopology(
    MeshLib::Mesh const& mesh, std::size_t const offset)
{
    std::vector<MeshLib::Element*> const& elements = mesh.getElements();
    std::vector<std::size_t> values;

    auto const push_cellnode_ids_to_vector =
        [&values, &offset](auto const& cell)
    {
        values |= ranges::actions::push_back(
            cell->nodes() | MeshLib::views::ids |
            ranges::views::transform([&offset](auto const node_id)
                                     { return node_id + offset; }));
    };

    auto const topology_type = getTopologyType(mesh);
    if (topology_type == ParentDataType::MIXED)
    {
        values.reserve(elements.size() * 2);  // each cell has at least two
                                              // numbers
        for (auto const& cell : elements)
        {
            auto const ogs_cell_type = cell->getCellType();
            auto const xdmf_cell_id = cellTypeOGS2XDMF(ogs_cell_type).id;
            values.push_back(to_underlying(xdmf_cell_id));
            push_cellnode_ids_to_vector(cell);
        }
    }
    else
    {
        values.reserve(elements.size() * elements[0]->getNumberOfNodes());
        for (auto const& cell : elements)
        {
            push_cellnode_ids_to_vector(cell);
        }
    }

    return {values, topology_type};
}

XdmfHdfData transformTopology(std::vector<std::size_t> const& values,
                              ParentDataType const parent_data_type,
                              unsigned int const n_files,
                              unsigned int const chunk_size_bytes)
{
    std::string const name = "topology";
    auto const tuple_size = ParentDataType2String(parent_data_type).second;
    HdfData const hdf = {values.data(),
                         values.size() / tuple_size,
                         tuple_size,
                         name,
                         MeshPropertyDataType::uint64,
                         n_files,
                         chunk_size_bytes};
    XdmfData const xdmf{values.size() / tuple_size,
                        tuple_size,
                        MeshPropertyDataType::uint64,
                        name,
                        std::nullopt,
                        3,
                        n_files,
                        parent_data_type};

    return XdmfHdfData{std::move(hdf), std::move(xdmf)};
}
}  // namespace MeshLib::IO
