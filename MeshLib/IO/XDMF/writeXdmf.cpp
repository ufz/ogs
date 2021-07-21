
#include "writeXdmf.h"

// ToDo (tm) Remove spdlogfmt includes with c++20
#include <spdlog/fmt/bundled/core.h>
#include <spdlog/fmt/bundled/format.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <list>
#include <map>
#include <optional>
#include <string_view>
#include <vector>

#include "BaseLib/cpp23.h"
#include "MeshLib/Location.h"

// naming convention for local function objects by postfix:
// _transform: functions that take data (mostly XDMF meta data) and return
// transformed data (XDMF string)
// _fn: functions that take a function and return a new function

using namespace fmt::literals;
using namespace BaseLib;
namespace MeshLib::IO
{
// similar definition in Location.h - XDMF uses capital letters
// Actually there is no correspondece for "IntegrationPoint" in XDMF Format
// specification
constexpr char const* mesh_item_type_strings[] = {"Node", "Edge", "Face",
                                                  "Cell", "IntegrationPoint"};

// Transforms MeshItemType into string
static std::string mesh_item_type_string(
    std::optional<MeshItemType> const& item_type)
{
    /// Char array names for all of MeshItemType values.
    if (item_type)
    {
        return mesh_item_type_strings[static_cast<int>(item_type.value())];
    }
    else
    {
        return "";
    }
}

// Transforms MeshPropertyDataType into string
// constexpr with literal type std::string not yet supported
static auto meshPropertyDatatypeString()
{
    std::array<std::string, to_underlying(MeshPropertyDataType::enum_length)>
        ogs_to_xdmf_type = {};
    ogs_to_xdmf_type[to_underlying(MeshPropertyDataType::float64)] = "Float";
    ogs_to_xdmf_type[to_underlying(MeshPropertyDataType::float32)] = "Float";
    ogs_to_xdmf_type[to_underlying(MeshPropertyDataType::int32)] = "Int";
    ogs_to_xdmf_type[to_underlying(MeshPropertyDataType::int64)] = "Int";
    ogs_to_xdmf_type[to_underlying(MeshPropertyDataType::uint32)] = "UInt";
    ogs_to_xdmf_type[to_underlying(MeshPropertyDataType::uint64)] = "UInt";
    ogs_to_xdmf_type[to_underlying(MeshPropertyDataType::int8)] = "Int";
    ogs_to_xdmf_type[to_underlying(MeshPropertyDataType::uint8)] = "UInt";
    return ogs_to_xdmf_type;
}

// Transform MeshPropertyDatatype into string )
static std::string getPropertyDataTypeString(
    MeshPropertyDataType const& ogs_data_type)
{
    auto ogs_to_xdmf_type = meshPropertyDatatypeString();
    return ogs_to_xdmf_type[to_underlying(ogs_data_type)];
}

// Returns MeshPropertyDataType-to-Size (in bytes) lookup-table
static constexpr auto meshPropertyDatatypeSize()
{
    std::array<int, to_underlying(MeshPropertyDataType::enum_length)>
        property_sizes{};
    property_sizes[to_underlying(MeshPropertyDataType::float64)] = 8;
    property_sizes[to_underlying(MeshPropertyDataType::float32)] = 4,
    property_sizes[to_underlying(MeshPropertyDataType::int32)] = 4,
    property_sizes[to_underlying(MeshPropertyDataType::int64)] = 8,
    property_sizes[to_underlying(MeshPropertyDataType::uint32)] = 4,
    property_sizes[to_underlying(MeshPropertyDataType::uint64)] = 8,
    property_sizes[to_underlying(MeshPropertyDataType::int8)] = 1,
    property_sizes[to_underlying(MeshPropertyDataType::uint8)] = 1;
    return property_sizes;
}

// Transform MeshPropertyDataType into type_sizes (in bytes)
static std::string getPropertyDataTypeSize(
    MeshPropertyDataType const& ogs_data_type)
{
    constexpr auto ogs_to_xdmf_type = meshPropertyDatatypeSize();
    return fmt::format("{}", ogs_to_xdmf_type[to_underlying(ogs_data_type)]);
}

// Collects all known data and return a function that takes all unknown data
// known: XdmfData of mesh to be written, unknown: timesteps
std::function<std::string(std::vector<double>)> write_xdmf(
    XdmfData const& geometry, XdmfData const& topology,
    std::vector<XdmfData> const& constant_attributes,
    std::vector<XdmfData> const& variable_attributes,
    std::string const& h5filename, std::string const& ogs_version)
{
    // Generates function that writes <DataItem>. Late binding needed because
    // maximum number of steps unknown. Time step and h5filename are captured to
    // return a unary function
    // _a suffix is a user defined literal for fmt named arguments
    auto const time_dataitem_genfn = [](int const time_step, int const max_step,
                                        std::string const h5filename) {
        return [time_step, max_step, h5filename](auto const& xdmfdata) {
            return fmt::format(
                "\n\t\t<DataItem DataType=\"{datatype}\" "
                "Dimensions=\"{local_dimensions}\" "
                "Format=\"HDF\" "
                "Precision=\"{precision}\">{filename}:data/{datasetname}|"
                "{time_step} {starts}:1 {strides}:1 "
                "{local_dimensions}:{max_step} "
                "{global_dimensions}</"
                "DataItem>",
                "datatype"_a = getPropertyDataTypeString(xdmfdata.data_type),
                "local_dimensions"_a =
                    fmt::join(xdmfdata.global_block_dims, " "),
                "precision"_a = getPropertyDataTypeSize(xdmfdata.data_type),
                "filename"_a = h5filename,
                "datasetname"_a = xdmfdata.name,
                "starts"_a = fmt::join(xdmfdata.starts, " "),
                "strides"_a = fmt::join(xdmfdata.strides, " "),
                "global_dimensions"_a =
                    fmt::join(xdmfdata.global_block_dims, " "),
                "time_step"_a = fmt::format("{}", time_step),
                "max_step"_a = fmt::format("{}", max_step));
        };
    };

    // string_join could be moved to ogs lib if used more
    auto const string_join_fn = [](auto const& collection) {
        return fmt::to_string(fmt::join(collection, ""));
    };

    // m_bind could be moved to ogs lib if used more
    auto const m_bind_fn = [](auto const& transform, auto const& join) {
        return [join, transform](auto const& collection) {
            std::vector<std::string> temp;
            temp.resize(collection.size());
            std::transform(collection.begin(), collection.end(),
                           std::back_inserter(temp), transform);
            return join(temp);
        };
    };

    // XDMF if part of the data that is already written in any previous step
    // subsequent steps can refer to this data collection of elements navigates
    // the xml tree (take first child, then in child take 2nd child ...)
    auto const pointer_transfrom = [](auto const& elements) {
        return fmt::format(
            "\n\t<xi:include xpointer=\"element(/{elements})\"/>",
            "elements"_a = fmt::join(elements, "/"));
    };

    // Defines content of <Attribute> in XDMF, child of Attribute is a single
    // <DataItem>
    auto const attribute_transform = [](XdmfData const& attribute,
                                        auto const& dataitem_transform) {
        return fmt::format(
            "\n\t<Attribute Center=\"{center}\" ElementCell=\"\" "
            "ElementDegree=\"0\" "
            "ElementFamily=\"\" ItemType=\"\" Name=\"{name}\" "
            "Type=\"None\">{dataitem}\n\t</Attribute>",
            "center"_a = mesh_item_type_string(attribute.attribute_center),
            "name"_a = attribute.name,
            "dataitem"_a = dataitem_transform(attribute));
    };

    // Define content of <Geometry> in XDMF, same as attribute_transform
    auto const geometry_transform = [](XdmfData const& geometry,
                                       auto const& dataitem_transform) {
        return fmt::format(
            "\n\t<Geometry Origin=\"\" Type=\"XYZ\">{dataitem}\n\t</Geometry>",
            "dataitem"_a = dataitem_transform(geometry));
    };

    // Define content of <Topology> in XDMF, same as attribute_transform
    auto const topology_transform = [](XdmfData const& topology,
                                       auto const& dataitem_transform) {
        return fmt::format(
            "\n\t<Topology Dimensions=\"{dimensions}\" "
            "Type=\"Mixed\">{dataitem}\n\t</Topology>",
            "dataitem"_a = dataitem_transform(topology),
            "dimensions"_a = fmt::join(topology.global_block_dims, " "));
    };

    // Defines content of <Grid> of a single time step, takes specific transform
    // functions for all of its child elements
    auto const grid_transform = [](double const& time_value,
                                   auto const& geometry, auto const& topology,
                                   auto const& constant_attributes,
                                   auto const& variable_attributes) {
        return fmt::format(
            "\n<Grid Name=\"Grid\">\n\t<Time "
            "Value=\"{time_value}\"/"
            ">{geometry}{topology}{fix_attributes}{variable_attributes}\n</"
            "Grid>",
            "time_value"_a = time_value, "geometry"_a = geometry,
            "topology"_a = topology, "fix_attributes"_a = constant_attributes,
            "variable_attributes"_a = variable_attributes);
    };

    // An attribute may change over time (variable) or stay constant
    enum class time_attribute
    {
        constant,
        variable
    };

    // Generates a function that either writes the data or points to existing
    // data
    auto const time_step_fn = [time_dataitem_genfn, pointer_transfrom,
                               h5filename](auto const& component_transform,
                                           int const time_step,
                                           int const max_step,
                                           time_attribute const variable) {
        return [component_transform, time_dataitem_genfn, pointer_transfrom,
                variable, time_step, max_step,
                h5filename](XdmfData const& attr) {
            // new data arrived
            bool changed =
                ((time_step > 0 && (variable == time_attribute::variable)) ||
                 time_step == 0);
            if (changed)
            {
                auto dataitem =
                    time_dataitem_genfn(time_step, max_step, h5filename);
                return component_transform(attr, dataitem);
            }
            else
            {
                std::vector<unsigned int> d = {1, 1, 2, 1, attr.index};
                return pointer_transfrom(d);
            };
        };
    };

    // Top-Level transform function (take spatial and temporal transform
    // functions) and return the time depended grid transform function
    // ToDo (tm) Remove capturing m_bind and string_join as helper function

    auto const time_grid_transform =
        [time_step_fn, m_bind_fn, string_join_fn, grid_transform,
         geometry_transform, topology_transform, attribute_transform](
            double const time, int const step, int const max_step,
            auto const& geometry, auto const& topology,
            auto const& constant_attributes, auto const& variable_attributes) {
            auto const time_step_geometry_transform = time_step_fn(
                geometry_transform, step, max_step, time_attribute::constant);
            auto const time_step_topology_transform = time_step_fn(
                topology_transform, step, max_step, time_attribute::constant);
            auto const time_step_const_attribute_fn = time_step_fn(
                attribute_transform, step, max_step, time_attribute::constant);
            auto const time_step_variable_attribute_fn = time_step_fn(
                attribute_transform, step, max_step, time_attribute::variable);

            // lift both functions from handling single attributes to work with
            // collection of attributes
            auto const variable_attributes_transform =
                m_bind_fn(time_step_variable_attribute_fn, string_join_fn);
            auto const constant_attributes_transform =
                m_bind_fn(time_step_const_attribute_fn, string_join_fn);

            return grid_transform(
                time, time_step_geometry_transform(geometry),
                time_step_topology_transform(topology),
                constant_attributes_transform(constant_attributes),
                variable_attributes_transform(variable_attributes));
        };

    // Function that combines all the data from spatial (geometry, topology,
    // attribute) and temporal domain (time) with the time domain
    // (time_step_fn). And writes <Grid CollectionType="Temporal">
    auto const temporal_grid_collection_transform =
        [time_grid_transform](
            auto const& times, auto const& geometry, auto const& topology,
            auto const& constant_attributes, auto const& variable_attributes) {
            // c++20 ranges zip missing
            auto const max_step = times.size();
            std::vector<std::string> grids;
            grids.resize(max_step);
            for (size_t time_step = 0; time_step < max_step; ++time_step)
            {
                grids.push_back(time_grid_transform(
                    times[time_step], time_step, max_step, geometry, topology,
                    constant_attributes, variable_attributes));
            }
            return fmt::format(
                "\n<Grid CollectionType=\"Temporal\" GridType=\"Collection\" "
                "Name=\"Collection\">{grids}\n</Grid>\n",
                "grids"_a = fmt::join(grids, ""));
        };

    // Generator function that return a function that represents complete xdmf
    // structure
    auto const xdmf_writer_fn = [temporal_grid_collection_transform](
                                    auto const& ogs_version,
                                    auto const& geometry, auto const& topology,
                                    auto const& constant_attributes,
                                    auto const& variable_attributes) {
        // This function will be the return of the surrounding named function
        // capture by copy - it is the function that will survive this scope!!
        // Use generalized lambda capture to move for optimization
        return [temporal_grid_collection_transform, ogs_version,
                geometry = std::move(geometry), topology = std::move(topology),
                constant_attributes = std::move(constant_attributes),
                variable_attributes = std::move(variable_attributes)](
                   std::vector<double> const& times) {
            return fmt::format(
                "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n<Xdmf "
                "xmlns:xi=\"http://www.w3.org/2001/XInclude\" "
                "Version=\"3.0\">\n<Domain>\n<Information Name=\"OGS_VERSION\" "
                "Value=\"{ogs_version}\"/>{grid_collection}</Domain>\n</Xdmf>",
                "ogs_version"_a = ogs_version,
                "grid_collection"_a = temporal_grid_collection_transform(
                    times, geometry, topology, constant_attributes,
                    variable_attributes));
        };
    };

    return xdmf_writer_fn(ogs_version, geometry, topology, constant_attributes,
                          variable_attributes);
}
}  // namespace MeshLib::IO