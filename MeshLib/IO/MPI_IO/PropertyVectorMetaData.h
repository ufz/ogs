/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <optional>
#include <string>

namespace MeshLib
{
namespace IO
{

struct PropertyVectorMetaData
{
    std::string property_name;
    /// is_int_type is true if the type of the components is an integer type, if
    /// it is a floating point number type the is_int_type is false
    bool is_int_type;
    /// if the component type is an integer number the flag is_data_type_signed
    /// signals if it has a sign or not
    bool is_data_type_signed;
    unsigned long data_type_size_in_bytes;
    unsigned long number_of_components;
    unsigned long number_of_tuples;

    template <typename T>
    void fillPropertyVectorMetaDataTypeInfo()
    {
        is_int_type = std::numeric_limits<T>::is_integer;
        is_data_type_signed = std::numeric_limits<T>::is_signed;
        data_type_size_in_bytes = sizeof(T);
    }
};

inline void writePropertyVectorMetaData(
    std::ostream& os, PropertyVectorMetaData const& pvmd)
{
    std::string::size_type s(pvmd.property_name.length());
    os.write(reinterpret_cast<char*>(&s), sizeof(std::string::size_type));

    os.write(
        const_cast<char*>(
            const_cast<PropertyVectorMetaData&>(pvmd).property_name.data()),
        s);
    os.write(reinterpret_cast<char*>(
                 &const_cast<PropertyVectorMetaData&>(pvmd).is_int_type),
             sizeof(bool));
    os.write(reinterpret_cast<char*>(&const_cast<PropertyVectorMetaData&>(
                 pvmd).is_data_type_signed),
             sizeof(bool));
    os.write(reinterpret_cast<char*>(&const_cast<PropertyVectorMetaData&>(
                 pvmd).data_type_size_in_bytes),
             sizeof(unsigned long));
    os.write(reinterpret_cast<char*>(&const_cast<PropertyVectorMetaData&>(
                 pvmd).number_of_components),
             sizeof(unsigned long));
    os.write(reinterpret_cast<char*>(
                 &const_cast<PropertyVectorMetaData&>(pvmd).number_of_tuples),
             sizeof(unsigned long));
}

inline void writePropertyVectorMetaData(PropertyVectorMetaData const& pvmd)
{
    DBUG(
        "name: '{:s}':\t is_int_data_type: {:d}, is_data_type_signed: "
        "{:d}, data_type_size_in_bytes: {:d}, number of components / "
        "tuples: {:d} / {:d}",
        pvmd.property_name, pvmd.is_int_type, pvmd.is_data_type_signed,
        pvmd.data_type_size_in_bytes, pvmd.number_of_components,
        pvmd.number_of_tuples);
}

inline std::optional<PropertyVectorMetaData> readPropertyVectorMetaData(
    std::istream& is)
{
    // read the size of the name of the PropertyVector
    std::string::size_type s = 0;
    if (!is.read(reinterpret_cast<char*>(&s), sizeof(std::string::size_type)))
    {
        return std::optional<PropertyVectorMetaData>();
    }

    PropertyVectorMetaData pvmd;
    char *dummy = new char[s];
    if (!is.read(dummy, s))
    {
        return std::nullopt;
    }
    pvmd.property_name = std::string(dummy, s);
    delete [] dummy;

    if(!is.read(reinterpret_cast<char*>(&pvmd.is_int_type), sizeof(bool)))
        return std::nullopt;
    if(!is.read(reinterpret_cast<char*>(&pvmd.is_data_type_signed), sizeof(bool)))
        return std::nullopt;
    if(!is.read(reinterpret_cast<char*>(&pvmd.data_type_size_in_bytes),
            sizeof(unsigned long)))
        return std::nullopt;
    if(!is.read(reinterpret_cast<char*>(&pvmd.number_of_components),
            sizeof(unsigned long)))
        return std::nullopt;
    if(!is.read(reinterpret_cast<char*>(&pvmd.number_of_tuples),
            sizeof(unsigned long)))
        return std::nullopt;
    return std::optional<PropertyVectorMetaData>(pvmd);
}

struct PropertyVectorPartitionMetaData
{
    unsigned long offset;
    unsigned long number_of_tuples;
};

inline void writePropertyVectorPartitionMetaData(
    std::ostream& os, PropertyVectorPartitionMetaData const& pvpmd)
{
    os.write(reinterpret_cast<char*>(
                 &const_cast<PropertyVectorPartitionMetaData&>(pvpmd)
                      .offset),
             sizeof(unsigned long));
    os.write(reinterpret_cast<char*>(
                 &const_cast<PropertyVectorPartitionMetaData&>(pvpmd)
                      .number_of_tuples),
             sizeof(unsigned long));
}

inline std::optional<PropertyVectorPartitionMetaData>
readPropertyVectorPartitionMetaData(std::istream& is)
{
    PropertyVectorPartitionMetaData pvpmd;
    if (!is.read(reinterpret_cast<char*>(&pvpmd.offset), sizeof(unsigned long)))
    {
        return std::optional<PropertyVectorPartitionMetaData>();
    }
    if (!is.read(reinterpret_cast<char*>(&pvpmd.number_of_tuples),
                 sizeof(unsigned long)))
    {
        return std::optional<PropertyVectorPartitionMetaData>();
    }
    return std::optional<PropertyVectorPartitionMetaData>(pvpmd);
}
}  // namespace IO
}  // namespace MeshLib
