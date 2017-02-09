/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <string>
#include <boost/optional.hpp>

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
};

inline void writePropertyVectorMetaDataBinary(
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
    DBUG("size of name: %d", pvmd.property_name.length());
    DBUG("name: '%s'", pvmd.property_name.c_str());
    DBUG("is_int_data_type: %d", pvmd.is_int_type);
    DBUG("is_data_type_signed: %d", pvmd.is_data_type_signed);
    DBUG("data_type_size_in_bytes: %d", pvmd.data_type_size_in_bytes);
    DBUG("number of components: i%d", pvmd.number_of_components);
    DBUG("number of tuples: %d", pvmd.number_of_tuples);
}

inline boost::optional<PropertyVectorMetaData> readPropertyVectorMetaData(
    std::istream& is)
{
    // read the size of the name of the PropertyVector
    std::string::size_type s = 0;
    if (!is.read(reinterpret_cast<char*>(&s), sizeof(std::string::size_type)))
        return boost::optional<PropertyVectorMetaData>();

    PropertyVectorMetaData pvmd;
    char *dummy = new char[s];
    if (!is.read(dummy, s))
        return boost::none;
    pvmd.property_name = std::string(dummy, s);
    delete [] dummy;

    if(!is.read(reinterpret_cast<char*>(&pvmd.is_int_type), sizeof(bool)))
        return boost::none;
    if(!is.read(reinterpret_cast<char*>(&pvmd.is_data_type_signed), sizeof(bool)))
        return boost::none;
    if(!is.read(reinterpret_cast<char*>(&pvmd.data_type_size_in_bytes),
            sizeof(unsigned long)))
        return boost::none;
    if(!is.read(reinterpret_cast<char*>(&pvmd.number_of_components),
            sizeof(unsigned long)))
        return boost::none;
    if(!is.read(reinterpret_cast<char*>(&pvmd.number_of_tuples),
            sizeof(unsigned long)))
        return boost::none;
    return boost::optional<PropertyVectorMetaData>(pvmd);
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

inline boost::optional<PropertyVectorPartitionMetaData>
readPropertyVectorPartitionMetaData(std::istream& is)
{
    PropertyVectorPartitionMetaData pvpmd;
    if (!is.read(reinterpret_cast<char*>(&pvpmd.offset),
                 sizeof(unsigned long)))
        return boost::optional<PropertyVectorPartitionMetaData>();
    if (!is.read(reinterpret_cast<char*>(&pvpmd.number_of_tuples),
                 sizeof(unsigned long)))
        return boost::optional<PropertyVectorPartitionMetaData>();
    return boost::optional<PropertyVectorPartitionMetaData>(pvpmd);
}
}
}
