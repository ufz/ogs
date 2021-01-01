/**
 * \file
 * \author Tobias Meisel
 * \date   2020-11-13
 * \brief  Definition of the data layer for writing Meshes
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <boost/shared_ptr.hpp>
#include <string>
#include <vector>

class XdmfAttributeCenter;
class XdmfArrayType;

namespace MeshLib::IO
{
using XdmfDimType = unsigned int;
using Hdf5DimType = unsigned long long;
struct Geometry final
{
    std::vector<double> flattend_values;
    std::vector<XdmfDimType> const starts;
    std::vector<XdmfDimType> const strides;
    std::vector<XdmfDimType> const vdims;
    std::vector<Hdf5DimType> const vldims;
};

struct Topology final
{
    std::vector<int> flattend_values;
    std::vector<XdmfDimType> const starts;
    std::vector<XdmfDimType> const strides;
    std::vector<XdmfDimType> const vdims;
    std::vector<Hdf5DimType> const vldims;
};

struct AttributeMeta final
{
    const void* data_start;
    std::string name;
    boost::shared_ptr<const XdmfAttributeCenter> attribute_center;
    boost::shared_ptr<const XdmfArrayType> data_type;
    std::vector<XdmfDimType> starts;
    std::vector<XdmfDimType> strides;
    std::vector<XdmfDimType> vdims;
    std::vector<Hdf5DimType> vldims;
};
}  // namespace MeshLib::IO