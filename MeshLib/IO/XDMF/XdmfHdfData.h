/**
 * \file
 * \author Tobias Meisel
 * \date   2020-11-13
 * \brief  Holds all data for the combined writing of xdmf and hdf5 file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <boost/shared_ptr.hpp>
#include <vector>

#include "HdfData.h"
#include "XdmfData.h"

// TODO (tm) This structs are the result of transformData.cpp. This structures
// can be eliminated.
namespace MeshLib::IO
{
struct Geometry final
{
    std::vector<double> flattened_values;
    HdfData hdf;
    XdmfData xdmf;
};

struct Topology final
{
    std::vector<int> flattened_values;
    HdfData hdf;
    XdmfData xdmf;
};

struct AttributeMeta final
{
    HdfData hdf;
    XdmfData xdmf;
};
}  // namespace MeshLib::IO