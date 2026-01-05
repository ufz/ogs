// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <boost/shared_ptr.hpp>
#include <vector>

#include "HdfData.h"
#include "XdmfData.h"

namespace MeshLib::IO
{
struct XdmfHdfData final
{
    HdfData hdf;
    XdmfData xdmf;
};
}  // namespace MeshLib::IO