/**
 * \file
 * \author Tobias Meisel
 * \date   2020-12-08
 * \brief  Dispatches functions specific to execution platform (w/o MPI)
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <utility>

namespace MeshLib::IO
{
std::pair<std::size_t, std::size_t> getPartitionInfo(std::size_t const size);
bool isFileManager();
}  // namespace MeshLib::IO
