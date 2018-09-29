/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file:   createStorageModel.h
 *
 * Created on August 16, 2016, 1:16 PM
 */

#pragma once

#include <memory>

namespace BaseLib
{
class ConfigTree;
}

namespace MaterialLib
{
namespace PorousMedium
{
class Storage;
/** Create a storage model
 *  @param config  ConfigTree object has a tag of `<storage>`
 */
std::unique_ptr<Storage> createStorageModel(BaseLib::ConfigTree const& config);

}  // end namespace
}  // end namespace
