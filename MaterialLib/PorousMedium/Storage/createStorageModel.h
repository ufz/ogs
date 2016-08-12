/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file:   createStorageModel.h
 *
 * Created on August 16, 2016, 1:16 PM
 */

#ifndef CREATESTORAGEMODEL_H
#define CREATESTORAGEMODEL_H

#include "Storage.h"

#include "BaseLib/ConfigTree.h"

namespace MaterialLib
{
namespace PorousMedium
{
/// Create a storage model
/// \param config  ConfigTree object has a tag of <storage>
Storage* createStorageModel(BaseLib::ConfigTree const* const config);

}  // end namespace
}  // end namespace

#endif /* CREATESTORAGEMODEL_H */
