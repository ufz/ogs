/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file:   createPorosityModel.h
 *
 * Created on August 16, 2016, 1:16 PM
 */

#ifndef CREATEPOROSITYMODEL_H
#define CREATEPOROSITYMODEL_H

#include <memory>

#include "Porosity.h"

#include "BaseLib/ConfigTree.h"

namespace MaterialLib
{
namespace PorousMedium
{
/// Create a porosity model
/// \param config  ConfigTree object has a tag of <porosity>
std::unique_ptr<Porosity> createPorosityModel(BaseLib::ConfigTree const& config);

}  // end namespace
}  // end namespace

#endif /* CREATEPOROSITYMODEL_H */
