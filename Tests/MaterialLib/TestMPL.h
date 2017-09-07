/**
 * \file
 * \author Norbert Grunwald
 * \date   Oct 22, 2018
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "Tests/TestTools.h"

#include "MaterialLib/MPL/mpMedium.h"
namespace MPL = MaterialPropertyLib;

MPL::Medium createTestMaterial(std::string const& xml);
