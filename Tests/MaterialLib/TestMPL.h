/**
 * \file
 * \author Norbert Grunwald
 * \date   Oct 22, 2018
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "Tests/TestTools.h"

#include "MaterialLib/MPL/Medium.h"
namespace MPL = MaterialPropertyLib;

std::unique_ptr<MPL::Medium> createTestMaterial(std::string const& xml);
