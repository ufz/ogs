/**
 * \file
 * \author Norbert Grunwald
 * \date   Sep 7, 2017
 *
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */
#pragma once

#include "MaterialLib/MPL/Component.h"

namespace MaterialPropertyLib
{
/// A class for Water derived from Component.
///
/// This class can holds material constants and default properties of ordinary
/// water.
struct Water final : public Component
{
    explicit Water(std::unique_ptr<PropertyArray>&& properties);
};
}  // namespace MaterialPropertyLib
