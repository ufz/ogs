/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "FormattingUtils.h"

namespace MathLib
{
Eigen::IOFormat const EigenIOFormat::full_precision{
    Eigen::FullPrecision, 0, " ", "\n", "", "", "", "", ' '};
}  // namespace MathLib
