// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "FormattingUtils.h"

namespace MathLib
{
Eigen::IOFormat const EigenIOFormat::full_precision{
    Eigen::FullPrecision, 0, " ", "\n", "", "", "", "", ' '};
}  // namespace MathLib
