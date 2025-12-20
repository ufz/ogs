// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "NonLinearFbar.h"

namespace ProcessLib
{
namespace NonLinearFbar
{

BarDetFType convertStringToDetFBarType(
    std::string_view const bar_det_f_type_name)
{
    if (boost::iequals(bar_det_f_type_name, "element_center_value"))
    {
        INFO(
            "Use F bar method with the element center value of F for the "
            "det(F) modification.");
        return BarDetFType::ELEMENT_CENTER_VALUE;
    }
    if (boost::iequals(bar_det_f_type_name, "element_average"))
    {
        INFO(
            "Use F bar method with the element average of F for the det(F) "
            "modification.");
        return BarDetFType::ELEMENT_AVERAGE;
    }

    return BarDetFType::NONE;
}
}  // namespace NonLinearFbar
}  // namespace ProcessLib