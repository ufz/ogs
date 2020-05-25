/**
 * \file
 * \author Karsten Rink
 * \date   2010-09-29
 * \brief  Definition of the OGSFilterInfo class.
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "VtkOGSFilter.h"
#include <string>

///Stores information about filters in VtkOGSFilter for access-routines from the GUI.
class OGSFilterInfo
{
public:
    enum class VtkTargetObject
    {
        POLYDATA         = 0,
        UNSTRUCTUREDGRID = 1,
        IMAGEDATA        = 3
    };

    OGSFilterInfo(std::string t, VtkOGSFilter::OGSVisFilter f,
                  VtkTargetObject v) : text_(t), filter_(f), target_(v) {}
    ~OGSFilterInfo() {}
    const std::string& text() const { return text_; }
    const VtkOGSFilter::OGSVisFilter& filter() const { return filter_; }
    const VtkTargetObject& target() const { return target_; }

private:
    std::string text_;
    VtkOGSFilter::OGSVisFilter filter_;
    VtkTargetObject target_;
};
