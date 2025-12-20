// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
                  VtkTargetObject v) : _text(t), _filter(f), _target(v) {}
    ~OGSFilterInfo() {}
    const std::string& text() const { return _text; }
    const VtkOGSFilter::OGSVisFilter& filter() const { return _filter; }
    const VtkTargetObject& target() const { return _target; }

private:
    std::string _text;
    VtkOGSFilter::OGSVisFilter _filter;
    VtkTargetObject _target;
};
