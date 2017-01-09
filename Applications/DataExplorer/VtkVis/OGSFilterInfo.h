/**
 * \file
 * \author Karsten Rink
 * \date   2010-09-29
 * \brief  Definition of the OGSFilterInfo class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef OGSFILTERINFO_H
#define OGSFILTERINFO_H

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

#endif // OGSFILTERINFO_H
