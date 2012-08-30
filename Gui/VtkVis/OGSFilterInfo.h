/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 * \file OGSFilterInfo.h
 *
 * Created on 2010-09-29 by Karsten Rink
 */

#ifndef OGSFILTERINFO_H
#define OGSFILTERINFO_H

#include "VtkOGSFilter.h"
#include <string>

///Stores information about filters in VtkOGSFilter for access-routines from the GUI.
class OGSFilterInfo
{
public:
	enum VtkTargetObject
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
