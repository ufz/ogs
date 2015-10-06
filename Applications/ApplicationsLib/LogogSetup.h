/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef APPLICATIONSLIB_LOGOGSETUP_H_
#define APPLICATIONSLIB_LOGOGSETUP_H_

#include <logog/include/logog.hpp>

#include "BaseLib/LogogSimpleFormatter.h"

namespace ApplicationsLib
{

/// Initialization and shutting down of the logog library.
class LogogSetup final
{
public:

	LogogSetup()
	{
		LOGOG_INITIALIZE();
		fmt = new BaseLib::LogogSimpleFormatter;
		logog_cout = new logog::Cout;
		logog_cout->SetFormatter(*fmt);
	}

	~LogogSetup()
	{
		delete fmt;
		delete logog_cout;
		LOGOG_SHUTDOWN();
	}

private:
	BaseLib::LogogSimpleFormatter* fmt;
	logog::Cout* logog_cout;
};

}	// ApplicationsLib

#endif  // APPLICATIONSLIB_LOGOGSETUP_H_
