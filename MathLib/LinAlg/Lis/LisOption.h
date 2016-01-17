/**
 * \file
 * \author Norihiro Watanabe
 * \date   2012-06-25
 * \brief  Definition of the LisOption class.
 *
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef LIS_OPTION_H_
#define LIS_OPTION_H_

#include <string>
#include <map>

#include <logog/include/logog.hpp>

#include "BaseLib/ConfigTreeNew.h"
#include "MathLib/LinAlg/LinearSolverOptions.h"

namespace MathLib
{
/** Option for Lis solver.
 *
 * Options are stored as a string obtained from the lis-tag of the passed
 * options.  The string will be passed as is to Lis.  In particular it will not
 * be checked by this method if options are valid.
 *
 * For possible values, please refer to the Lis User Guide:
 * http://www.ssisc.org/lis/lis-ug-en.pdf
 *
 * Note: Option -omp_num_threads cannot be used with this class since Lis
 * currently (version 1.5.57) only sets the number of threads in
 * \c lis_initialize(). Refer to the Lis source code for details.
 */
struct LisOption
{
	LisOption(BaseLib::ConfigTreeNew const* const options)
	{
		if (options)
		{
			ignoreOtherLinearSolvers(*options, "lis");
			if (auto s = options->getConfParamOptional<std::string>("lis")) {
				if (!s->empty()) _option_string += " " + *s;
			}
		}
#ifndef NDEBUG
		INFO("Lis options: \"%s\"", _option_string.c_str());
#endif
	}

	std::string _option_string = "-initxzeros 0";

};
}
#endif  // LIS_OPTION_H_
