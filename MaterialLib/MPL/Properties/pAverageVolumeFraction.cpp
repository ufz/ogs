/**
 * \author Norbert Grunwald
 * \date   13.09.2017
 * \brief  
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "pAverageVolumeFraction.h"

namespace MaterialPropertyLib
{
/// This constructor initializes the private medium pointer
/// using the provided argument.
AverageVolumeFraction::AverageVolumeFraction(Medium* m)
: _medium (m){};
/// This constructor throws an error, since the property is not
/// implemented on the phase scale.
AverageVolumeFraction::AverageVolumeFraction(Phase*)
: _medium (0)
{
	OGS_FATAL("The 'AverageVolumeFraction' property is not "
			"available on the Phase scale");
};
/// This constructor throws an error, since the property is not
/// implemented on the component scale.
AverageVolumeFraction::AverageVolumeFraction(Component*)
: _medium (0)
{
	OGS_FATAL("The 'AverageVolumeFraction' property is not "
			"available on the Component scale");
}

/**
 * This method computes an average property by
 * \f$\xi=\sum_{\alpha=0}^{n}\phi_\alpha\xi_\alpha\f$
 * where \f$\xi\f$ is an arbitrary property, \f$\alpha\f$ is a phase and
 * \f$\phi_\alpha\f$ is the volume fraction of that phase.
 */
PropertyDataType AverageVolumeFraction::value (VariableArray const&)
{
	/// \todo: implementation of AverageVolumeFraction.
	return _value;
}
} //MaterialPropertyLib

