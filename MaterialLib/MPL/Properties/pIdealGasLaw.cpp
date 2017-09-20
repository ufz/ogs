/**
 * \author Norbert Grunwald
 * \date   18.09.2017
 * \brief  
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "pIdealGasLaw.h"

namespace MaterialPropertyLib
{
IdealGasLaw::IdealGasLaw(Medium*)
: _phase (0),
  _component (0)
{
    OGS_FATAL("The 'AverageMoleFraction' property is currently not "
            "implemented on the Medium scale");
}

IdealGasLaw::IdealGasLaw(Phase* p)
: _phase (p),
  _component (0)
{};

IdealGasLaw::IdealGasLaw(Component* c)
: _phase (0),
  _component (c){};

/**
 */
PropertyDataType IdealGasLaw::value (VariableArray const&)
{
    /// \todo: implementation of IdealGasLaw.
    return _value;
}
} //MaterialPropertyLib

