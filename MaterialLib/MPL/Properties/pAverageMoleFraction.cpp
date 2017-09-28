/**
 * \author Norbert Grunwald
 * \date   12.09.2017
 * \brief
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "pAverageMoleFraction.h"

namespace MaterialPropertyLib
{
AverageMoleFraction::AverageMoleFraction(Medium*) : _phase(0)
{
    notImplemented("AverageMoleFraction", "Medium");
}

AverageMoleFraction::AverageMoleFraction(Phase* p) : _phase(p){};

AverageMoleFraction::AverageMoleFraction(Component*) : _phase(0)
{
    notImplemented("AverageMoleFraction", "Component");
}
/**
 * This method computes an average phase property by
 * \f$\xi^\alpha=\sum_{\zeta=0}^{n}x_{n\zeta}^\alpha\xi_\alpha\f$
 * where \f$\xi\f$ is an arbitrary property, \f$\alpha\f$ is the phase,
 * \f$x_{n\zeta}^\alpha\f$ is the mole fraction of the component \f$\zeta\f$.
 */
PropertyDataType AverageMoleFraction::value(VariableArray const&)
{
    if (isUpdated())
        return _value;

    /// \todo: implementation of AverageMoleFraction.
    return _value;
}
}  // MaterialPropertyLib
