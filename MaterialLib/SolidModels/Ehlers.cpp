/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "Ehlers.h"
#include "Ehlers-impl.h"

namespace MaterialLib
{
namespace Solids
{
namespace Ehlers
{
template class SolidEhlers<2>;
template class SolidEhlers<3>;

template <>
ProcessLib::KelvinMatrixType<3> sOdotS<3>(
    ProcessLib::KelvinVectorType<3> const& v)
{
    ProcessLib::KelvinMatrixType<3> result;

    result(0, 0) = v(0) * v(0);
    result(0, 1) = result(1, 0) = v(3) * v(3) / 2.;
    result(0, 2) = result(2, 0) = v(5) * v(5) / 2.;
    result(0, 3) = result(3, 0) = v(0) * v(3);
    result(0, 4) = result(4, 0) = v(3) * v(5) / std::sqrt(2.);
    result(0, 5) = result(5, 0) = v(0) * v(5);

    result(1, 1) = v(1) * v(1);
    result(1, 2) = result(2, 1) = v(4) * v(4) / 2.;
    result(1, 3) = result(3, 1) = v(3) * v(1);
    result(1, 4) = result(4, 1) = v(1) * v(4);
    result(1, 5) = result(5, 1) = v(3) * v(4) / std::sqrt(2.);

    result(2, 2) = v(2) * v(2);
    result(2, 3) = result(3, 2) = v(5) * v(4) / std::sqrt(2.);
    result(2, 4) = result(4, 2) = v(4) * v(2);
    result(2, 5) = result(5, 2) = v(5) * v(2);

    result(3, 3) = v(0) * v(1) + v(3) * v(3) / 2.;
    result(3, 4) = result(4, 3) =
        v(3) * v(4) / 2. + v(5) * v(1) / std::sqrt(2.);
    result(3, 5) = result(5, 3) =
        v(0) * v(4) / std::sqrt(2.) + v(3) * v(5) / 2.;

    result(4, 4) = v(1) * v(2) + v(4) * v(4) / 2.;
    result(4, 5) = result(5, 4) =
        v(3) * v(2) / std::sqrt(2.) + v(5) * v(4) / 2.;

    result(5, 5) = v(0) * v(2) + v(5) * v(5) / 2.;
    return result;
}

template <>
ProcessLib::KelvinMatrixType<2> sOdotS<2>(
    ProcessLib::KelvinVectorType<2> const& v)
{
    ProcessLib::KelvinMatrixType<2> result;

    result(0, 0) = v(0) * v(0);
    result(0, 1) = result(1, 0) = v(3) * v(3) / 2.;
    result(0, 2) = result(2, 0) = 0;
    result(0, 3) = result(3, 0) = v(0) * v(3);

    result(1, 1) = v(1) * v(1);
    result(1, 2) = result(2, 1) = 0;
    result(1, 3) = result(3, 1) = v(3) * v(1);

    result(2, 2) = v(2) * v(2);
    result(2, 3) = result(3, 2) = 0;

    result(3, 3) = v(0) * v(1) + v(3) * v(3) / 2.;

    return result;
}

}  // namespace Ehlers
}  // namespace Solids
}  // namespace MaterialLib
