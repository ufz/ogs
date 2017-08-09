/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "PhaseFieldProcess.h"

namespace ProcessLib
{
namespace PhaseField
{
template class PhaseFieldProcess<2>;
template class PhaseFieldProcess<3>;

/**
 * \brief A class to simulate mechanical fracturing process
 * using phase-field approach in solids described by
 *
 * \f[
 *     \mathrm{div} \left[ \left(d^2 + k \right) \boldsymbol{\sigma}_0^+
 *      + \boldsymbol{\sigma}_0^- \right] +
 *        \varrho \boldsymbol{b} = \boldsymbol{0}
 * \f]
 * \f[
 *     2d H(\boldsymbol{\epsilon}_\mathrm{el})
 *      - \frac{1 - d}{2 \varepsilon} g_\mathrm{c}
 *      - 2 \varepsilon g_\mathrm{c} \mathrm{div}(\mathrm{grad} d) = 0
 * \f]
 * where
 *    \f{eqnarray*}{
 *       &d:&                    \mbox{order parameter,}\\
 *       &\varrho:&              \mbox{density,}\\
 *       &H:&                    \mbox{history field,}\\
 *       &g_\mathrm{c}:&         \mbox{fracture energy,}\\
 *       &\varepsilon:&          \mbox{length scale}\\
 *    \f}
 *
 * Detailed model description can refer
 * <a href="Miao_Biot2017.pdff" target="_blank"><b>Phase field method</b></a>
 */

}  // namespace PhaseField
}  // namespace ProcessLib
