/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef NUMLIB_EXTRAPOLATOR_H
#define NUMLIB_EXTRAPOLATOR_H

#include <vector>

#include <Eigen/Eigen>
#include "NumLib/NumericsConfig.h"
#include "ExtrapolatableElementCollection.h"

namespace NumLib
{

//! Interface for classes that extrapolate integration point values to nodal
//! values.
class Extrapolator
{
public:
    //! Extrapolates the given \c property from the given local assemblers.
    virtual void extrapolate(
            ExtrapolatableElementCollection const& extrapolatables) = 0;

    /*! Computes residuals from the extrapolation of the given \c property.
     *
     * The residuals are computed as element values.
     *
     * \pre extrapolate() must have been called before with the same arguments.
     */
    virtual void calculateResiduals(
            ExtrapolatableElementCollection const& extrapolatables) = 0;

    //! Returns the extrapolated nodal values.
    //! \todo Maybe write directly to a MeshProperty.
    virtual GlobalVector const& getNodalValues() const = 0;

    //! Returns the extrapolation residuals.
    //! \todo Maybe write directly to a MeshProperty.
    virtual GlobalVector const& getElementResiduals() const = 0;

    virtual ~Extrapolator() = default;
};

} // namespace ProcessLib

#endif // NUMLIB_EXTRAPOLATOR_H
