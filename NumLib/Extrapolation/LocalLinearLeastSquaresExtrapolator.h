/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef NUMLIB_LOCAL_LLSQ_EXTRAPOLATOR_H
#define NUMLIB_LOCAL_LLSQ_EXTRAPOLATOR_H

#include <map>

#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "NumLib/DOF/GlobalMatrixProviders.h"
#include "Extrapolator.h"

namespace NumLib
{
/*! Extrapolator doing a linear least squares extrapolation locally in each
 *  element.
 *
 * The results of the element-wise least squares are averaged at each node.
 *
 * \note
 * \parblock
 * The described procedure is not least-squares optimal on the whole mesh! But
 * the residuals are computed from the actual interpolation result that is being
 * returned by this class.
 *
 * Currently this class only supports interpolating single component variables.
 *
 * Furthermore, the number of integration points in each element must be greater
 * than or equal to the number of nodes of that element. This restriction is due to
 * the use of the least squares which requires an exact or overdetermined equation
 * system.
 * \endparblock
 */
class LocalLinearLeastSquaresExtrapolator : public Extrapolator
{
public:
    /*! Constructs a new instance.
     *
     * \note
     * The \c dof_table must point to a d.o.f. table for one single-component
     * variable.
     */
    explicit LocalLinearLeastSquaresExtrapolator(
        NumLib::LocalToGlobalIndexMap const& dof_table);

    void extrapolate(
            ExtrapolatableElementCollection const& extrapolatables) override;

    /*! \copydoc Extrapolator::calculateResiduals()
     *
     * The computed residuals are root-mean-square of the difference between
     * the integration point values obtained from the local assemblers and the
     * extrapolation results when interpolated back to the integration points
     * again.
     */
    void calculateResiduals(
            ExtrapolatableElementCollection const& extrapolatables) override;

    GlobalVector const& getNodalValues() const override
    {
        return _nodal_values;
    }

    GlobalVector const& getElementResiduals() const override
    {
        return _residuals;
    }

    ~LocalLinearLeastSquaresExtrapolator()
    {
        NumLib::GlobalVectorProvider::provider.releaseVector(
            _nodal_values);
    }

private:
    //! Extrapolate one element.
    void extrapolateElement(
        std::size_t const element_index,
        ExtrapolatableElementCollection const& extrapolatables,
        GlobalVector& counts);

    //! Compute the residuals for one element
    void calculateResidualElement(
        std::size_t const element_index,
        ExtrapolatableElementCollection const& extrapolatables);

    GlobalVector& _nodal_values;  //!< extrapolated nodal values
    GlobalVector _residuals;      //!< extrapolation residuals

    //! DOF table used for writing to global vectors.
    NumLib::LocalToGlobalIndexMap const& _local_to_global;

    //! Avoids frequent reallocations.
    std::vector<double> _integration_point_values_cache;

    //! Stores a QR decomposition of some matrix.
    struct CachedData {
        explicit CachedData(Eigen::Ref<const Eigen::RowVectorXd> N_0_)
            : N_0(N_0_)
        {
        }

        //! Zeroth shape matrix. Used for assertion only.
        Eigen::RowVectorXd const N_0;

        //! Moore-Penrose pseudo-inverse of the nodal value to integration point
        //! value interpolation matrix.
        Eigen::MatrixXd p_inv;
    };

    /*! Maps (#nodes, #int_pts) to (N_0, QR decomposition),
     * where N_0 is the shape matrix of the first integration point.
     *
     * \note It is assumed that the pair (#nodes, #int_pts) uniquely identifies
     * the set of all shape matrices N for a mesh element (i.e., only N, not
     * dN/dx etc.).
     *
     * \todo Add the element dimension as identifying criterion, or change to
     * typeid.
     */
    std::map<std::pair<unsigned, unsigned>, CachedData> _qr_decomposition_cache;
};

}  // namespace NumLib

#endif  // NUMLIB_LOCAL_LLSQ_EXTRAPOLATOR_H
