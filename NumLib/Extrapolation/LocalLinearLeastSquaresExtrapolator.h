/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

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
 * The number of integration points in each element must be greater than or
 * equal to the number of nodes of that element. This restriction is due
 * to the use of the least squares which requires an exact or overdetermined
 * equation system.
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

    void extrapolate(const unsigned num_components,
                     ExtrapolatableElementCollection const& extrapolatables,
                     const double t,
                     GlobalVector const& current_solution,
                     LocalToGlobalIndexMap const& dof_table) override;

    /*! \copydoc Extrapolator::calculateResiduals()
     *
     * The computed residuals are root-mean-square of the difference between
     * the integration point values obtained from the local assemblers and the
     * extrapolation results when interpolated back to the integration points
     * again.
     */
    void calculateResiduals(
        const unsigned num_components,
        ExtrapolatableElementCollection const& extrapolatables,
        const double t,
        GlobalVector const& current_solution,
        LocalToGlobalIndexMap const& dof_table) override;

    GlobalVector const& getNodalValues() const override
    {
        return *_nodal_values;
    }

    GlobalVector const& getElementResiduals() const override
    {
        return *_residuals;
    }

private:
    //! Extrapolate one element.
    void extrapolateElement(
        std::size_t const element_index, const unsigned num_components,
        ExtrapolatableElementCollection const& extrapolatables, const double t,
        GlobalVector const& current_solution,
        LocalToGlobalIndexMap const& dof_table, GlobalVector& counts);

    //! Compute the residuals for one element
    void calculateResidualElement(
        std::size_t const element_index,
        const unsigned num_components,
        ExtrapolatableElementCollection const& extrapolatables,
        const double t,
        GlobalVector const& current_solution,
        LocalToGlobalIndexMap const& dof_table);

    std::unique_ptr<GlobalVector> _nodal_values;  //!< extrapolated nodal values
    std::unique_ptr<GlobalVector> _residuals;     //!< extrapolation residuals

    //! DOF table used for writing to global vectors.
    NumLib::LocalToGlobalIndexMap const& _dof_table_single_component;

    //! Avoids frequent reallocations.
    std::vector<double> _integration_point_values_cache;

    //! Stores a matrix and its Moore-Penrose pseudo-inverse.
    struct CachedData
    {
        //! The matrix A.
        Eigen::MatrixXd A;

        //! Moore-Penrose pseudo-inverse of A.
        Eigen::MatrixXd A_pinv;
    };

    /*! Maps (\#nodes, \#int_pts) to (N_0, QR decomposition),
     * where N_0 is the shape matrix of the first integration point.
     *
     * \note It is assumed that the pair (\#nodes, \#int_pts) uniquely
     * identifies the set of all shape matrices N for a mesh element (i.e., only
     * N, not dN/dx etc.).
     *
     * \todo Add the element dimension as identifying criterion, or change to
     * typeid.
     */
    std::map<std::pair<unsigned, unsigned>, CachedData> _qr_decomposition_cache;
};

}  // namespace NumLib
