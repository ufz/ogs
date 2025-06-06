/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include "AbstractJacobianAssembler.h"

namespace ProcessLib
{
//! Assembles the Jacobian matrix using forward differences.
class ForwardDifferencesJacobianAssembler final
    : public AbstractJacobianAssembler
{
public:
    //! Constructs a new instance.
    //!
    //! \param absolute_epsilons perturbations of the components of the local
    //! solution vector used for evaluating the finite differences.
    //!
    //! \note The size of \c absolute_epsilons defines the "number of
    //! components" of the local solution vector (This is not the number of
    //! elements of the vector!). Therefore the size of the local solution
    //! vector must be divisible by the size of \c absolute_epsilons. This is
    //! the only consistency check performed. It is not checked whether said
    //! "number of components" is sensible. E.g., one could pass one epsilon per
    //! node, which would be valid but would not make sense at all.
    explicit ForwardDifferencesJacobianAssembler(
        std::vector<double>&& absolute_epsilons);

    //! Assembles the Jacobian, the matrices \f$M\f$ and \f$K\f$, and the vector
    //! \f$b\f$.
    //! For the assembly the assemble() method of the given \c local_assembler
    //! is called several times and the Jacobian is built from finite
    //! differences.
    //! The number of calls of the assemble() method is \f$2N+1\f$ if \f$N\f$ is
    //! the size of \c local_x.
    //!
    //! \attention It is assumed that the local vectors and matrices are ordered
    //! by component.
    void assembleWithJacobian(LocalAssemblerInterface& local_assembler,
                              double const t, double const dt,
                              std::vector<double> const& local_x_data,
                              std::vector<double> const& local_x_prev_data,
                              std::vector<double>& local_b_data,
                              std::vector<double>& local_Jac_data) override;

    std::unique_ptr<AbstractJacobianAssembler> copy() const override;

private:
    std::vector<double> const _absolute_epsilons;

    // temporary data only stored here in order to avoid frequent memory
    // reallocations.
    std::vector<double> _local_M_data;
    std::vector<double> _local_K_data;
    std::vector<double> _local_b_data;
    std::vector<double> _local_x_perturbed_data;
};

}  // namespace ProcessLib
