/**
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "MathLib/FormattingUtils.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Node.h"
#include "NumLib/Fem/InitShapeMatrices.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"

namespace ApplicationUtils
{
/**
 * Models the computation of natural coordinates inside a given mesh element for
 * some given real coordinates as a root finding problem.
 *
 * That root finding problem is designed such that it can be easily used/solved
 * by the NumLib::NewtonRaphson solver.
 */
template <typename ShapeFunction>
class ComputeNaturalCoordsRootFindingProblem
{
    static constexpr int Dim = ShapeFunction::DIM;
    using ShapeMatricesType = ShapeMatrixPolicyType<ShapeFunction, Dim>;

public:
    using LocalJacobianMatrix =
        Eigen::Matrix<double, Dim, Dim, Eigen::RowMajor>;
    //! Used for the residuum ("res"), the solution vector ("r", natural
    //! coordinates) and real coordinates ("x").
    using LocalResidualVector = Eigen::Matrix<double, Dim, 1>;

    ComputeNaturalCoordsRootFindingProblem(
        MeshLib::Element const& e, LocalResidualVector const& x_expected)
        : e_{e},
          node_coords_{getNodeCoords(e)},
          x_expected_{x_expected},
          // TODO 2D elements with 3D coordinates etc.
          // TODO check all dims
          sm_{ShapeFunction::DIM, e.getDimension(), ShapeFunction::NPOINTS}
    {
        updateShp();
    }

    void updateJacobian(LocalJacobianMatrix& J) const
    {
        auto const& dNdr = sm_.dNdr;

        for (int comp = 0; comp < Dim; ++comp)
        {
            J.row(comp) = node_coords_.row(comp) * dNdr.transpose();
        }
    }

    void updateResidual(LocalResidualVector& res) const
    {
        res = node_coords_ * sm_.N.transpose() - x_expected_;
    }

    //! Updates the current guess.
    //!
    //! \param delta_r the solution increment (natural coordinates)
    void updateSolution(LocalResidualVector const& delta_r)
    {
        r_ += delta_r;
        updateShp();
    }

    LocalResidualVector const& getNaturalCoordinates() const { return r_; }

private:
    using NodeCoordsMatrix =
        Eigen::Matrix<double, Dim, ShapeFunction::NPOINTS, Eigen::RowMajor>;

    void updateShp()
    {
        auto const fe =
            NumLib::createIsoparametricFiniteElement<ShapeFunction,
                                                     ShapeMatricesType>(e_);

        constexpr auto GlobalDim = Dim;
        fe.template computeShapeFunctions<NumLib::ShapeMatrixType::N_J>(
            r_.data(), sm_, GlobalDim, false);
    }

    static NodeCoordsMatrix getNodeCoords(MeshLib::Element const& e)
    {
        NodeCoordsMatrix node_coords;

        for (std::size_t n = 0; n < ShapeFunction::NPOINTS; ++n)
        {
            node_coords.col(n) = Eigen::Map<const Eigen::Vector<double, Dim>>(
                e.getNode(n)->data());
        }

        return node_coords;
    }

    static constexpr LocalResidualVector initialGuess()
    {
        return LocalResidualVector::Constant(0.5);
    }

    MeshLib::Element const& e_;
    NodeCoordsMatrix const node_coords_;
    LocalResidualVector const x_expected_;

    typename ShapeMatricesType::ShapeMatrices sm_;
    LocalResidualVector r_ = initialGuess();
};
}  // namespace ApplicationUtils
