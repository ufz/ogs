// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <Eigen/LU>
#include <array>
#include <optional>

#include "ComputeNaturalCoordsRootFindingProblem.h"
#include "NumLib/NewtonRaphson.h"

namespace ApplicationUtils
{
/**
 *
 * Computes natural coordinates inside a given mesh element for some passed real
 * coordinates.
 *
 * This base class has implementations templated for each shape function, see
 * ComputeNaturalCoordsSolverImplementation below.
 *
 */
class ComputeNaturalCoordsSolverInterface
{
public:
    virtual Eigen::Vector3d solve(MeshLib::Element const& e,
                                  Eigen::Vector3d const& real_coords,
                                  int const max_iter,
                                  double const real_coords_tolerance) const = 0;

    virtual ~ComputeNaturalCoordsSolverInterface() = default;
};

template <typename ShapeFunction>
class ComputeNaturalCoordsSolverImplementation
    : public ComputeNaturalCoordsSolverInterface
{
    static constexpr int ElementDim = ShapeFunction::DIM;
    using ShapeMatricesType = ShapeMatrixPolicyType<ShapeFunction, ElementDim>;

public:
    Eigen::Vector3d solve(MeshLib::Element const& e,
                          Eigen::Vector3d const& real_coords,
                          int const max_iter,
                          double const real_coords_tolerance) const override
    {
        using Problem = ComputeNaturalCoordsRootFindingProblem<ShapeFunction>;
        using LJM = typename Problem::LocalJacobianMatrix;
        using LRV = typename Problem::LocalResidualVector;

        Problem problem(e, Eigen::Map<const LRV>(real_coords.data()));

        Eigen::PartialPivLU<LJM> linear_solver(ElementDim);
        LJM jacobian;

        const double increment_tolerance = 0;
        auto const newton_solver = NumLib::makeNewtonRaphson(
            linear_solver,
            [&problem](LJM& J) { problem.updateJacobian(J); },
            [&problem](LRV& res) { problem.updateResidual(res); },
            [&problem](auto& delta_r) { problem.updateSolution(delta_r); },
            {max_iter, real_coords_tolerance, increment_tolerance});

        auto const opt_iter = newton_solver.solve(jacobian);

        if (opt_iter)
        {
            DBUG("Newton solver succeeded after {} iterations", *opt_iter);

            Eigen::Vector3d natural_coords = Eigen::Vector3d::Zero();
            natural_coords.head<ElementDim>() = problem.getNaturalCoordinates();
            return natural_coords;
        }

        OGS_FATAL(
            "Newton solver failed. Please consider increasing the error "
            "tolerance or the max. number of Newton iterations.");
    }
};
}  // namespace ApplicationUtils
