/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>
#include <vector>

#include "ProcessLib/LocalAssemblerInterface.h"
#include "NumLib/Extrapolation/ExtrapolatableElement.h"
#include "TESAssemblyParams.h"
#include "TESLocalAssemblerInner-fwd.h"

namespace ProcessLib
{
namespace TES
{
class TESLocalAssemblerInterface
    : public ProcessLib::LocalAssemblerInterface,
      public NumLib::ExtrapolatableElement
{
public:
    virtual ~TESLocalAssemblerInterface() = default;

    virtual bool checkBounds(std::vector<double> const& local_x,
                             std::vector<double> const& local_x_prev_ts) = 0;

    virtual std::vector<double> const& getIntPtSolidDensity(
        std::vector<double>& /*cache*/) const = 0;

    virtual std::vector<double> const& getIntPtLoading(
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtReactionDampingFactor(
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtReactionRate(
        std::vector<double>& /*cache*/) const = 0;

    virtual std::vector<double> const& getIntPtDarcyVelocityX(
        std::vector<double>& /*cache*/) const = 0;

    virtual std::vector<double> const& getIntPtDarcyVelocityY(
        std::vector<double>& /*cache*/) const = 0;

    virtual std::vector<double> const& getIntPtDarcyVelocityZ(
        std::vector<double>& /*cache*/) const = 0;
};

template <typename ShapeFunction_, typename IntegrationMethod_,
          unsigned GlobalDim>
class TESLocalAssembler final
    : public TESLocalAssemblerInterface
{
public:
    using ShapeFunction = ShapeFunction_;
    using ShapeMatricesType = ShapeMatrixPolicyType<ShapeFunction, GlobalDim>;
    using ShapeMatrices = typename ShapeMatricesType::ShapeMatrices;

    TESLocalAssembler(MeshLib::Element const& e,
                      std::size_t const local_matrix_size,
                      bool is_axially_symmetric,
                      unsigned const integration_order,
                      AssemblyParams const& asm_params);

    void assemble(double const t, std::vector<double> const& local_x,
                  std::vector<double>& local_M_data,
                  std::vector<double>& local_K_data,
                  std::vector<double>& local_b_data) override;

    Eigen::Map<const Eigen::RowVectorXd> getShapeMatrix(
        const unsigned integration_point) const override
    {
        auto const& N = _shape_matrices[integration_point].N;

        // assumes N is stored contiguously in memory
        return Eigen::Map<const Eigen::RowVectorXd>(N.data(), N.size());
    }

    bool checkBounds(std::vector<double> const& local_x,
                     std::vector<double> const& local_x_prev_ts) override;

    std::vector<double> const& getIntPtSolidDensity(
        std::vector<double>& /*cache*/) const override;

    std::vector<double> const& getIntPtLoading(
        std::vector<double>& cache) const override;

    std::vector<double> const& getIntPtReactionDampingFactor(
        std::vector<double>& cache) const override;

    std::vector<double> const& getIntPtReactionRate(
        std::vector<double>& /*cache*/) const override;

    std::vector<double> const& getIntPtDarcyVelocityX(
        std::vector<double>& /*cache*/) const override;

    std::vector<double> const& getIntPtDarcyVelocityY(
        std::vector<double>& /*cache*/) const override;

    std::vector<double> const& getIntPtDarcyVelocityZ(
        std::vector<double>& /*cache*/) const override;
private:
    IntegrationMethod_ const _integration_method;

    std::vector<ShapeMatrices> _shape_matrices;

    using LAT = LocalAssemblerTraits<ShapeMatricesType, ShapeFunction::NPOINTS,
                                     NODAL_DOF, GlobalDim>;

    TESLocalAssemblerInner<LAT> _d;

    using NodalMatrixType = typename LAT::LocalMatrix;
    using NodalVectorType = typename LAT::LocalVector;

    static_assert(
        std::is_same<NodalMatrixType, typename LAT::LocalMatrix>::value,
        "local matrix and data traits matrix do not coincide");
    static_assert(
        std::is_same<NodalVectorType, typename LAT::LocalVector>::value,
        "local vector and data traits vector do not coincide");

};

}  // namespace TES
}  // namespace ProcessLib

#include "TESLocalAssembler-impl.h"
