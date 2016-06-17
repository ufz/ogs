/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_TES_FEM_H_
#define PROCESS_LIB_TES_FEM_H_

#include <memory>
#include <vector>

#include "TESAssemblyParams.h"
#include "TESLocalAssemblerInner-fwd.h"

#include "NumLib/Extrapolation/Extrapolator.h"

namespace ProcessLib
{
namespace TES
{
class TESLocalAssemblerInterface
    : public NumLib::Extrapolatable<TESIntPtVariables>
{
public:
    virtual ~TESLocalAssemblerInterface() = default;

    virtual void assemble(double const t,
                          std::vector<double> const& local_x) = 0;

    virtual void addToGlobal(
        NumLib::LocalToGlobalIndexMap::RowColumnIndices const&,
        GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b) const = 0;

    virtual bool checkBounds(std::vector<double> const& local_x,
                             std::vector<double> const& local_x_prev_ts) = 0;
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
                      unsigned const integration_order,
                      AssemblyParams const& asm_params);

    void assemble(double const t, std::vector<double> const& local_x) override;

    void addToGlobal(
        NumLib::LocalToGlobalIndexMap::RowColumnIndices const& indices,
        GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b) const override;

    Eigen::Map<const Eigen::VectorXd> getShapeMatrix(
        const unsigned integration_point) const override
    {
        auto const& N = _shape_matrices[integration_point].N;

        // assumes N is stored contiguously in memory
        return Eigen::Map<const Eigen::VectorXd>(N.data(), N.size());
    }

    bool checkBounds(std::vector<double> const& local_x,
                     std::vector<double> const& local_x_prev_ts) override;

    std::vector<double> const& getIntegrationPointValues(
        TESIntPtVariables const var, std::vector<double>& cache) const override;

private:
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

    // TODO Change VectorMatrixAssembler s.t. these can be omitted.
    NodalMatrixType _local_M;
    NodalMatrixType _local_K;
    NodalVectorType _local_b;

    // TODO Use the value from Process
    unsigned const _integration_order;
};

}  // namespace TES
}  // namespace ProcessLib

#include "TESLocalAssembler-impl.h"

#endif  // PROCESS_LIB_TES_FEM_H_
