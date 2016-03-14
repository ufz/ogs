/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_TES_FEM_H_
#define PROCESS_LIB_TES_FEM_H_

#include <memory>
#include <vector>

#include "TESFEM-data-fwd.h"
#include "TESProcess-notpl.h"

#include "NumLib/Extrapolation/LocalNodalDOF.h"

namespace ProcessLib
{

namespace TES
{

class Extrapolatable
{
public:
    virtual Eigen::Map<const Eigen::VectorXd>
    getShapeMatrix(const unsigned integration_point) const = 0;

    virtual std::vector<double> const&
    getIntegrationPointValues(SecondaryVariables var, NumLib::LocalNodalDOF& nodal_dof) const = 0;
};


template <typename GlobalMatrix, typename GlobalVector>
class LocalAssemblerDataInterface
        : public Extrapolatable
{
public:
    virtual ~LocalAssemblerDataInterface() = default;

    virtual void init(MeshLib::Element const& e,
                      std::size_t const local_matrix_size,
                      unsigned const integration_order,
                      TESProcessInterface* process) = 0;

    virtual void assemble(double const t, std::vector<double> const& local_x) = 0;

    virtual void addToGlobal(AssemblerLib::LocalToGlobalIndexMap::RowColumnIndices const&,
            GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b) const = 0;

    virtual bool checkBounds(std::vector<double> const& localX,
                             std::vector<double> const& localX_pts) = 0;
};



template <typename ShapeFunction_,
          typename IntegrationMethod_,
          typename GlobalMatrix,
          typename GlobalVector,
          unsigned GlobalDim>
class LocalAssemblerData final
        : public LocalAssemblerDataInterface<GlobalMatrix, GlobalVector>
{
public:
    using ShapeFunction = ShapeFunction_;
    using ShapeMatricesType = ShapeMatrixPolicyType<ShapeFunction, GlobalDim>;
    using ShapeMatrices     = typename ShapeMatricesType::ShapeMatrices;

    void
    init(MeshLib::Element const& e,
         std::size_t const local_matrix_size,
         unsigned const integration_order,
         TESProcessInterface* process) override;

    void assemble(double const t, std::vector<double> const& local_x) override;

    void addToGlobal(AssemblerLib::LocalToGlobalIndexMap::RowColumnIndices const& indices,
            GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b) const override;

    Eigen::Map<const Eigen::VectorXd>
    getShapeMatrix(const unsigned integration_point) const override {
        auto const& N = _shape_matrices[integration_point].N;

        // assumes N is stored contiguously in memory
        return Eigen::Map<const Eigen::VectorXd>(N.data(), N.size());
    }

    bool checkBounds(std::vector<double> const& localX,
                     std::vector<double> const& localX_pts);

    std::vector<double> const&
    getIntegrationPointValues(SecondaryVariables var, NumLib::LocalNodalDOF& nodal_dof) const override;

private:
    std::vector<ShapeMatrices> _shape_matrices;

    using DT = DataTraits<ShapeMatricesType, ShapeFunction::NPOINTS, NODAL_DOF, GlobalDim>;

    LADataNoTpl<DT> _data;

    using NodalMatrixType = typename DT::LocalMatrix;
    using NodalVectorType = typename DT::LocalVector;

    static_assert(std::is_same<NodalMatrixType, typename DT::LocalMatrix>::value,
                  "local matrix and data traits matrix do not coincide");
    static_assert(std::is_same<NodalVectorType, typename DT::LocalVector>::value,
                  "local vector and data traits vector do not coincide");
    NodalMatrixType _localA;
    NodalVectorType _localRhs;

    unsigned _integration_order = 2;

    std::unique_ptr<std::vector<double> > _integration_point_values_cache;
};


}   // namespace TES
}   // namespace ProcessLib


#include "TESFEM-impl.h"


#endif  // PROCESS_LIB_TES_FEM_H_
