/**
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_GROUNDWATERFLOW_FEM_H_
#define PROCESS_LIB_GROUNDWATERFLOW_FEM_H_

#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"

namespace ProcessLib
{

namespace GroundwaterFlow
{

template <typename GlobalMatrix, typename GlobalVector>
class LocalAssemblerDataInterface
{
public:
    virtual ~LocalAssemblerDataInterface() = default;

    virtual void init(MeshLib::Element const& e,
            std::size_t const local_matrix_size,
            double const hydraulic_conductivity) = 0;

    virtual void assemble() = 0;

    virtual void addToGlobal(GlobalMatrix& A, GlobalVector& rhs,
            AssemblerLib::LocalToGlobalIndexMap::RowColumnIndices const&) const = 0;
};

template <typename ShapeFunction_,
         typename IntegrationMethod_,
         unsigned IntegrationOrder_,
         typename GlobalMatrix,
         typename GlobalVector>
class LocalAssemblerData : public LocalAssemblerDataInterface<GlobalMatrix, GlobalVector>
{
public:
    using ShapeFunction = ShapeFunction_;
    using NodalMatrixType = typename ShapeMatrixPolicyType<ShapeFunction>::NodalMatrixType;
    using NodalVectorType = typename ShapeMatrixPolicyType<ShapeFunction>::NodalVectorType;

    using ShapeMatrices = typename ShapeMatrixPolicyType<ShapeFunction>::ShapeMatrices;


    static unsigned constexpr integration_order = IntegrationOrder_;
    static unsigned constexpr n_integration_points =
        MathLib::pow(integration_order, ShapeFunction::DIM);
    using IntegrationMethod = IntegrationMethod_;

    /// The hydraulic_conductivity factor is directly integrated into the local
    /// element matrix.
    void
    init(MeshLib::Element const& e,
        std::size_t const local_matrix_size,
        double const hydraulic_conductivity)
    {
        using FemType = NumLib::TemplateIsoparametric<
            ShapeFunction, ShapeMatrices>;

        FemType fe(*static_cast<const typename ShapeFunction::MeshElement*>(&e));


        IntegrationMethod _integration_method(integration_order);

        for (std::size_t ip(0); ip < n_integration_points; ip++) {
            fe.computeShapeFunctions(
                    _integration_method.getWeightedPoint(ip).getCoords(),
                    _shape_matrices[ip]);
        }

        _hydraulic_conductivity = hydraulic_conductivity;

        localA.reset(new NodalMatrixType(local_matrix_size, local_matrix_size));
        localRhs.reset(new NodalVectorType(local_matrix_size));
    }

    void assemble()
    {
        localA->setZero();
        localRhs->setZero();

        IntegrationMethod _integration_method(integration_order);

        for (std::size_t ip(0); ip < n_integration_points; ip++) {
            auto const& sm = _shape_matrices[ip];
            auto const& wp = _integration_method.getWeightedPoint(ip);
            *localA += sm.dNdx.transpose() * _hydraulic_conductivity *
                        sm.dNdx * sm.detJ * wp.getWeight();
        }
    }

    void addToGlobal(GlobalMatrix& A, GlobalVector& rhs,
            AssemblerLib::LocalToGlobalIndexMap::RowColumnIndices const& indices) const
    {
        A.add(indices, *localA);
        rhs.add(indices.rows, *localRhs);
    }

private:
    std::array<ShapeMatrices, n_integration_points> _shape_matrices;
    double _hydraulic_conductivity;

    std::unique_ptr<NodalMatrixType> localA;
    std::unique_ptr<NodalVectorType> localRhs;
};


}   // namespace GroundwaterFlow
}   // namespace ProcessLib

#endif  // PROCESS_LIB_GROUNDWATERFLOW_FEM_H_
