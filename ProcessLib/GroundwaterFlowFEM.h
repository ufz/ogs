/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_GROUNDWATERFLOW_FEM_H_
#define PROCESS_LIB_GROUNDWATERFLOW_FEM_H_

#include <memory>
#include <vector>

#include "logog/include/logog.hpp"

#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"

#include "GroundwaterFlowMaterialProperty.h"

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
            const GroundwaterFlowMaterialProperty &mat,
            unsigned const integration_order) = 0;

    virtual void assemble() = 0;

    virtual void addToGlobal(GlobalMatrix& A, GlobalVector& rhs,
            AssemblerLib::LocalToGlobalIndexMap::RowColumnIndices const&) const = 0;
};

template <typename ShapeFunction_,
         typename IntegrationMethod_,
         typename GlobalMatrix,
         typename GlobalVector>
class LocalAssemblerData : public LocalAssemblerDataInterface<GlobalMatrix, GlobalVector>
{
public:
    using ShapeFunction = ShapeFunction_;
    using NodalMatrixType = typename ShapeMatrixPolicyType<ShapeFunction>::NodalMatrixType;
    using NodalVectorType = typename ShapeMatrixPolicyType<ShapeFunction>::NodalVectorType;

    using ShapeMatricesType = ShapeMatrixPolicyType<ShapeFunction>;
    using ShapeMatrices = typename ShapeMatricesType::ShapeMatrices;

    /// The hydraulic_conductivity factor is directly integrated into the local
    /// element matrix.
    void
    init(MeshLib::Element const& e,
        std::size_t const local_matrix_size,
		const GroundwaterFlowMaterialProperty &mat,
        unsigned const integration_order)
    {
        using FemType = NumLib::TemplateIsoparametric<
            ShapeFunction, ShapeMatricesType>;

        FemType fe(*static_cast<const typename ShapeFunction::MeshElement*>(&e));


        _integration_order = integration_order;
        IntegrationMethod_ integration_method(_integration_order);
        std::size_t const n_integration_points = integration_method.getNPoints();

        _shape_matrices.resize(n_integration_points);
        for (std::size_t ip(0); ip < n_integration_points; ip++) {
            _shape_matrices[ip].resize(ShapeFunction::DIM, ShapeFunction::NPOINTS);
            fe.computeShapeFunctions(
                    integration_method.getWeightedPoint(ip).getCoords(),
                    _shape_matrices[ip]);
        }

        // mat_group_id could be get by using MeshLib::Element const& e in future.
        std::size_t const mat_group_id =0;
        _hydraulic_conductivity = mat.getConductivity(mat_group_id);
        _storage = mat.getStorage(mat_group_id);

        _localA.reset(new NodalMatrixType(local_matrix_size, local_matrix_size));
        _localRhs.reset(new NodalVectorType(local_matrix_size));
    }

    void assemble()
    {
        _localA->setZero();
        _localRhs->setZero();

        IntegrationMethod_ integration_method(_integration_order);

        switch(_hydraulic_conductivity->getType())
        {
            case MaterialLib::PermeabilityType::ISOTROPIC:
            {
                using IsoHydraulicConductivity
                          = MaterialLib::TensorParameter<MaterialLib::PermeabilityType,
				                             MaterialLib::ConstantScalarModel, double>;
                const IsoHydraulicConductivity *scalarK
                          = reinterpret_cast<IsoHydraulicConductivity*>(_hydraulic_conductivity);

                assemblyLaplacian(integration_method, *scalarK);
            }
            break;
            case MaterialLib::PermeabilityType::ANISOTROPIC:
            {
                using AnistropicHydraulicConductivity
                          = MaterialLib::TensorParameter<MaterialLib::PermeabilityType,
                                                         MaterialLib::ConstantTensor<Matrix>, Matrix>;
                const AnistropicHydraulicConductivity *anisK
                          = reinterpret_cast<AnistropicHydraulicConductivity*>(_hydraulic_conductivity);

                assemblyLaplacian(integration_method, *anisK);
            }
            break;
            default:
            {
                WARN("No hydraulic conductivity model is provided for this element.");
            }
        }
    }

    void addToGlobal(GlobalMatrix& A, GlobalVector& rhs,
            AssemblerLib::LocalToGlobalIndexMap::RowColumnIndices const& indices) const
    {
        A.add(indices, *_localA);
        rhs.add(indices.rows, *_localRhs);
    }

private:
    std::vector<ShapeMatrices> _shape_matrices;

    using Permeability = MaterialLib::ParameterBase<MaterialLib::PermeabilityType>;
    Permeability* _hydraulic_conductivity;

    using Storage = MaterialLib::ParameterBase<MaterialLib::StorageType>;
    Storage* _storage;

    std::unique_ptr<NodalMatrixType> _localA;
    std::unique_ptr<NodalVectorType> _localRhs;

    unsigned _integration_order = 2;

    template<typename LaplaceCoefficient>
        void assemblyLaplacian(IntegrationMethod_ &integration_method,
                               const LaplaceCoefficient & laplace_coefficient)
    {
        unsigned const n_integration_points = integration_method.getNPoints();

        for (std::size_t ip(0); ip < n_integration_points; ip++) {
            auto const& sm = _shape_matrices[ip];
            auto const& wp = integration_method.getWeightedPoint(ip);
            _localA->noalias() += sm.dNdx.transpose() * laplace_coefficient.getValue() *
                                  sm.dNdx * sm.detJ * wp.getWeight();
        }
    }
};

}   // namespace GroundwaterFlow
}   // namespace ProcessLib

#endif  // PROCESS_LIB_GROUNDWATERFLOW_FEM_H_
