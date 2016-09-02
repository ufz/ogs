/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_HEATCONDUCTION_FEM_H_
#define PROCESS_LIB_HEATCONDUCTION_FEM_H_

#include <vector>

#include "NumLib/Extrapolation/ExtrapolatableElement.h"
#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "ProcessLib/LocalAssemblerInterface.h"
#include "ProcessLib/LocalAssemblerTraits.h"
#include "ProcessLib/Parameter/Parameter.h"
#include "ProcessLib/Utils/InitShapeMatrices.h"
#include "HeatConductionProcessData.h"

namespace ProcessLib
{
namespace HeatConduction
{
const unsigned NUM_NODAL_DOF = 1;

class HeatConductionLocalAssemblerInterface
    : public ProcessLib::LocalAssemblerInterface,
      public NumLib::ExtrapolatableElement
{
public:
    virtual std::vector<double> const& getIntPtHeatFluxX(
        std::vector<double>& /*cache*/) const = 0;

    virtual std::vector<double> const& getIntPtHeatFluxY(
        std::vector<double>& /*cache*/) const = 0;

    virtual std::vector<double> const& getIntPtHeatFluxZ(
        std::vector<double>& /*cache*/) const = 0;
};

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
class LocalAssemblerData : public HeatConductionLocalAssemblerInterface
{
    using ShapeMatricesType = ShapeMatrixPolicyType<ShapeFunction, GlobalDim>;
    using ShapeMatrices = typename ShapeMatricesType::ShapeMatrices;

    using LocalAssemblerTraits = ProcessLib::LocalAssemblerTraits<
        ShapeMatricesType, ShapeFunction::NPOINTS, NUM_NODAL_DOF, GlobalDim>;

    using NodalMatrixType = typename LocalAssemblerTraits::LocalMatrix;
    using NodalVectorType = typename LocalAssemblerTraits::LocalVector;
    using GlobalDimVectorType = typename ShapeMatricesType::GlobalDimVectorType;

public:
    /// The thermal_conductivity factor is directly integrated into the local
    /// element matrix.
    LocalAssemblerData(MeshLib::Element const& element,
                       std::size_t const local_matrix_size,
                       unsigned const integration_order,
                       HeatConductionProcessData const& process_data)
        : _element(element),
          _process_data(process_data),
          _localK(local_matrix_size,
                  local_matrix_size),  // TODO narrowing conversion
          _localM(local_matrix_size, local_matrix_size),
          _localRhs(local_matrix_size),
          _integration_method(integration_order),
          _shape_matrices(initShapeMatrices<ShapeFunction, ShapeMatricesType,
                                            IntegrationMethod, GlobalDim>(
              element, _integration_method))
    {
        // This assertion is valid only if all nodal d.o.f. use the same shape
        // matrices.
        assert(local_matrix_size == ShapeFunction::NPOINTS * NUM_NODAL_DOF);
    }

    void assembleConcrete(
        double const t, std::vector<double> const& local_x,
        NumLib::LocalToGlobalIndexMap::RowColumnIndices const& indices,
        GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b) override
    {
        _localK.setZero();
        _localM.setZero();
        _localRhs.setZero();

        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        SpatialPosition pos;
        pos.setElementID(_element.getID());

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            pos.setIntegrationPoint(ip);
            auto const& sm = _shape_matrices[ip];
            auto const& wp = _integration_method.getWeightedPoint(ip);
            auto const k = _process_data.thermal_conductivity(t, pos)[0];
            auto const heat_capacity = _process_data.heat_capacity(t, pos)[0];
            auto const density = _process_data.density(t, pos)[0];

            _localK.noalias() +=
                sm.dNdx.transpose() * k * sm.dNdx * sm.detJ * wp.getWeight();
            _localM.noalias() += sm.N.transpose() * density * heat_capacity *
                                 sm.N * sm.detJ * wp.getWeight();
            // heat flux only computed for output.
            GlobalDimVectorType const heat_flux =
                -k * sm.dNdx * Eigen::Map<const NodalVectorType>(
                                    local_x.data(), ShapeFunction::NPOINTS);

            for (unsigned d = 0; d < GlobalDim; ++d)
            {
                _heat_fluxes[d][ip] = heat_flux[d];
            }
        }

        K.add(indices, _localK);
        M.add(indices, _localM);
        b.add(indices.rows, _localRhs);
    }

    Eigen::Map<const Eigen::RowVectorXd> getShapeMatrix(
        const unsigned integration_point) const override
    {
        auto const& N = _shape_matrices[integration_point].N;

        // assumes N is stored contiguously in memory
        return Eigen::Map<const Eigen::RowVectorXd>(N.data(), N.size());
    }

    std::vector<double> const& getIntPtHeatFluxX(
        std::vector<double>& /*cache*/) const override
    {
        assert(_heat_fluxes.size() > 0);
        return _heat_fluxes[0];
    }

    std::vector<double> const& getIntPtHeatFluxY(
        std::vector<double>& /*cache*/) const override
    {
        assert(_heat_fluxes.size() > 1);
        return _heat_fluxes[1];
    }

    std::vector<double> const& getIntPtHeatFluxZ(
        std::vector<double>& /*cache*/) const override
    {
        assert(_heat_fluxes.size() > 2);
        return _heat_fluxes[2];
    }

private:
    MeshLib::Element const& _element;
    HeatConductionProcessData const& _process_data;

    NodalMatrixType _localK;
    NodalMatrixType _localM;
    NodalVectorType _localRhs;

    IntegrationMethod const _integration_method;
    std::vector<ShapeMatrices> _shape_matrices;

    std::vector<std::vector<double>> _heat_fluxes =
        std::vector<std::vector<double>>(
            GlobalDim, std::vector<double>(ShapeFunction::NPOINTS));
};

}  // namespace HeatConduction
}  // namespace ProcessLib

#endif  // PROCESS_LIB_HEATCONDUCTION_FEM_H_
