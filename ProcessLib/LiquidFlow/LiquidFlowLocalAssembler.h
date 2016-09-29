/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   LiquidFlowLocalAssembler.h
 *
 * Created on August 19, 2016, 2:28 PM
 */

#ifndef OGS_LIQUIDFLOWLOCALASSEMBLER_H
#define OGS_LIQUIDFLOWLOCALASSEMBLER_H

#include <vector>

#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "NumLib/Extrapolation/ExtrapolatableElement.h"
#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "ProcessLib/LocalAssemblerInterface.h"
#include "ProcessLib/LocalAssemblerTraits.h"
#include "ProcessLib/Parameter/Parameter.h"
#include "ProcessLib/Utils/InitShapeMatrices.h"

#include "LiquidFlowMaterialProperties.h"

namespace ProcessLib
{
namespace LiquidFlow
{
const unsigned NUM_NODAL_DOF = 1;

class LiquidFlowLocalAssemblerInterface
    : public ProcessLib::LocalAssemblerInterface,
      public NumLib::ExtrapolatableElement
{
public:
    virtual std::vector<double> const& getIntPtDarcyVelocityX(
        std::vector<double>& /*cache*/) const = 0;

    virtual std::vector<double> const& getIntPtDarcyVelocityY(
        std::vector<double>& /*cache*/) const = 0;

    virtual std::vector<double> const& getIntPtDarcyVelocityZ(
        std::vector<double>& /*cache*/) const = 0;
};

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
class LiquidFlowLocalAssembler : public LiquidFlowLocalAssemblerInterface
{
    using ShapeMatricesType = ShapeMatrixPolicyType<ShapeFunction, GlobalDim>;
    using ShapeMatrices = typename ShapeMatricesType::ShapeMatrices;

    using LocalAssemblerTraits = ProcessLib::LocalAssemblerTraits<
        ShapeMatricesType, ShapeFunction::NPOINTS, NUM_NODAL_DOF, GlobalDim>;

    using NodalMatrixType = typename LocalAssemblerTraits::LocalMatrix;
    using NodalVectorType = typename LocalAssemblerTraits::LocalVector;
    using GlobalDimVectorType = typename ShapeMatricesType::GlobalDimVectorType;

public:
    LiquidFlowLocalAssembler(MeshLib::Element const& element,
                             std::size_t const /*local_matrix_size*/,
                             bool const is_axially_symmetric,
                             unsigned const integration_order,
                             int const gravitational_axis_id,
                             double const gravitational_acceleration,
                             LiquidFlowMaterialProperties& material_propertries)
        : _element(element),
          _integration_method(integration_order),
          _shape_matrices(initShapeMatrices<ShapeFunction, ShapeMatricesType,
                                            IntegrationMethod, GlobalDim>(
              element, is_axially_symmetric, _integration_method)),
          _gravitational_axis_id(gravitational_axis_id),
          _gravitational_acceleration(gravitational_acceleration),
          _material_properties(material_propertries)
    {
    }

    void assemble(double const t, std::vector<double> const& local_x,
                  std::vector<double>& local_M_data,
                  std::vector<double>& local_K_data,
                  std::vector<double>& local_b_data) override;

    void computeSecondaryVariable(std::vector<double> const& local_x) override;

    Eigen::Map<const Eigen::RowVectorXd> getShapeMatrix(
        const unsigned integration_point) const override
    {
        auto const& N = _shape_matrices[integration_point].N;

        // assumes N is stored contiguously in memory
        return Eigen::Map<const Eigen::RowVectorXd>(N.data(), N.size());
    }

    std::vector<double> const& getIntPtDarcyVelocityX(
        std::vector<double>& /*cache*/) const override
    {
        assert(_darcy_velocities.size() > 0);
        return _darcy_velocities[0];
    }

    std::vector<double> const& getIntPtDarcyVelocityY(
        std::vector<double>& /*cache*/) const override
    {
        assert(_darcy_velocities.size() > 1);
        return _darcy_velocities[1];
    }

    std::vector<double> const& getIntPtDarcyVelocityZ(
        std::vector<double>& /*cache*/) const override
    {
        assert(_darcy_velocities.size() > 2);
        return _darcy_velocities[2];
    }

private:
    MeshLib::Element const& _element;

    IntegrationMethod const _integration_method;
    std::vector<ShapeMatrices> _shape_matrices;

    std::vector<std::vector<double>> _darcy_velocities =
        std::vector<std::vector<double>>(
            GlobalDim,
            std::vector<double>(_integration_method.getNumberOfPoints()));

    /**
     *  Calculator of the Laplacian and the gravity term for anisotropic
     *  permeability tensor
     */
    struct AnisotropicLaplacianAndGravityTermCalculator
    {
        static void calculate(Eigen::Map<NodalMatrixType>& local_K,
                              Eigen::Map<NodalVectorType>& local_b,
                              ShapeMatrices const& sm,
                              Eigen::MatrixXd const& perm,
                              double const integration_factor, double const mu,
                              double const rho_g,
                              int const gravitational_axis_id);
    };

    /**
     *  Calculator of the Laplacian and the gravity term for isotropic
     *  permeability tensor
     */
    struct IsotropicLaplacianAndGravityTermCalculator
    {
        static void calculate(Eigen::Map<NodalMatrixType>& local_K,
                              Eigen::Map<NodalVectorType>& local_b,
                              ShapeMatrices const& sm,
                              Eigen::MatrixXd const& perm,
                              double const integration_factor, double const mu,
                              double const rho_g,
                              int const gravitational_axis_id);
    };

    template <typename LaplacianAndGravityTermCalculator>
    void local_assemble(double const t, std::vector<double> const& local_x,
                        std::vector<double>& local_M_data,
                        std::vector<double>& local_K_data,
                        std::vector<double>& local_b_data,
                        SpatialPosition const& pos,
                        Eigen::MatrixXd const& perm);

    const int _gravitational_axis_id;
    const double _gravitational_acceleration;
    LiquidFlowMaterialProperties& _material_properties;
    double _temperature;
};

}  // end of namespace
}  // end of namespace

#include "LiquidFlowLocalAssembler-impl.h"

#endif /* LIQUIDFLOWLOCALASSEMBLER_H */
