/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   LiquidFlowLocalAssembler.h
 *
 * Created on August 19, 2016, 2:28 PM
 */

#pragma once

#include <map>
#include <unordered_map>
#include <vector>
#include <typeindex>

#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "NumLib/Extrapolation/ExtrapolatableElement.h"
#include "NumLib/DOF/DOFTableUtil.h"
#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "ProcessLib/LocalAssemblerInterface.h"
#include "ProcessLib/LocalAssemblerTraits.h"
#include "ProcessLib/Parameter/Parameter.h"
#include "ProcessLib/Utils/InitShapeMatrices.h"

#include "LiquidFlowMaterialProperties.h"

namespace ProcessLib
{
struct CoupledSolutionsForStaggeredScheme;

namespace LiquidFlow
{
const unsigned NUM_NODAL_DOF = 1;

class LiquidFlowLocalAssemblerInterface
    : public ProcessLib::LocalAssemblerInterface,
      public NumLib::ExtrapolatableElement
{
public:
    LiquidFlowLocalAssemblerInterface(
        CoupledSolutionsForStaggeredScheme* const coupling_term)
        : _coupling_term(coupling_term)
    {
    }

    void setCoupledSolutionsForStaggeredScheme(std::size_t const /*mesh_item_id*/,
                                  CoupledSolutionsForStaggeredScheme* const coupling_term)
    {
        _coupling_term = coupling_term;
    }

    virtual std::vector<double> const& getIntPtDarcyVelocity(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& /*cache*/) const = 0;

protected:
    // TODO: remove _coupling_term or move integration point data from local
    // assembler class to a new class to make local assembler unique for each
    //process.
    /// Pointer that is set from a Process class.
    CoupledSolutionsForStaggeredScheme* _coupling_term;
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

    using MatrixOfVelocityAtIntegrationPoints = Eigen::Map<
        Eigen::Matrix<double, GlobalDim, Eigen::Dynamic, Eigen::RowMajor>>;

public:
    LiquidFlowLocalAssembler(
        MeshLib::Element const& element,
        std::size_t const /*local_matrix_size*/,
        bool const is_axially_symmetric,
        unsigned const integration_order,
        int const gravitational_axis_id,
        double const gravitational_acceleration,
        double const reference_temperature,
        LiquidFlowMaterialProperties const& material_propertries,
        CoupledSolutionsForStaggeredScheme* coupling_term)
        : LiquidFlowLocalAssemblerInterface(coupling_term),
          _element(element),
          _integration_method(integration_order),
          _shape_matrices(initShapeMatrices<ShapeFunction, ShapeMatricesType,
                                            IntegrationMethod, GlobalDim>(
              element, is_axially_symmetric, _integration_method)),
          _gravitational_axis_id(gravitational_axis_id),
          _gravitational_acceleration(gravitational_acceleration),
          _reference_temperature(reference_temperature),
          _material_properties(material_propertries)
    {
    }

    void assemble(double const t, std::vector<double> const& local_x,
                  std::vector<double>& local_M_data,
                  std::vector<double>& local_K_data,
                  std::vector<double>& local_b_data) override;

    void assembleWithCoupledTerm(
        double const t, std::vector<double> const& local_x,
        std::vector<double>& local_M_data, std::vector<double>& local_K_data,
        std::vector<double>& local_b_data,
        LocalCoupledSolutions const& coupled_term) override;

    Eigen::Map<const Eigen::RowVectorXd> getShapeMatrix(
        const unsigned integration_point) const override
    {
        auto const& N = _shape_matrices[integration_point].N;

        // assumes N is stored contiguously in memory
        return Eigen::Map<const Eigen::RowVectorXd>(N.data(), N.size());
    }

    std::vector<double> const& getIntPtDarcyVelocity(
        const double t,
        GlobalVector const& current_solution,
        NumLib::LocalToGlobalIndexMap const& dof_table,
        std::vector<double>& veloctiy_cache) const override;

private:
    MeshLib::Element const& _element;

    IntegrationMethod const _integration_method;
    std::vector<ShapeMatrices, Eigen::aligned_allocator<ShapeMatrices>>
        _shape_matrices;

    /**
     *  Calculator of the Laplacian and the gravity term for anisotropic
     *  permeability tensor
     */
    struct AnisotropicCalculator
    {
        static void calculateLaplacianAndGravityTerm(
            Eigen::Map<NodalMatrixType>& local_K,
            Eigen::Map<NodalVectorType>& local_b, ShapeMatrices const& sm,
            Eigen::MatrixXd const& permeability,
            double const integration_factor, double const mu,
            double const rho_g, int const gravitational_axis_id);

        static void calculateVelocity(
            unsigned const ip, Eigen::Map<const NodalVectorType> const& local_p,
            ShapeMatrices const& sm, Eigen::MatrixXd const& permeability,
            double const mu, double const rho_g,
            int const gravitational_axis_id,
            MatrixOfVelocityAtIntegrationPoints& darcy_velocity_at_ips);
    };

    /**
     *  Calculator of the Laplacian and the gravity term for isotropic
     *  permeability tensor
     */
    struct IsotropicCalculator
    {
        static void calculateLaplacianAndGravityTerm(
            Eigen::Map<NodalMatrixType>& local_K,
            Eigen::Map<NodalVectorType>& local_b, ShapeMatrices const& sm,
            Eigen::MatrixXd const& permeability,
            double const integration_factor, double const mu,
            double const rho_g, int const gravitational_axis_id);

        static void calculateVelocity(
            unsigned const ip, Eigen::Map<const NodalVectorType> const& local_p,
            ShapeMatrices const& sm, Eigen::MatrixXd const& permeability,
            double const mu, double const rho_g,
            int const gravitational_axis_id,
            MatrixOfVelocityAtIntegrationPoints& darcy_velocity_at_ips);
    };

    template <typename LaplacianGravityVelocityCalculator>
    void assembleMatrixAndVector(const int material_id, double const t,
                                 std::vector<double> const& local_x,
                                 std::vector<double>& local_M_data,
                                 std::vector<double>& local_K_data,
                                 std::vector<double>& local_b_data,
                                 SpatialPosition const& pos,
                                 Eigen::MatrixXd const& permeability);

    template <typename LaplacianGravityVelocityCalculator>
    void assembleWithCoupledWithHeatTransport(
        const int material_id, double const t, double const dt,
        std::vector<double> const& local_x, std::vector<double> const& local_T0,
        std::vector<double> const& local_T1, std::vector<double>& local_M_data,
        std::vector<double>& local_K_data, std::vector<double>& local_b_data,
        SpatialPosition const& pos, Eigen::MatrixXd const& permeability);

    void computeDarcyVelocity(
        Eigen::MatrixXd const& permeability,
        std::vector<double> const& local_x,
        MatrixOfVelocityAtIntegrationPoints& darcy_velocity_at_ips) const;

    void computeDarcyVelocityWithCoupling(
        Eigen::MatrixXd const& permeability, std::vector<double> const& local_x,
        std::unordered_map<std::type_index, const std::vector<double>> const&
            coupled_local_solutions,
        MatrixOfVelocityAtIntegrationPoints& darcy_velocity_at_ips) const;

    template <typename LaplacianGravityVelocityCalculator>
    void computeDarcyVelocityLocal(
        std::vector<double> const& local_x,
        Eigen::MatrixXd const& permeability,
        MatrixOfVelocityAtIntegrationPoints& darcy_velocity_at_ips) const;

    template <typename LaplacianGravityVelocityCalculator>
    void computeDarcyVelocityCoupledWithHeatTransportLocal(
        std::vector<double> const& local_x,
        std::vector<double> const& local_T,
        Eigen::MatrixXd const& permeability,
        MatrixOfVelocityAtIntegrationPoints& darcy_velocity_at_ips) const;

    const int _gravitational_axis_id;
    const double _gravitational_acceleration;
    const double _reference_temperature;
    const LiquidFlowMaterialProperties& _material_properties;
};

}  // end of namespace
}  // end of namespace

#include "LiquidFlowLocalAssembler-impl.h"
