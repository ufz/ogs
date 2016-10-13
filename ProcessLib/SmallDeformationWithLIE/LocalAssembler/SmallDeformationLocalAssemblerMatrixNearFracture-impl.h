/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_SMALLDEFORMATION_WITH_LIE_SMALLDEFORMATIONLOCALASSEMBLER_MATRIX_IMPL_H_
#define PROCESSLIB_SMALLDEFORMATION_WITH_LIE_SMALLDEFORMATIONLOCALASSEMBLER_MATRIX_IMPL_H_

#include "SmallDeformationLocalAssemblerMatrixNearFracture.h"

#include <valarray>
#include <vector>

#include <Eigen/Eigen>

#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "MathLib/Point3d.h"

#include "MeshLib/Node.h"

#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"


#include "ProcessLib/Deformation/BMatrixPolicy.h"
#include "ProcessLib/Deformation/LinearBMatrix.h"
#include "ProcessLib/Utils/InitShapeMatrices.h"

#include "ProcessLib/SmallDeformationWithLIE/Common/LevelSetFunction.h"
#include "ProcessLib/SmallDeformationWithLIE/Common/Utils.h"

#include "IntegrationPointDataMatrix.h"
#include "SecondaryData.h"
#include "SmallDeformationLocalAssemblerInterface.h"

namespace ProcessLib
{
namespace SmallDeformationWithLIE
{
template <typename ShapeFunction, typename IntegrationMethod,
          int DisplacementDim>
SmallDeformationLocalAssemblerMatrixNearFracture<ShapeFunction, IntegrationMethod,
                               DisplacementDim>::
    SmallDeformationLocalAssemblerMatrixNearFracture(
        MeshLib::Element const& e,
        std::size_t const n_variables,
        std::size_t const /*local_matrix_size*/,
        std::vector<unsigned> const& dofIndex_to_localIndex,
        bool is_axially_symmetric,
        unsigned const integration_order,
        SmallDeformationProcessData<DisplacementDim>& process_data)
    : SmallDeformationLocalAssemblerInterface(
          n_variables * ShapeFunction::NPOINTS * DisplacementDim,
          dofIndex_to_localIndex),
      _process_data(process_data),
      _integration_method(integration_order),
      _shape_matrices(
          initShapeMatrices<ShapeFunction, ShapeMatricesType,
                            IntegrationMethod, DisplacementDim>(
              e, is_axially_symmetric, _integration_method)),
      _element(e)
{
    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    _ip_data.reserve(n_integration_points);
    _secondary_data.N.resize(n_integration_points);

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        _ip_data.emplace_back(*_process_data._material);
        auto& ip_data = _ip_data[ip];
        auto const& sm = _shape_matrices[ip];
        ip_data._detJ = sm.detJ;
        ip_data._integralMeasure = sm.integralMeasure;
        ip_data._b_matrices.resize(
            KelvinVectorDimensions<DisplacementDim>::value,
            ShapeFunction::NPOINTS * DisplacementDim);

        auto const x_coord =
            interpolateXCoordinate<ShapeFunction, ShapeMatricesType>(e,
                                                                     sm.N);
        LinearBMatrix::computeBMatrix<DisplacementDim,
                                      ShapeFunction::NPOINTS>(
            sm.dNdx, ip_data._b_matrices, is_axially_symmetric, sm.N,
            x_coord);

        ip_data._sigma.resize(KelvinVectorDimensions<DisplacementDim>::value);
        ip_data._sigma_prev.resize(
            KelvinVectorDimensions<DisplacementDim>::value);
        ip_data._eps.resize(KelvinVectorDimensions<DisplacementDim>::value);
        ip_data._eps_prev.resize(
            KelvinVectorDimensions<DisplacementDim>::value);
        ip_data._C.resize(KelvinVectorDimensions<DisplacementDim>::value,
                          KelvinVectorDimensions<DisplacementDim>::value);

        _secondary_data.N[ip] = sm.N;
    }
}


template <typename ShapeFunction, typename IntegrationMethod,
          int DisplacementDim>
void SmallDeformationLocalAssemblerMatrixNearFracture<ShapeFunction, IntegrationMethod,
                                    DisplacementDim>::
assembleWithJacobian(
    double const t,
    Eigen::VectorXd const& local_u,
    Eigen::VectorXd& local_b,
    Eigen::MatrixXd& local_J)
{
    assert (_element.getDimension() == DisplacementDim);

    auto const n_dof_per_var = ShapeFunction::NPOINTS * DisplacementDim;

    auto localRhs_ru = local_b.segment(0, n_dof_per_var);
    auto localRhs_du = local_b.segment(n_dof_per_var, n_dof_per_var);

    auto localA_uu = local_J.block(0, 0, n_dof_per_var, n_dof_per_var);
    auto localA_ug = local_J.block(0, n_dof_per_var, n_dof_per_var, n_dof_per_var);
    auto localA_gu = local_J.block(n_dof_per_var, 0, n_dof_per_var, n_dof_per_var);
    auto localA_gg = local_J.block(n_dof_per_var, n_dof_per_var, n_dof_per_var, n_dof_per_var);

    auto const nodal_ru = local_u.segment(0, n_dof_per_var);
    auto nodal_du = local_u.segment(n_dof_per_var, n_dof_per_var);

    auto const &fracture_props = *_process_data._fracture_property;

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    SpatialPosition x_position;
    x_position.setElementID(_element.getID());

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);

        auto const& sm = _shape_matrices[ip];
        auto &ip_data = _ip_data[ip];
        auto const& wp = _integration_method.getWeightedPoint(ip);
        auto const& detJ = ip_data._detJ;
        auto const& integralMeasure = ip_data._integralMeasure;

        // levelset functions
        auto const ip_physical_coords = computePhysicalCoordinates(_element, sm.N);
        double levelsets = calculateLevelSetFunction(fracture_props, ip_physical_coords.getCoords());

        // nodal displacement = u^hat + levelset* [u]
        NodalDisplacementVectorType nodal_u = nodal_ru + levelsets * nodal_du;

        // strain, stress
        auto const& B = ip_data._b_matrices;
        auto const& eps_prev = ip_data._eps_prev;
        auto const& sigma_prev = ip_data._sigma_prev;

        auto& eps = ip_data._eps;
        auto& sigma = ip_data._sigma;
        auto& C = ip_data._C;
        auto& material_state_variables = *ip_data._material_state_variables;

        eps.noalias() = B * nodal_u;

        if (!ip_data._solid_material.computeConstitutiveRelation(
                t, x_position, _process_data.dt, eps_prev, eps, sigma_prev,
                sigma, C, material_state_variables))
            OGS_FATAL("Computation of local constitutive relation failed.");

        // r_ru = B^T Stress = B^T C B (u + phi*[u])
        localRhs_ru.noalias() -= B.transpose() * sigma * detJ * wp.getWeight() * integralMeasure;

        // r_[u] = (phi*B)^T Stress = (phi*B)^T C B (u + phi*[u])
        localRhs_du.noalias() -= levelsets * B.transpose() * sigma * detJ * wp.getWeight() * integralMeasure;

        // J_uu += B^T * C * B
        localA_uu.noalias() += B.transpose() * C * B * detJ * wp.getWeight() * integralMeasure;

        // J_u[u] += B^T * C * B * levelset
        localA_ug.noalias() += B.transpose() * C * levelsets * B * detJ * wp.getWeight() * integralMeasure;

        // J_[u]u += (levelset B)^T * C * B
        localA_gu.noalias() += levelsets * B.transpose() * C * B * detJ * wp.getWeight() * integralMeasure;

        // J_[u][u] += (levelset B)^T * C * (levelset B)
        localA_gg.noalias() += levelsets * B.transpose() * C * levelsets * B * detJ * wp.getWeight() * integralMeasure;
    }
}


template <typename ShapeFunction, typename IntegrationMethod,
          int DisplacementDim>
void SmallDeformationLocalAssemblerMatrixNearFracture<ShapeFunction, IntegrationMethod,
                                    DisplacementDim>::
postTimestepConcrete(std::vector<double> const& /*local_x*/)
{
    const int n = 3;
    std::valarray<double> ele_stress(0.0, n);
    std::valarray<double> ele_strain(0.0, n);
    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        auto& ip_data = _ip_data[ip];

        ele_stress[0] += ip_data._sigma[0];
        ele_stress[1] += ip_data._sigma[1];
        ele_stress[2] += ip_data._sigma[3];

        ele_strain[0] += ip_data._eps[0];
        ele_strain[1] += ip_data._eps[1];
        ele_strain[2] += ip_data._eps[3];
    }
    ele_stress /= n_integration_points;
    ele_strain /= n_integration_points;
    (*_process_data._mesh_prop_stress_xx)[_element.getID()] = ele_stress[0];
    (*_process_data._mesh_prop_stress_yy)[_element.getID()] = ele_stress[1];
    (*_process_data._mesh_prop_stress_xy)[_element.getID()] = ele_stress[2];
    (*_process_data._mesh_prop_strain_xx)[_element.getID()] = ele_strain[0];
    (*_process_data._mesh_prop_strain_yy)[_element.getID()] = ele_strain[1];
    (*_process_data._mesh_prop_strain_xy)[_element.getID()] = ele_strain[2];
}

}  // namespace SmallDeformationWithLIE
}  // namespace ProcessLib

#endif // PROCESSLIB_SMALLDEFORMATION_WITH_LIE_SMALLDEFORMATIONLOCALASSEMBLER_MATRIX_IMPL_H_
