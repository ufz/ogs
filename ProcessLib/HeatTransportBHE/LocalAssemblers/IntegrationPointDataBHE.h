/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <Eigen/Eigen>

#include "ProcessLib/HeatTransportBHE/BHE/BHEAbstract.h"

namespace ProcessLib
{
namespace HeatTransportBHE
{
using namespace BHE;

template <typename ShapeMatrixType>
struct IntegrationPointDataBHE final
{
    explicit IntegrationPointDataBHE(BHEAbstract& bhe_instance)
        : _bhe_instance(bhe_instance)
    {
        // depending on the type of BHE
        const int unknown_size = _bhe_instance.getNumUnknowns();
        // initialization
        _vec_mass_coefficients.resize(unknown_size);
        for (int i = 0; i < unknown_size; i++)
        {
            Eigen::MatrixXd mat_laplace(3, 3);
            mat_laplace.setZero();
            _vec_mat_Laplace.push_back(mat_laplace);
            Eigen::VectorXd vec_advection(3);
            vec_advection.setZero();
            _vec_Advection_vectors.push_back(vec_advection);
        }

        // set parameter values
        for (int j = 0; j < unknown_size; j++)
        {
            // mass matrix coefficients
            _vec_mass_coefficients[j] = _bhe_instance.getMassCoeff(j);
            // laplace matrix
            _bhe_instance.getLaplaceMatrix(j, _vec_mat_Laplace[j]);
            // advection vector
            _bhe_instance.getAdvectionVector(j, _vec_Advection_vectors[j]);
        }
    }

    BHEAbstract& _bhe_instance;

    typename ShapeMatrixType::NodalRowVectorType N;
    typename ShapeMatrixType::GlobalDimNodalMatrixType dNdx;
    double integration_weight;

    // mass coefficients, length depending on the type of BHE
    std::vector<double> _vec_mass_coefficients;

    // Laplace matrices
    std::vector<Eigen::MatrixXd> _vec_mat_Laplace;

    // Advection vectors
    std::vector<Eigen::VectorXd> _vec_Advection_vectors;

    void pushBackState() {}

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
