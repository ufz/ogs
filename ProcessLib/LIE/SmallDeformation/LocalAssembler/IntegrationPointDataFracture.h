/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_LIE_SMALLDEFORMATION_INTEGRATIONPOINTDATAFRACTURE_H_
#define PROCESSLIB_LIE_SMALLDEFORMATION_INTEGRATIONPOINTDATAFRACTURE_H_

#include <Eigen/Eigen>

#include "MaterialLib/FractureModels/FractureModelBase.h"

namespace ProcessLib
{
namespace LIE
{
namespace SmallDeformation
{

template <typename HMatricesType, int DisplacementDim>
struct IntegrationPointDataFracture final
{
    explicit IntegrationPointDataFracture(
        MaterialLib::Fracture::FractureModelBase<DisplacementDim>& fracture_material)
        : _fracture_material(fracture_material)
    {
    }

#if defined(_MSC_VER) && _MSC_VER < 1900
    // The default generated move-ctor is correctly generated for other
    // compilers.
    explicit IntegrationPointDataFracture(IntegrationPointDataFracture&& other)
        : _h_matrices(std::move(other._h_matrices)),
          _sigma(std::move(other._sigma)),
          _sigma_prev(std::move(other._sigma_prev)),
          _w(std::move(other._w)),
          _w_prev(std::move(other._w_prev)),
          _aperture(std::move(other._aperture)),
          _aperture_prev(std::move(other._aperture_prev)),
          _aperture0(std::move(other._aperture0)),
          _fracture_material(other._fracture_material),
          _C(std::move(other._C)),
          _detJ(std::move(other._detJ))
    {
    }
#endif  // _MSC_VER

    typename HMatricesType::HMatrixType _h_matrices;
    typename HMatricesType::ForceVectorType _sigma, _sigma_prev;
    typename HMatricesType::ForceVectorType _w, _w_prev;
    double _aperture = 0.0;
    double _aperture_prev = 0.0;
    double _aperture0 = 0.0;

    MaterialLib::Fracture::FractureModelBase<DisplacementDim>& _fracture_material;

    Eigen::MatrixXd _C;
    double _detJ = 0.0;
    double _integralMeasure;


    void pushBackState()
    {
        _w_prev = _w;
        _sigma_prev = _sigma;
        _aperture_prev = _aperture;
    }
};

}  // namespace SmallDeformation
}  // namespace LIE
}  // namespace ProcessLib

#endif // PROCESSLIB_LIE_SMALLDEFORMATION_INTEGRATIONPOINTDATAFRACTURE_H_
