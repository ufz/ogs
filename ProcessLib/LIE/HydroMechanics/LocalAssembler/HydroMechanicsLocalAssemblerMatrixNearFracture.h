/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_LIE_HYDROMECHANICS_HYDROMECHANICSLOCALASSEMBLER_MATRIX_NEAR_FRACTUER_H_
#define PROCESSLIB_LIE_HYDROMECHANICS_HYDROMECHANICSLOCALASSEMBLER_MATRIX_NEAR_FRACTUER_H_

#include <vector>

#include "HydroMechanicsLocalAssemblerMatrix.h"

namespace ProcessLib
{
namespace LIE
{
namespace HydroMechanics
{
template <typename ShapeFunctionDisplacement,
          typename ShapeFunctionPressure,
          typename IntegrationMethod,
          unsigned GlobalDim>
class HydroMechanicsLocalAssemblerMatrixNearFracture
    : public HydroMechanicsLocalAssemblerMatrix<ShapeFunctionDisplacement,
                                                ShapeFunctionPressure,
                                                IntegrationMethod,
                                                GlobalDim>
{
    using Base = HydroMechanicsLocalAssemblerMatrix<ShapeFunctionDisplacement,
                                                    ShapeFunctionPressure,
                                                    IntegrationMethod,
                                                    GlobalDim>;

public:
    HydroMechanicsLocalAssemblerMatrixNearFracture(HydroMechanicsLocalAssemblerMatrixNearFracture const&) =
        delete;
    HydroMechanicsLocalAssemblerMatrixNearFracture(HydroMechanicsLocalAssemblerMatrixNearFracture&&) = delete;

    HydroMechanicsLocalAssemblerMatrixNearFracture(
        MeshLib::Element const& e,
        std::size_t const n_variables,
        std::size_t const local_matrix_size,
        std::vector<unsigned> const& dofIndex_to_localIndex,
        bool const is_axially_symmetric,
        unsigned const integration_order,
        HydroMechanicsProcessData<GlobalDim>& process_data);

private:
    void assembleWithJacobianConcrete(
        double const t,
        Eigen::VectorXd const& local_u,
        Eigen::VectorXd const& local_udot,
        Eigen::VectorXd& local_b,
        Eigen::MatrixXd& local_J) override;

    void computeSecondaryVariableConcreteWithVector(
        double const t, Eigen::VectorXd const& local_x) override;

    using typename Base::ShapeMatricesTypeDisplacement;
    using typename Base::BMatricesType;
    using Base::_process_data;
    using Base::_ip_data;
    using Base::_element;
    using Base::pressure_index;
    using Base::pressure_size;
    using Base::displacement_index;
    using Base::displacement_size;
    using Base::kelvin_vector_size;

    static const int displacement_jump_index =
        displacement_index + displacement_size;
};

}  // namespace HydroMechanics
}  // namespace LIE
}  // namespace ProcessLib

#include "HydroMechanicsLocalAssemblerMatrixNearFracture-impl.h"

#endif // PROCESSLIB_LIE_HYDROMECHANICS_HYDROMECHANICSLOCALASSEMBLER_MATRIX_H_
