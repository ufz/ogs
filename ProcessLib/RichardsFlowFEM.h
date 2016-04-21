/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_RICHARDSFLOW_FEM_H_
#define PROCESS_LIB_RICHARDSFLOW_FEM_H_

#include <vector>

#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "NumLib/Function/Interpolation.h"

#include "Parameter.h"
#include "ProcessUtil.h"
#include "Richards_materialproperty.h"
#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
namespace ProcessLib
{

namespace RichardsFlow
{

// TODO now this interface is basically the same for all processes that assemble a
//      FirstOrderImplicitQuasiLinear ODE system.
template <typename GlobalMatrix, typename GlobalVector>
class LocalAssemblerDataInterface
{
public:
    virtual ~LocalAssemblerDataInterface() = default;

    virtual void assemble(double const t, std::vector<double> const& local_x) = 0;

    virtual void addToGlobal(AssemblerLib::LocalToGlobalIndexMap::RowColumnIndices const&,
            GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b) const = 0;
};

template <typename ShapeFunction_,
         typename IntegrationMethod_,
         typename GlobalMatrix,
         typename GlobalVector,
         unsigned GlobalDim>
class LocalAssemblerData : public LocalAssemblerDataInterface<GlobalMatrix, GlobalVector>
{
public:
    using ShapeFunction = ShapeFunction_;
    using ShapeMatricesType = ShapeMatrixPolicyType<ShapeFunction, GlobalDim>;
    using NodalMatrixType = typename ShapeMatricesType::NodalMatrixType;
    using NodalVectorType = typename ShapeMatricesType::NodalVectorType;
    using ShapeMatrices = typename ShapeMatricesType::ShapeMatrices;

    /// The intrinsic_permeability factor is directly integrated into the local
    /// element matrix.
    LocalAssemblerData(MeshLib::Element const& e,
                       std::size_t const local_matrix_size,
                       unsigned const integration_order,
                       Parameter<double, MeshLib::Element const&> const&
						intrinsic_permeability,
						Parameter<double, MeshLib::Element const&> const&
						porosity,
						Parameter<double, MeshLib::Element const&> const&
						viscosity,
						bool const& gravity,
						std::map<std::string,
							std::unique_ptr<MathLib::PiecewiseLinearInterpolation >> const&
						curves)
        : _shape_matrices(
              initShapeMatrices<ShapeFunction, ShapeMatricesType, IntegrationMethod_, GlobalDim>(
                  e, integration_order))
        , _intrinsic_permeability([&intrinsic_permeability, &e]()
          {
              return intrinsic_permeability(e);
          })
		, _porosity([&porosity, &e]()
		  {
			  return porosity(e);
		  })
		, _viscosity([&viscosity, &e]()
		  {
			  return viscosity(e);
		  })
		, _has_gravity(gravity)
		, _curves(curves)
        // TODO narrowing conversion
		, _localM(local_matrix_size, local_matrix_size)
        , _localA(local_matrix_size, local_matrix_size)
        , _localRhs(local_matrix_size)
        , _integration_order(integration_order)
    {
		const MeshLib::CoordinateSystem coordsystem(e);
		const MeshLib::ElementCoordinatesMappingLocal ele_local_coord(e, coordsystem);// wp.getCoords());
		
	}

    void assemble(double const /*t*/, std::vector<double> const& local_x) override
    {
        _localA.setZero();
		_localM.setZero();
        _localRhs.setZero();
		const double rho_w = 1000.;//water density
		const double storage = 0.0;
		double Sw(0.0);//water saturation
		double Pc(0.0);//capillary pressure
		double k_rel = 0.0;  // relative permeability
		double drhow_dp(0.0);
		double dSwdPc(0.0);

		Eigen::MatrixXd mass_mat_coeff = Eigen::MatrixXd::Zero(1, 1);
		Eigen::MatrixXd K_mat_coeff = Eigen::MatrixXd::Zero(1, 1);
		
		MathLib::PiecewiseLinearInterpolation const&  interP_Pc = *_curves.at("curveA");
		MathLib::PiecewiseLinearInterpolation const&  interP_Kr = *_curves.at("curveB");

		//const bool hasGravityEffect = false;

        IntegrationMethod_ integration_method(_integration_order);
        unsigned const n_integration_points = integration_method.getNPoints();//retuen gauss point number

		double P_int_pt = 0.0;
		std::array<double*, 1> const int_pt_array = { &P_int_pt };

		for (std::size_t ip(0); ip < n_integration_points; ip++) {
            auto const& sm = _shape_matrices[ip];
            auto const& wp = integration_method.getWeightedPoint(ip);

			NumLib::shapeFunctionInterpolate(local_x, sm.N, int_pt_array);
			
			Pc = -P_int_pt;
			
			//Sw = getSwbyPc_van(Pc);
			Sw = interP_Pc.getValue(Pc);//read from Pc-S curve
			//dSwdPc = getdSwdPc_van(Pc);
			dSwdPc = interP_Pc.getSlope(Pc);//read from slope of Pc-S curve
			//k_rel = getKrelbySw_van(Sw,0);
			k_rel = interP_Kr.getValue(Sw);//read from S-Kr curve

			mass_mat_coeff(0, 0) = storage * Sw + _porosity() * Sw * drhow_dp - _porosity() * dSwdPc;
			K_mat_coeff(0, 0) = _intrinsic_permeability()*k_rel / _viscosity();

            _localA.noalias() += sm.dNdx.transpose() *
									K_mat_coeff(0,0) * sm.dNdx *
                                 sm.detJ * wp.getWeight();
			_localM.noalias() += sm.N *
				mass_mat_coeff(0, 0) * sm.N.transpose() *
				sm.detJ * wp.getWeight();//Eigen::Map<Eigen::VectorXd>

			if (_has_gravity) {

				Eigen::Vector3d const vec_g(0, 0, -9.81);
				// since no primary vairable involved
				// directly assemble to the Right-Hand-Side
				// F += dNp^T * K * gz
				_localRhs.noalias() += sm.dNdx.transpose() * K_mat_coeff(0, 0) * rho_w * vec_g * sm.detJ * wp.getWeight();
			} // end of if hasGravityEffect

        }
    }

    void addToGlobal(AssemblerLib::LocalToGlobalIndexMap::RowColumnIndices const& indices,
        GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b)
        const override
    {
		M.add(indices, _localM);
        K.add(indices, _localA);
        b.add(indices.rows, _localRhs);
    }

private:
    std::vector<ShapeMatrices> _shape_matrices;
    std::function<double(void)> _intrinsic_permeability;
	std::function<double(void)> _porosity;
	std::function<double(void)> _viscosity;
	std::map<std::string,
		std::unique_ptr<MathLib::PiecewiseLinearInterpolation >> const&
		_curves;
	bool const _has_gravity;

    NodalMatrixType _localA;
	NodalMatrixType _localM;
    NodalVectorType _localRhs;

    unsigned const _integration_order;


};


}   // namespace RichardsFlow
}   // namespace ProcessLib

#endif  // PROCESS_LIB_RICHARDSFLOW_FEM_H_
