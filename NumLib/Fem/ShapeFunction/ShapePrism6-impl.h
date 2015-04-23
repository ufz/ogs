/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

namespace NumLib
{

template <class T_X, class T_N>
void ShapePrism6::computeShapeFunction(const T_X &x, T_N &N)
{
	double L1 = x[0];
	double L2 = x[1];
	double t = x[2];
	N[0] = 0.5 * (1.0 - L1 - L2) * (1.0 - t);
	N[1] = 0.5 * L1 * (1.0 - t);
	N[2] = 0.5 * L2 * (1.0 - t);
	N[3] = 0.5 * (1.0 - L1 - L2) * (1.0 + t);
	N[4] = 0.5 * L1 * (1.0 + t);
	N[5] = 0.5 * L2 * (1.0 + t);
}

template <class T_X, class T_N>
void ShapePrism6::computeGradShapeFunction(const T_X &x, T_N &dN)
{
	double L1 = x[0];
	double L2 = x[1];
	double t = x[2];
	//  dN/dL1
	dN[0] = -0.5 * (1.0 - t);
	dN[1] = 0.5 * (1.0 - t);
	dN[2] = 0.0;
	dN[3] = -0.5 * (1.0 + t);
	dN[4] = 0.5 * (1.0 + t);
	dN[5] = 0.0;
	//  dN/dL2
	dN[6] = -0.5 * (1.0 - t);
	dN[7] = 0.0;
	dN[8] = 0.5 * (1.0 - t);
	dN[9] = -0.5 * (1.0 + t);
	dN[10] = 0.0;
	dN[11] = 0.5 * (1.0 + t);
	//  dN/dt
	dN[12] = -0.5 * (1.0 - L1 - L2);
	dN[13] = -0.5 * L1;
	dN[14] = -0.5 * L2;
	dN[15] = 0.5 * (1.0 - L1 - L2);
	dN[16] = 0.5 * L1;
	dN[17] = 0.5 * L2;
}

}

