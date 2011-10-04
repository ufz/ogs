/*
 * LinearInterpolation.h
 *
 *  Created on: Sep 7, 2010
 *      Author: TF
 */

#ifndef LINEARINTERPOLATION_H_
#define LINEARINTERPOLATION_H_

#include <vector>
#include <limits>

namespace MathLib {

class LinearInterpolation {
public:
	LinearInterpolation(const std::vector<double>& supporting_points, const std::vector<double>& values_at_supp_pnts);
	LinearInterpolation(const std::vector<double>& supporting_points, const std::vector<double>& values_at_supp_pnts, const std::vector<double>& points_to_interpolate, std::vector<double>& values_at_interpol_pnts);
	virtual ~LinearInterpolation();

	double getValue ( double pnt_to_interpolate );

private:
	const std::vector<double>& _supporting_points;
	const std::vector<double>& _values_at_supp_pnts;
};

} // end namespace MathLib

#endif /* LINEARINTERPOLATION_H_ */
