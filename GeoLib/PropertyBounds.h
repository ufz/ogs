/**
 * \file PropertyBounds.h
 * 18/01/2010 KR Initial implementation
 */

#ifndef PROPERTYBOUNDS_H
#define PROPERTYBOUNDS_H

#include <string>

class PropertyBounds
{
public:
	PropertyBounds(std::string pname, double minVal, double maxVal ) :  _name(pname), _minVal(minVal), _maxVal(maxVal) {}

	std::string getName() const { return _name; }
	double getMin() const { return _minVal; }
	double getMax() const { return _maxVal; }

	void setMin(double val) { _minVal = val; }
	void setMax(double val) { _maxVal = val; }

private:
	std::string _name;
	double _minVal;
	double _maxVal;
};

#endif //PROPERTYBOUNDS_H
