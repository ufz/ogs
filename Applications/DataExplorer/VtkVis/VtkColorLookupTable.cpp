/**
 * \file
 * \author Karsten Rink
 * \date   2010-04-23
 * \brief  Implementation of the VtkColorLookupTable class.
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "VtkColorLookupTable.h"

#include <cmath>
#include <sstream>

// ThirdParty/logog
#include "logog/include/logog.hpp"

#include <Color.h>
#include <vtkObjectFactory.h>

vtkStandardNewMacro(VtkColorLookupTable);

VtkColorLookupTable::VtkColorLookupTable()
	: _type(VtkColorLookupTable::LUTType::LINEAR)
{
}

VtkColorLookupTable::~VtkColorLookupTable()
{
	for (std::map<double, unsigned char*>::const_iterator it = _dict.begin(); it != _dict.end();
	     ++it)
		delete it->second;
}

unsigned char VtkColorLookupTable::linInterpolation(unsigned char a, unsigned char b,
                                                    double p) const
{
	return static_cast<unsigned char>(a * (1 - p) + b * p);
}

unsigned char VtkColorLookupTable::expInterpolation(unsigned char a,
                                                    unsigned char b,
                                                    double gamma,
                                                    double p) const
{
	assert (gamma > 0 && gamma < 4);
	return static_cast<unsigned char>((b - a) * pow(p,gamma) + a);
}

void VtkColorLookupTable::Build()
{
	double range[2];
	this->GetTableRange(range);
	const double interval = range[1]-range[0];
	this->SetNumberOfTableValues(ceil(interval)+1);
//	const vtkIdType nColours = this->GetNumberOfTableValues();
	if (!_dict.empty())
	{
		// make sure that color map starts with the first color in the dictionary
		unsigned char startcolor[4] = { 0, 0 , 0 , 0 };
		std::pair<size_t, unsigned char*> lastValue(0, startcolor);
		size_t nextIndex(0);

		for (std::map<double, unsigned char*>::const_iterator it = _dict.begin(); it != _dict.end(); ++it)
		{
			double val = (it->first < range[0]) ? range[0] : ((it->first > range[1]) ? range[1] : it->first);
			nextIndex = static_cast<size_t>( std::floor(val-range[0]) );

			this->SetTableValue(nextIndex, it->second);

			if ( nextIndex - lastValue.first > 0 )
				for (size_t i = lastValue.first + 1; i < nextIndex; i++)
				{
					unsigned char int_rgba[4];
					double pos = (i - lastValue.first) / (static_cast<double>(nextIndex - lastValue.first));

					if (_type == VtkColorLookupTable::LUTType::LINEAR)
						for (size_t j = 0; j < 4; j++)
							int_rgba[j] = linInterpolation( (lastValue.second)[j], (it->second)[j], pos);
					else if (_type == VtkColorLookupTable::LUTType::EXPONENTIAL)
						for (size_t j = 0; j < 4; j++)
							int_rgba[j] = expInterpolation((lastValue.second)[j], (it->second)[j], 0.2, pos);
					else // no interpolation
						for (size_t j = 0; j < 4; j++)
							int_rgba[j] = (lastValue.second)[j];

					this->SetTableValue(i, int_rgba);
				}

			lastValue.first = nextIndex;
			lastValue.second = it->second;
		}
	}
	else
		vtkLookupTable::Build();
}

void VtkColorLookupTable::writeToFile(const std::string &filename)
{
	std::stringstream strout;
	strout << "Writing color table to " << filename << " ... ";
	std::ofstream out( filename.c_str(), std::ios::out );

	size_t nColors = this->GetNumberOfTableValues();
	for (size_t i = 0; i < nColors; i++)
	{
		unsigned char rgba[4];
		this->GetTableValue(i, rgba);
		out << i << "\t" << rgba[0] << "\t" << rgba[1] << "\t" << rgba[2] << "\n";
	}

	strout << " done." << std::endl;
	INFO("%s", strout.str().c_str());
	out.close();
}

void VtkColorLookupTable::SetTableValue(vtkIdType idx, unsigned char rgba[4])
{
	double value[4];

	for (unsigned i=0; i<4; ++i)
		value[i] = rgba[i]/255.0;
	vtkLookupTable::SetTableValue(idx, value);
}

void VtkColorLookupTable::GetTableValue(vtkIdType idx, unsigned char rgba[4])
{
	double value[4];
	vtkLookupTable::GetTableValue(idx, value);

	for (unsigned i=0; i<4; ++i)
		rgba[i] = value[i]*255.0;
}

void VtkColorLookupTable::setColor(double pos, unsigned char rgba[4])
{
	unsigned char* dict_rgba = new unsigned char[4];
	for (size_t i = 0; i < 4; i++)
		dict_rgba[i] = rgba[i];
	_dict.insert( std::pair<double, unsigned char*>(pos, dict_rgba) );
}

void VtkColorLookupTable::getColor(vtkIdType indx, unsigned char rgba[4]) const
{
	indx =
	        ((indx < this->TableRange[0])
				? static_cast<vtkIdType>(this->TableRange[0])
				: (indx >=this->TableRange[1] ? static_cast<vtkIdType>(this->TableRange[1]) - 1 : indx));
	indx =
	        static_cast<size_t>( std::floor( (indx - this->TableRange[0]) *
	                                    (this->NumberOfColors / (this->TableRange[1] - this->TableRange[0])) ) );

	unsigned char* _rgba;
	_rgba = this->Table->GetPointer(indx * 4);
	for (size_t i = 0; i < 4; i++)
		rgba[i] = _rgba[i];
}
