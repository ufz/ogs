/**
 * \file    ElementErrorCode.h
 * \author  Karsten Rink
 * \date    2014-02-21
 * \brief   Definition of ElementErrorCodes.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef ELEMENTERRORCODES_H
#define ELEMENTERRORCODES_H

#include <bitset>
#include <string>


/// Possible error flags for mesh elements
enum class ElementErrorFlag 
{
	ZeroVolume,
	NonCoplanar,
	NonConvex,
	NodeOrder,
	//... add other error flags here
	MaxValue // this needs to be last to set the bitset size correctly!
};

/// Collects error flags for mesh elements
class ElementErrorCode
{
public:
	ElementErrorCode() {}
	ElementErrorCode(std::bitset< static_cast<std::size_t>(ElementErrorFlag::MaxValue) > error_flags) 
		: _errors(error_flags) {}

	~ElementErrorCode() {}

	/// Get value for a specific flag
	bool get(ElementErrorFlag e) const { return _errors[static_cast<std::size_t>(e)]; }
	/// Set a specific flag
	void set(ElementErrorFlag e) { _errors[static_cast<std::size_t>(e)] = true; }
	/// Reset a specific flag
	void reset(ElementErrorFlag e) { _errors[static_cast<std::size_t>(e)] = false; }
	
	/// Returns true if ANY flag is set, false otherwise
	bool any() const { return _errors.any(); }
	/// Returns true if ALL flag is set, false otherwise
	bool all() const { return _errors.all(); }
	/// Returns true if NO flags are set, false otherwise
	bool none() const { return _errors.none(); }

	/// Returns the size of the error flag set	
	std::size_t size() const { return _errors.size(); }

	/// Returns the underlying bitset
	const std::bitset< static_cast<std::size_t>(ElementErrorFlag::MaxValue) >& bitset() const { return _errors; }

	//inline bool operator[](const ElementErrorFlag e) { return _errors[static_cast<std::size_t>(e)]; }
	//inline const bool operator[](const ElementErrorFlag e) const { return _errors[static_cast<std::size_t>(e)]; }

	inline ElementErrorCode operator|=(ElementErrorCode e) { return _errors |= e.bitset(); }

	/// Returns a string output for a specific error flag
	static std::string toString(const ElementErrorFlag e)
	{
		if (e == ElementErrorFlag::ZeroVolume)
			return "zero volume";
		else if (e == ElementErrorFlag::NonCoplanar)
			return "non coplanar nodes";
		else if (e == ElementErrorFlag::NonConvex)
			return "non-convex geometry";
		else if (e == ElementErrorFlag::NodeOrder)
			return "wrong node order";
		return "nonspecified error"; 
	}

private:
	/// The bit set collecting the error flags
    std::bitset< static_cast<std::size_t>(ElementErrorFlag::MaxValue) > _errors;
};


#endif //ELEMENTERRORCODES_H
