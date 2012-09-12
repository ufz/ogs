/*
 * GMSHLine.h
 *
 *  Created on: Mar 22, 2012
 *      Author: fischeth
 */

#ifndef GMSHLINE_H_
#define GMSHLINE_H_

#include <ostream>

namespace FileIO {

class GMSHLine {
public:
	GMSHLine(std::size_t start_point_id, std::size_t end_point_id);
	virtual ~GMSHLine();
	void write(std::ostream &os, std::size_t id) const;
	void resetLineData(std::size_t start_point_id, std::size_t end_point_id);

private:
	std::size_t _start_pnt_id;
	std::size_t _end_pnt_id;
};

}

#endif /* GMSHLINE_H_ */
