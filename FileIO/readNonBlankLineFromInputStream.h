/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 * \file readNonBlankLineFromInputStream.h
 *
 *  Created on 2011-04-19 by Thomas Fischer
 */

#ifndef READNONBLANKLINEFROMINPUTSTREAM_H_
#define READNONBLANKLINEFROMINPUTSTREAM_H_

#include <istream>
#include <string>

/**
 * read a non blank line from given input stream
 * @param in the input stream
 * @return read line into a string
 */
std::string readNonBlankLineFromInputStream(std::istream & in);

#endif /* READNONBLANKLINEFROMINPUTSTREAM_H_ */
