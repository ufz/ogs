/**
 * \file
 * \author Karsten Rink
 * \date   2012-11-26
 * \brief  Definition of the zLibDataCompressor class.
 *         Based on the vtkZLibDataCompressor-class in VTK 5.6
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Based on the vtkZLibDataCompressor-class in VTK 5.6
 */

#ifndef ZLIBDATACOMPRESSOR_H
#define ZLIBDATACOMPRESSOR_H

class zLibDataCompressor
{
public:
	// Compression method required by vtkDataCompressor.
	static unsigned long CompressBuffer(const unsigned char* uncompressedData,
	                                    unsigned long uncompressedSize,
	                                    unsigned char* compressedData,
	                                    unsigned long compressionSpace);

	// Decompression method required by vtkDataCompressor.
	static unsigned long UncompressBuffer(const unsigned char* compressedData,
	                                      unsigned long compressedSize,
	                                      unsigned char* uncompressedData,
	                                      unsigned long uncompressedSize);
};

#endif //ZLIBDATACOMPRESSOR_H

