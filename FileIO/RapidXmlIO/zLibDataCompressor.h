/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 * \file zLibDataCompressor.h
 *
 * Created on 2012-11-26 by Karsten Rink
 * Based on the vtkZLibDataCompressor-class in VTK 5.6
 */

#ifndef ZLIBDATACOMPRESSOR_H
#define ZLIBDATACOMPRESSOR_H



class zLibDataCompressor
{
public:
	unsigned long GetMaximumCompressionSpace(unsigned long size);

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
private:

};

#endif //ZLIBDATACOMPRESSOR_H

