/**
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 *
 * \file zLibDataCompressor.cpp
 *
 * Created on 2012-11-26 by Karsten Rink
 * Based on the vtkZLibDataCompressor-class in VTK 5.6
 */

#include <cstddef>
#include <iostream>

// ThirdParty/logog
#include "logog/include/logog.hpp"

#include "zlib/zlib.h"

#include "zLibDataCompressor.h"

unsigned long zLibDataCompressor::CompressBuffer(const unsigned char* uncompressedData,
                                                 unsigned long uncompressedSize,
                                                 unsigned char* compressedData,
                                                 unsigned long compressionSpace)
{
	int CompressionLevel = Z_DEFAULT_COMPRESSION;
	unsigned long compressedSize = compressionSpace;
	Bytef* cd = reinterpret_cast<Bytef*>(compressedData);
	const Bytef* ud = reinterpret_cast<const Bytef*>(uncompressedData);

	// Call zlib's compress function.
	if(compress2(cd, &compressedSize, ud, uncompressedSize, CompressionLevel) != Z_OK)
	{
		ERR("zLibDataCompressor::CompressBuffer(): Zlib error while compressing data.");
		return 0;
	}

	return compressedSize;
}

unsigned long zLibDataCompressor::UncompressBuffer(const unsigned char* compressedData,
                                                   unsigned long compressedSize,
                                                   unsigned char* uncompressedData,
                                                   unsigned long uncompressedSize)
{
	unsigned long decSize = uncompressedSize;
	Bytef* ud = reinterpret_cast<Bytef*>(uncompressedData);
	const Bytef* cd = reinterpret_cast<const Bytef*>(compressedData);

	// Call zlib's uncompress function.
	if(uncompress(ud, &decSize, cd, compressedSize) != Z_OK)
	{
		ERR("zLibDataCompressor::CompressBuffer(): Zlib error while uncompressing data.");
		return 0;
	}

	// Make sure the output size matched that expected.
	if(decSize != uncompressedSize)
	{
		WARN("Decompression produced incorrect size. Expected %d and got %d.",
		     uncompressedSize, decSize);
		return 0;
	}

	return decSize;
}
