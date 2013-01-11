/**
 * \file
 * \author Karsten Rink
 * \date   2012-11-26
 * \brief  Implementation of the zLibDataCompressor class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "zLibDataCompressor.h"
#include <cstddef>
#include <iostream>
#include "zlib/zlib.h"


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
		std::cout << "Zlib error while compressing data." << std::endl;
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
		std::cout << "Zlib error while uncompressing data." << std::endl;
		return 0;
	}

	// Make sure the output size matched that expected.
	if(decSize != uncompressedSize)
	{
		std::cout << "Decompression produced incorrect size. Expected "
			      << uncompressedSize << " and got " << decSize << std::endl;
		return 0;
	}

	return decSize;
}
