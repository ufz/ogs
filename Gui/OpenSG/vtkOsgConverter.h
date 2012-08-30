/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 * \file vtkOsgConverter.h
 *
 * Created on 2011-07-27 by Lars Bilke
 * Derived from class vtkOsgActor from Bjoern Zehner
 */

#ifndef VTKOSGCONVERTER_H
#define VTKOSGCONVERTER_H

#include <OpenSG/OSGChunkMaterial.h>
#include <OpenSG/OSGNode.h>
#include <OpenSG/OSGRefPtr.h>
#include <OpenSG/OSGTextureChunk.h>
#include <OpenSG/OSGTransform.h>

class vtkActor;
class vtkMapper;
class vtkTexture;

/// @brief Converts a vtkActor to an OpenSG-Node.
/// Example usage:
///
/// @code
/// vtkOsgConverter* osgConverter = new vtkOsgConverter(aVtkActor);
/// if(osgConverter->WriteAnActor())
/// {
///   beginEditCP(rootNode);
///   rootNode->addChild(osgConverter->GetOsgNode());
///   endEditCP(rootNode);
/// }
/// @endcode
class vtkOsgConverter
{
public:
	vtkOsgConverter(vtkActor* actor);
	virtual ~vtkOsgConverter();

	bool WriteAnActor();
	void SetVerbose(bool value);
	OSG::NodePtr GetOsgNode();

protected:

private:
	vtkActor* _actor;

	enum {NOT_GIVEN, PER_VERTEX, PER_CELL};
	bool _verbose;

	//For the translation to OpenSG
	OSG::RefPtr<OSG::NodePtr> _osgRoot;
	OSG::RefPtr<OSG::TransformPtr> _osgTransform;

	OSG::TextureChunkPtr CreateTexture(vtkTexture* vtkTexture);
	OSG::ChunkMaterialPtr CreateMaterial(bool lit, bool hasTexCoords);
};

#endif // VTKOSGCONVERTER_H