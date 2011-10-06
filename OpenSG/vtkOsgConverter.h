/**
 * \file vtkOsgConverter.h
 * 27/07/2011 LB Initial implementation
 * Derived from class vtkOsgActor from Bjoern Zehner
 */

#ifndef VTKOSGCONVERTER_H
#define VTKOSGCONVERTER_H

#include <OpenSG/OSGRefPtr.h>
#include <OpenSG/OSGNode.h>
#include <OpenSG/OSGTransform.h>
#include <OpenSG/OSGTextureChunk.h>
#include <OpenSG/OSGChunkMaterial.h>

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
  vtkMapper* _mapper;

  enum {NOT_GIVEN, PER_VERTEX, PER_CELL};
  bool _verbose;


  //For the translation to OpenSG
  OSG::RefPtr<OSG::NodePtr> _osgRoot;
  OSG::RefPtr<OSG::TransformPtr> _osgTransform;

  OSG::TextureChunkPtr CreateTexture(vtkTexture* vtkTexture);
  OSG::ChunkMaterialPtr CreateMaterial(bool lit, bool hasTexCoords);
};

#endif // VTKOSGCONVERTER_H