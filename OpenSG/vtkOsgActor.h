/////////////////////////////////////////////////////////////////////////////////////////
//                       By Bjoern Zehner 2006/2010                                    //
//			UFZ-Helmholtz Center for Environmental Research Leipzig                    //
//			      Permoser Str. 15, D-04318 Leizig, Germany                            //
/////////////////////////////////////////////////////////////////////////////////////////
// Class for translating VTK Objects into OpenSG representations, derives from vtkActor//
// See Project/file VtkOsgActorTestXXX for an example for the use                      //

#pragma once
#include "vtkOpenGLActor.h"

#include "vtkObjectFactory.h"
#include "vtkTimeStamp.h"

#include <vtkActor.h>
#include <vtkPolyDataMapper.h>
#include <vtkMapper.h>
#include <vtkPolyData.h>
#include <vtkDataSet.h>
#include <vtkPointData.h>
#include <vtkProperty.h>
#include <vtkCellData.h>
#include <vtkCellArray.h>
#include <vtkTexture.h>
#include <vtkImageData.h>

#include <OpenSG/OSGRefPtr.h>
#include <OpenSG/OSGNode.h>
#include <OpenSG/OSGTransform.h>
#include <OpenSG/OSGMatrix.h>
#include <OpenSG/OSGSimpleGeometry.h>
#include <OpenSG/OSGSimpleMaterial.h>
#include <OpenSG/OSGChunkMaterial.h>
#include <OpenSG/OSGMaterialChunk.h>
#include <OpenSG/OSGTextureChunk.h>
#include <OpenSG/OSGGeoFunctions.h>
#include <OpenSG/OSGGroup.h>
#include <OpenSG/OSGTwoSidedLightingChunk.h>
#include <OpenSG/OSGPolygonChunk.h>
#include <OpenSG/OSGGeoFunctions.h>
#include <OpenSG/OSGTransform.h>
#include <OpenSG/OSGImage.h>

OSG_USING_NAMESPACE

class vtkOsgActor : public vtkOpenGLActor
{
protected:
	vtkOsgActor(void);
	virtual ~vtkOsgActor(void);
	void InitOpenSG(void);

public: //VTK related
	vtkTypeMacro(vtkOsgActor,vtkOpenGLActor);
	static vtkOsgActor* New();
	void PrintSelf(ostream& os, vtkIndent indent);
	void Render(vtkRenderer *ren, vtkMapper *mapper);

public: //The added value - OpenSG related stuff
	void UpdateOsg();
	void ClearOsg();
	void SetVerbose(bool value);
	void SetTexture(vtkTexture *vtkTex);
	NodePtr GetOsgRoot();

private:
	vtkDataArray			*m_pvtkNormals;
	vtkDataArray			*m_pvtkTexCoords;
	vtkUnsignedCharArray	*m_pvtkColors;
	vtkTexture				*m_pvtkTexture;
	bool					m_bTextureHasChanged;

	enum {NOT_GIVEN, PER_VERTEX, PER_CELL};
	int						m_iColorType;
	int						m_iNormalType;
	bool					m_bVerbose;

	int						m_iNumPoints;
	int						m_iNumNormals;
	int						m_iNumColors;
	int						m_iNumGLPoints;
	int						m_iNumGLLineStrips;
	int						m_iNumGLPolygons;
	int						m_iNumGLTriStrips;
	int						m_iNumGLPrimitives;

	//For the translation to OpenSG
	RefPtr<NodePtr> m_posgRoot;
	RefPtr<TransformPtr> m_posgTransform;
	RefPtr<NodePtr> m_posgGeomNode;
	RefPtr<GeometryPtr> m_posgGeometry;
	RefPtr<ChunkMaterialPtr> m_posgMaterial;
	RefPtr<MaterialChunkPtr> m_posgMaterialChunk;
	RefPtr<TextureChunkPtr> m_posgTextureChunk;
	RefPtr<PolygonChunkPtr> m_posgPolygonChunk;
	RefPtr<ImagePtr> m_posgImage;

	RefPtr<GeoPTypesPtr> m_posgTypes;
	RefPtr<GeoPLengthsPtr> m_posgLengths;
	RefPtr<GeoIndicesUI32Ptr> m_posgIndices;
	RefPtr<GeoPositions3fPtr> m_posgPoints;
	RefPtr<GeoColors3fPtr> m_posgColors;
	RefPtr<GeoNormals3fPtr> m_posgNormals;
	RefPtr<GeoTexCoords2dPtr> m_posgTexCoords;

private:
  vtkOsgActor(const vtkOsgActor&);  // Not implemented.
  void operator=(const vtkActor&);  // Not implemented.

private:
	//for the translation to OpenSG
	void LookForNormals();
	void LookForColors();
	void LookForTexCoords();
	void LookForArraySizes();

	void CreateTexture();
	ChunkMaterialPtr CreateMaterial();

	//Can use OpenSG simple indexed geometry
	NodePtr ProcessGeometryNormalsAndColorsPerVertex();

	//Can't use indexing and so requires a lot of storage space
	NodePtr ProcessGeometryNonIndexedCopyAttributes(int gl_primitive_type);

	NodePtr GetNodePtr();
};
