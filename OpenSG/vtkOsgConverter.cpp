/**
 * \file vtkOsgConverter.cpp
 * 27/07/2011 LB Initial implementation
 *
 * Implementation of vtkOsgConverter class
 */

// ** INCLUDES **
#include "vtkOsgConverter.h"

#include <vtkActor.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkCompositeDataGeometryFilter.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkDataSetMapper.h>
#include <vtkGeometryFilter.h>
#include <vtkImageData.h>
#include <vtkMapper.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkSmartPointer.h>
#include <vtkTexture.h>
#include <vtkUnsignedCharArray.h>

#include <OpenSG/OSGGeoFunctions.h>
#include <OpenSG/OSGGroup.h>
#include <OpenSG/OSGImage.h>
#include <OpenSG/OSGLineChunk.h>
#include <OpenSG/OSGMaterialChunk.h>
#include <OpenSG/OSGMatrix.h>
#include <OpenSG/OSGPointChunk.h>
#include <OpenSG/OSGPolygonChunk.h>
#include <OpenSG/OSGSimpleGeometry.h>
#include <OpenSG/OSGTwoSidedLightingChunk.h>

OSG_USING_NAMESPACE

vtkOsgConverter::vtkOsgConverter(vtkActor* actor) :
	_actor(actor),
	_verbose(true),
	_osgRoot(NullFC),
	_osgTransform(NullFC)
{
	TransformPtr tptr;
	_osgRoot = makeCoredNode<osg::Transform>(&tptr);
	_osgTransform = tptr;
	_mapper = _actor->GetMapper();
}

vtkOsgConverter::~vtkOsgConverter(void)
{
	_osgRoot = NullFC;
}

bool vtkOsgConverter::WriteAnActor()
{
	// see if the actor has a mapper. it could be an assembly
	if (_actor->GetMapper() == NULL)
		return false;
	// dont export when not visible
	if (_actor->GetVisibility() == 0)
		return false;

	vtkDataObject* inputDO = _actor->GetMapper()->GetInputDataObject(0, 0);
	if (inputDO == NULL)
		return false;

	// Get PolyData. Convert if necessary becasue we only want polydata
	vtkSmartPointer<vtkPolyData> pd;
	if(inputDO->IsA("vtkCompositeDataSet"))
	{
		vtkCompositeDataGeometryFilter* gf = vtkCompositeDataGeometryFilter::New();
		gf->SetInput(inputDO);
		gf->Update();
		pd = gf->GetOutput();
		gf->Delete();
	}
	else if(inputDO->GetDataObjectType() != VTK_POLY_DATA)
	{
		vtkGeometryFilter* gf = vtkGeometryFilter::New();
		gf->SetInput(inputDO);
		gf->Update();
		pd = gf->GetOutput();
		gf->Delete();
	}
	else
		pd = static_cast<vtkPolyData*>(inputDO);

	// Copy mapper to a new one
	vtkPolyDataMapper* pm = vtkPolyDataMapper::New();
	pm->SetInput(pd);
	pm->SetScalarRange(_actor->GetMapper()->GetScalarRange());
	pm->SetScalarVisibility(_actor->GetMapper()->GetScalarVisibility());
	pm->SetLookupTable(_actor->GetMapper()->GetLookupTable());
	pm->SetScalarMode(_actor->GetMapper()->GetScalarMode());

	if(pm->GetScalarMode() == VTK_SCALAR_MODE_USE_POINT_FIELD_DATA ||
	   pm->GetScalarMode() == VTK_SCALAR_MODE_USE_CELL_FIELD_DATA )
	{
		if(_actor->GetMapper()->GetArrayAccessMode() == VTK_GET_ARRAY_BY_ID )
			pm->ColorByArrayComponent(_actor->GetMapper()->GetArrayId(),
			                          _actor->GetMapper()->GetArrayComponent());
		else
			pm->ColorByArrayComponent(_actor->GetMapper()->GetArrayName(),
			                          _actor->GetMapper()->GetArrayComponent());
	}

	_mapper = pm;
	// vtkPoints* points = pd->GetPoints();
	vtkPointData* pntData = pd->GetPointData();
	bool hasTexCoords = false;
	vtkUnsignedCharArray* vtkColors  = pm->MapScalars(1.0);

	// ARRAY SIZES
	vtkIdType m_iNumPoints = pd->GetNumberOfPoints();
	if (m_iNumPoints == 0)
		return false;
	vtkIdType m_iNumGLPoints = pd->GetVerts()->GetNumberOfCells();
	vtkIdType m_iNumGLLineStrips = pd->GetLines()->GetNumberOfCells();
	vtkIdType m_iNumGLPolygons = pd->GetPolys()->GetNumberOfCells();
	vtkIdType m_iNumGLTriStrips = pd->GetStrips()->GetNumberOfCells();
	vtkIdType m_iNumGLPrimitives = m_iNumGLPoints + m_iNumGLLineStrips + m_iNumGLPolygons +
	                               m_iNumGLTriStrips;
	bool lit = !(m_iNumGLPolygons == 0 && m_iNumGLTriStrips == 0);

	if (_verbose)
	{
		std::cout << "Array sizes:" << std::endl;
		std::cout << "  number of vertices: " << m_iNumPoints << std::endl;
		std::cout << "  number of GL_POINTS: " << m_iNumGLPoints << std::endl;
		std::cout << "  number of GL_LINE_STRIPS: " << m_iNumGLLineStrips << std::endl;
		std::cout << "  number of GL_POLYGON's: " << m_iNumGLPolygons << std::endl;
		std::cout << "  number of GL_TRIANGLE_STRIPS: " << m_iNumGLTriStrips << std::endl;
		std::cout << "  number of primitives: " << m_iNumGLPrimitives << std::endl;
	}

	_mapper->Update();

	// NORMALS
	vtkDataArray* vtkNormals = NULL;
	int m_iNormalType = NOT_GIVEN;
	if (_actor->GetProperty()->GetInterpolation() == VTK_FLAT)
	{
		vtkNormals = pd->GetCellData()->GetNormals();
		if (vtkNormals != NULL)
			m_iNormalType = PER_CELL;
	}
	else
	{
		vtkNormals = pntData->GetNormals();
		if (vtkNormals != NULL)
			m_iNormalType = PER_VERTEX;
	}
	if (_verbose)
	{
		std::cout << "Normals:" << std::endl;
		if (m_iNormalType != NOT_GIVEN)
		{
			std::cout << "  number of normals: " << vtkNormals->GetNumberOfTuples() <<
			std::endl;
			std::cout << "  normals are given: ";
			std::cout <<
			((m_iNormalType == PER_VERTEX) ? "per vertex" : "per cell") << std::endl;
		}
		else
			std::cout << "  no normals are given" << std::endl;
	}

	// COLORS
	int m_iColorType = NOT_GIVEN;
	if(pm->GetScalarVisibility())
	{
		int iScalarMode = pm->GetScalarMode();
		if(vtkColors == NULL)
		{
			m_iColorType = NOT_GIVEN;
			std::cout << "WARNING: MapScalars(1.0) did not return array!" << std::endl;
		}
		else if(iScalarMode == VTK_SCALAR_MODE_USE_CELL_DATA)
			m_iColorType = PER_CELL;
		else if(iScalarMode == VTK_SCALAR_MODE_USE_POINT_DATA)
			m_iColorType = PER_VERTEX;
		else if(iScalarMode == VTK_SCALAR_MODE_USE_CELL_FIELD_DATA)
		{
			std::cout <<
			"WARNING TO BE REMOVED: Can not process colours with scalar mode using cell field data!"
			          << std::endl;
			m_iColorType = PER_CELL;
		}
		else if(iScalarMode == VTK_SCALAR_MODE_USE_POINT_FIELD_DATA)
		{
			std::cout <<
			"WARNING TO BE REMOVED: Can not process colours with scalar mode using point field data!"
			          << std::endl;
			m_iColorType = PER_VERTEX;
		}
		else if(iScalarMode == VTK_SCALAR_MODE_DEFAULT)
		{
			//Bummer, we do not know what it is. may be we can make a guess
			int numColors = vtkColors->GetNumberOfTuples();
			if (numColors == 0)
			{
				m_iColorType = NOT_GIVEN;
				std::cout << "WARNING: No colors found!" << std::endl;
			}
			else if (numColors == m_iNumPoints)
				m_iColorType = PER_VERTEX;
			else if (numColors == m_iNumGLPrimitives)
				m_iColorType = PER_CELL;
			else
			{
				m_iColorType = NOT_GIVEN;
				std::cout <<
				"WARNING: Number of colors do not match number of points / cells!"
				          << std::endl;
			}
		}
	}
	if (_verbose)
	{
		std::cout << "Colors:" << std::endl;
		if (m_iColorType != NOT_GIVEN)
		{
			std::cout << "  number of colors: " << vtkColors->GetNumberOfTuples() <<
			std::endl;
			std::cout << "  colors are given: " <<
			((m_iColorType == PER_VERTEX) ? "per vertex" : "per cell") << std::endl;
		}
		else
			std::cout << "  no colors are given" << std::endl;
	}

	// TEXCOORDS
	vtkDataArray* vtkTexCoords = pntData->GetTCoords();
	if (_verbose)
	{
		std::cout << "Tex-coords:" << std::endl;
		if (vtkTexCoords)
		{
			std::cout << "  Number of tex-coords: " <<
			vtkTexCoords->GetNumberOfTuples() << std::endl;
			hasTexCoords = true;
		}
		else
			std::cout << "  No tex-coords where given" << std::endl;
	}

	// TRANSFORMATION
	double scaling[3];
	double translation[3];
	// double rotation[3];

	_actor->GetPosition(translation);
	_actor->GetScale(scaling);
	//_actor->GetRotation(rotation[0], rotation[1], rotation[2]);

	if (_verbose)
		std::cout << "set scaling: " << scaling[0] << " " << scaling[1] << " " <<
		scaling[2] << std::endl;

	osg::Matrix m;
	m.setIdentity();
	m.setTranslate(translation[0], translation[1], translation[2]);
	m.setScale(scaling[0], scaling[1], scaling[2]);
	// TODO QUATERNION m.setRotate(rotation[0], rotation[1], rotation[2])
	beginEditCP(_osgTransform);
	_osgTransform->setMatrix(m);
	endEditCP(_osgTransform);

	_mapper->Update();

	// Get the converted OpenSG node
	NodePtr osgGeomNode = Node::create();
	GeometryPtr osgGeometry = Geometry::create();
	beginEditCP(osgGeomNode);
	osgGeomNode->setCore(osgGeometry);
	endEditCP(osgGeomNode);

	bool osgConversionSuccess = false;

	GeoPTypesPtr osgTypes = GeoPTypesUI8::create();
	GeoPLengthsPtr osgLengths = GeoPLengthsUI32::create();
	GeoIndicesUI32Ptr osgIndices = GeoIndicesUI32::create();
	GeoPositions3fPtr osgPoints = GeoPositions3f::create();
	GeoNormals3fPtr osgNormals = GeoNormals3f::create();
	GeoColors3fPtr osgColors = GeoColors3f::create();
	GeoTexCoords2dPtr osgTexCoords = GeoTexCoords2d::create();

	//Rendering with OpenSG simple indexed geometry
	if (((m_iNormalType == PER_VERTEX) || (m_iNormalType == NOT_GIVEN))  &&
	    ((m_iColorType == PER_VERTEX) || (m_iColorType == NOT_GIVEN)))
	{
		if (_verbose)
			std::cout << "Start ProcessGeometryNormalsAndColorsPerVertex()" <<
			std::endl;

		//getting the vertices:
		beginEditCP(osgPoints);
		{
			for (int i = 0; i < m_iNumPoints; i++)
			{
				double* aVertex = pd->GetPoint(i);
				osgPoints->addValue(Vec3f(aVertex[0], aVertex[1], aVertex[2]));
			}
		} endEditCP(osgPoints);

		//possibly getting the normals
		if (m_iNormalType == PER_VERTEX)
		{
			vtkIdType iNumNormals = vtkNormals->GetNumberOfTuples();
			beginEditCP(osgNormals);
			{
				double* aNormal;
				for (int i = 0; i < iNumNormals; i++)
				{
					aNormal = vtkNormals->GetTuple(i);
					osgNormals->addValue(Vec3f(aNormal[0], aNormal[1],
					                           aNormal[2]));
				}
			} endEditCP(osgNormals);
			if (iNumNormals != m_iNumPoints)
			{
				std::cout <<
				"WARNING: CVtkActorToOpenSG::ProcessGeometryNormalsAndColorsPerVertex() number of normals"
				          << std::endl;
				std::cout << "should equal the number of vertices (points)!" <<
				std::endl << std::endl;
			}
		}

		//possibly getting the colors
		if (m_iColorType == PER_VERTEX)
		{
			vtkIdType iNumColors = vtkColors->GetNumberOfTuples();
			beginEditCP(osgColors);
			{
				unsigned char aColor[4];
				for (int i = 0; i < iNumColors; i++)
				{
					vtkColors->GetTupleValue(i, aColor);
					float r = ((float) aColor[0]) / 255.0f;
					float g = ((float) aColor[1]) / 255.0f;
					float b = ((float) aColor[2]) / 255.0f;
					osgColors->addValue(Color3f(r, g, b));
				}
			} endEditCP(osgColors);
			if (iNumColors != m_iNumPoints)
			{
				std::cout <<
				"WARNING: CVtkActorToOpenSG::ProcessGeometryNormalsAndColorsPerVertex() number of colors"
				          << std::endl;
				std::cout << "should equal the number of vertices (points)!" <<
				std::endl << std::endl;
			}
		}

		//possibly getting the texture coordinates. These are alwary per vertex
		if (vtkTexCoords != NULL)
		{
			int numTuples = vtkTexCoords->GetNumberOfTuples();
			for (int i = 0; i < numTuples; i++)
			{
				double texCoords[3];
				vtkTexCoords->GetTuple(i, texCoords);
				osgTexCoords->addValue(Vec2f(texCoords[0], texCoords[1]));
			}
		}

		//getting the cells
		beginEditCP(osgTypes);
		beginEditCP(osgLengths);
		beginEditCP(osgIndices);
		{
			vtkCellArray* pCells;
			vtkIdType npts, * pts;
			int prim;

			prim = 0;
			pCells = pd->GetVerts();
			if (pCells->GetNumberOfCells() > 0)
				for (pCells->InitTraversal(); pCells->GetNextCell(npts, pts);
				     prim++)
				{
					osgLengths->addValue(npts);
					osgTypes->addValue(GL_POINTS);
					for (int i = 0; i < npts; i++)
						osgIndices->addValue(pts[i]);
				}

			prim = 0;
			pCells = pd->GetLines();
			if (pCells->GetNumberOfCells() > 0)
				for (pCells->InitTraversal(); pCells->GetNextCell(npts, pts);
				     prim++)
				{
					osgLengths->addValue(npts);
					osgTypes->addValue(GL_LINE_STRIP);
					for (int i = 0; i < npts; i++)
						osgIndices->addValue(pts[i]);
				}

			prim = 0;
			pCells = pd->GetPolys();
			if (pCells->GetNumberOfCells() > 0)
				for (pCells->InitTraversal(); pCells->GetNextCell(npts, pts);
				     prim++)
				{
					osgLengths->addValue(npts);
					osgTypes->addValue(GL_POLYGON);
					for (int i = 0; i < npts; i++)
						osgIndices->addValue(pts[i]);
				}

			prim = 0;
			pCells = pd->GetStrips();
			if (pCells->GetNumberOfCells() > 0)
				for (pCells->InitTraversal(); pCells->GetNextCell(npts, pts);
				     prim++)
				{
					osgLengths->addValue(npts);
					osgTypes->addValue(GL_TRIANGLE_STRIP);
					for (int i = 0; i < npts; i++)
						osgIndices->addValue(pts[i]);
				}
		} endEditCP(osgIndices);
		endEditCP(osgLengths);
		endEditCP(osgTypes);

		ChunkMaterialPtr material = CreateMaterial(lit, hasTexCoords);
		beginEditCP(osgGeometry);
		{
			osgGeometry->setPositions(osgPoints);
			osgGeometry->setTypes(osgTypes);
			osgGeometry->setLengths(osgLengths);
			osgGeometry->setIndices(osgIndices);
			osgGeometry->setMaterial(material);

			if (m_iNormalType == PER_VERTEX)
				osgGeometry->setNormals(osgNormals);
			if (m_iColorType == PER_VERTEX)
				osgGeometry->setColors(osgColors);
			if (osgTexCoords->getSize() > 0)
				osgGeometry->setTexCoords(osgTexCoords);
		} endEditCP(osgGeometry);

		osgConversionSuccess = true;

		if (_verbose)
			std::cout << "    End ProcessGeometryNormalsAndColorsPerVertex()" <<
			std::endl;
	}
	else
	{
		//Rendering with OpenSG non indexed geometry by copying a lot of attribute data
		if (_verbose)
			std::cout <<
			"Start ProcessGeometryNonIndexedCopyAttributes(int gl_primitive_type)" <<
			std::endl;
		int gl_primitive_type = -1;
		if(m_iNumGLPolygons > 0)
		{
			if(m_iNumGLPolygons != m_iNumGLPrimitives)
				std::cout <<
				"WARNING: vtkActor contains different kind of primitives" <<
				std::endl;
			gl_primitive_type = GL_POLYGON;
			//osgConversionSuccess = this->ProcessGeometryNonIndexedCopyAttributes(GL_POLYGON, pd, osgGeometry);
		}
		else if(m_iNumGLLineStrips > 0)
		{
			if (m_iNumGLLineStrips != m_iNumGLPrimitives)
				std::cout <<
				"WARNING: vtkActor contains different kind of primitives" <<
				std::endl;
			gl_primitive_type = GL_LINE_STRIP;
			//osgConversionSuccess = this->ProcessGeometryNonIndexedCopyAttributes(GL_LINE_STRIP, pd osgGeometry);
		}
		else if(m_iNumGLTriStrips > 0)
		{
			if (m_iNumGLTriStrips != m_iNumGLPrimitives)
				std::cout <<
				"WARNING: vtkActor contains different kind of primitives" <<
				std::endl;
			gl_primitive_type = GL_TRIANGLE_STRIP;
			//osgConversionSuccess = this->ProcessGeometryNonIndexedCopyAttributes(GL_TRIANGLE_STRIP, pd osgGeometry);
		}
		else if (m_iNumGLPoints > 0)
		{
			if (m_iNumGLPoints != m_iNumGLPrimitives)
				std::cout <<
				"WARNING: vtkActor contains different kind of primitives" <<
				std::endl;
			gl_primitive_type = GL_POINTS;
			//osgConversionSuccess = this->ProcessGeometryNonIndexedCopyAttributes(GL_POINTS, pd osgGeometry);
		}
		if(gl_primitive_type != -1)
		{
			vtkCellArray* pCells;
			if (gl_primitive_type == GL_POINTS)
				pCells = pd->GetVerts();
			else if (gl_primitive_type == GL_LINE_STRIP)
				pCells = pd->GetLines();
			else if (gl_primitive_type == GL_POLYGON)
				pCells = pd->GetPolys();
			else if (gl_primitive_type == GL_TRIANGLE_STRIP)
				pCells = pd->GetStrips();
			else
			{
				std::cout <<
				"CVtkActorToOpenSG::ProcessGeometryNonIndexedCopyAttributes(int gl_primitive_type)"
				          << std::endl;
				std::cout <<
				"  was called with non implemented gl_primitive_type!" << std::endl;
			}

			beginEditCP(osgTypes);
			beginEditCP(osgLengths);
			beginEditCP(osgPoints);
			beginEditCP(osgColors);
			beginEditCP(osgNormals);
			{
				int prim = 0;
				if (pCells->GetNumberOfCells() > 0)
				{
					vtkIdType npts, * pts;
					for (pCells->InitTraversal(); pCells->GetNextCell(npts, pts);
					     prim++)
					{
						osgLengths->addValue(npts);
						osgTypes->addValue(GL_POLYGON);
						for (int i = 0; i < npts; i++)
						{
							double* aVertex;
							double* aNormal;
							unsigned char aColor[4];

							aVertex = pd->GetPoint(pts[i]);
							osgPoints->addValue(Vec3f(aVertex[0],
							                          aVertex[1],
							                          aVertex[2]));

							if (m_iNormalType == PER_VERTEX)
							{
								aNormal =
								        vtkNormals->GetTuple(pts[i]);
								osgNormals->addValue(Vec3f(aNormal[
								                                   0
								                           ],
								                           aNormal[1], aNormal[2]));
							}
							else if (m_iNormalType == PER_CELL)
							{
								aNormal = vtkNormals->GetTuple(prim);
								osgNormals->addValue(Vec3f(aNormal[
								                                   0
								                           ],
								                           aNormal[1], aNormal[2]));
							}

							if (m_iColorType == PER_VERTEX)
							{
								vtkColors->GetTupleValue(pts[i],
								                         aColor);
								float r = ((float) aColor[0]) /
								          255.0f;
								float g = ((float) aColor[1]) /
								          255.0f;
								float b = ((float) aColor[2]) /
								          255.0f;
								osgColors->addValue(Color3f(r, g, b));
							}
							else if (m_iColorType == PER_CELL)
							{
								vtkColors->GetTupleValue(prim,
								                         aColor);
								float r = ((float) aColor[0]) /
								          255.0f;
								float g = ((float) aColor[1]) /
								          255.0f;
								float b = ((float) aColor[2]) /
								          255.0f;
								osgColors->addValue(Color3f(r, g, b));
							}
						}
					}
				}
			} endEditCP(osgTypes);
			endEditCP(osgLengths);
			endEditCP(osgPoints);
			endEditCP(osgColors);
			endEditCP(osgNormals);

			//possibly getting the texture coordinates. These are always per vertex
			vtkPoints* points = pd->GetPoints();
			if ((vtkTexCoords != NULL) && (points != NULL))
			{
				int numPoints = points->GetNumberOfPoints();
				int numTexCoords = vtkTexCoords->GetNumberOfTuples();
				if (numPoints == numTexCoords)
				{
					beginEditCP(osgTexCoords);
					{
						int numTuples = vtkTexCoords->GetNumberOfTuples();
						for (int i = 0; i < numTuples; i++)
						{
							double texCoords[3];
							vtkTexCoords->GetTuple(i, texCoords);
							osgTexCoords->addValue(Vec2f(texCoords[0],
							                             texCoords[1]));
						}
					} endEditCP(osgTexCoords);
				}
			}

			ChunkMaterialPtr material = CreateMaterial(lit, hasTexCoords);
			//GeometryPtr geo = Geometry::create();
			beginEditCP(osgGeometry);
			{
				osgGeometry->setPositions(osgPoints);
				osgGeometry->setTypes(osgTypes);
				osgGeometry->setLengths(osgLengths);
				osgGeometry->setMaterial(material);

				if (m_iNormalType != NOT_GIVEN)
					osgGeometry->setNormals(osgNormals);
				if (m_iColorType != NOT_GIVEN)
					osgGeometry->setColors(osgColors);
				if (osgTexCoords->getSize() > 0)
					osgGeometry->setTexCoords(osgTexCoords);
				//geo->setMaterial(getDefaultMaterial());
			} endEditCP(osgGeometry);

			osgConversionSuccess = true;
		}
		if (_verbose)
			std::cout <<
			"    End ProcessGeometryNonIndexedCopyAttributes(int gl_primitive_type)" <<
			std::endl;
	}

	if(!osgConversionSuccess)
	{
		std::cout << "OpenSG converter was not able to convert this actor." << std::endl;
		return false;
	}

	if(m_iNormalType == NOT_GIVEN)
	{
		//GeometryPtr newGeometryPtr = GeometryPtr::dcast(newNodePtr->getCore());
		if((osgGeometry != NullFC) && (m_iColorType == PER_VERTEX))
		{
			std::cout <<
			"WARNING: Normals are missing in the vtk layer, calculating normals per vertex!"
			          << std::endl;
			calcVertexNormals(osgGeometry);
		}
		else if ((osgGeometry != NullFC) && (m_iColorType == PER_CELL))
		{
			std::cout <<
			"WARNING: Normals are missing in the vtk layer, calculating normals per face!"
			          << std::endl;
			calcFaceNormals(osgGeometry);
		}
		else if (osgGeometry != NullFC)
		{
			std::cout <<
			"WARNING: Normals are missing in the vtk layer, calculating normals per vertex!"
			          << std::endl;
			calcVertexNormals(osgGeometry);
		}
	}

	std::cout << "Conversion finished." << std::endl;

	// Add node to root
	beginEditCP(_osgRoot);
	_osgRoot->addChild(osgGeomNode);
	endEditCP(_osgRoot);

	return true;
}

void vtkOsgConverter::SetVerbose(bool value)
{
	_verbose = value;
}

NodePtr vtkOsgConverter::GetOsgNode()
{
	return _osgRoot;
}

TextureChunkPtr vtkOsgConverter::CreateTexture(vtkTexture* vtkTexture)
{
	if(_verbose)
		std::cout << "Calling CreateTexture()" << std::endl;

	// if(! m_bTextureHasChanged)
	// {
	//   if (_verbose)
	//   {
	//     std::cout << "    ... nothing to do" << std::endl;
	//     std::cout << "    End CreateTexture()" << std::endl;
	//   }
	//   //can we check if the actual image has been updated, even if the texture is the same?
	//   return;
	// }
	// else if(m_pvtkTexture == NULL)
	// {
	//   if (_verbose)
	//   {
	//     std::cout << "    ... texture is (still ?) NULL" << std::endl;
	//     std::cout << "    EndCreateTexture()" << std::endl;
	//   }
	//   //the texture has changed but it is NULL now. We should remove the texture from the material
	//   return;
	// }
	// m_bTextureHasChanged = false;

	vtkImageData* imgData = vtkTexture->GetInput();
	int iImgDims[3];
	imgData->GetDimensions(iImgDims);

	vtkPointData* imgPointData = imgData->GetPointData();
	vtkCellData* imgCellData = imgData->GetCellData();

	vtkDataArray* data = NULL;

	if (imgPointData != NULL)
		if (NULL != imgPointData->GetScalars())
		{
			data = imgPointData->GetScalars();
			if (_verbose)
				std::cout << "    Found texture data in point data" << std::endl;
		}

	if (imgCellData != NULL)
		if (NULL != imgCellData->GetScalars())
		{
			data = imgCellData->GetScalars();
			if (_verbose)
				std::cout << "    Found texture data in cell data" << std::endl;
		}

	if (data == NULL)
	{
		std::cout << "    could not load texture data" << std::endl;
		return NullFC;
	}

	int iImgComps = data->GetNumberOfComponents();
	int iImgPixels = data->GetNumberOfTuples();
	if (iImgPixels != (iImgDims[0] * iImgDims[1] * iImgDims[2]))
	{
		std::cout << "Number of pixels in data array seems to be wrong!" << std::endl;
		return NullFC;
	}

	UInt8* newImageData = new UInt8[iImgDims[0] * iImgDims[1] * iImgDims[2] * iImgComps];
	vtkUnsignedCharArray* oldImageUChar = NULL;
	oldImageUChar = dynamic_cast<vtkUnsignedCharArray*>(data);
	int ucharCounter = 0;
	if (oldImageUChar != NULL)
		for (int i = 0; i < iImgPixels; i++)
		{
			unsigned char pixel[4];
			oldImageUChar->GetTupleValue(i, pixel);
			for (int j = 0; j < iImgComps; j++)
				newImageData[ucharCounter + j] = pixel[j];
			ucharCounter += iImgComps;
		}
	else
		std::cout << "Pixel data come in unsupported vtk type" << std::endl;

	ImagePtr osgImage = Image::create();
	beginEditCP(osgImage);
	{
		osgImage->setWidth(iImgDims[0]);
		osgImage->setHeight(iImgDims[1]);
		osgImage->setDepth(iImgDims[2]);
		osgImage->setDataType(Image::OSG_UINT8_IMAGEDATA);
		if (iImgComps == 1)
			osgImage->setPixelFormat(Image::OSG_L_PF);
		else if (iImgComps == 2)
			osgImage->setPixelFormat(Image::OSG_LA_PF);
		else if (iImgComps == 3)
			osgImage->setPixelFormat(Image::OSG_RGB_PF);
		else if (iImgComps == 4)
			osgImage->setPixelFormat(Image::OSG_RGBA_PF);
		else
		{
			std::cout << "Unsupported image type!" << std::endl;
			delete [] newImageData;
			return NullFC;
		}
		osgImage->setData(newImageData);
	} endEditCP(osgImage);

	TextureChunkPtr osgTextureChunk = TextureChunk::create();
	beginEditCP(osgTextureChunk);
	{
		osgTextureChunk->setImage(osgImage);
	} endEditCP(osgTextureChunk);

	if (_verbose)
	{
		std::cout << "    Loading image with " << iImgDims[0] << " x " << iImgDims[1] <<
		" x " << iImgDims[2] << " pixels." << std::endl;
		std::cout << "    Components: " << iImgComps << std::endl;
		std::cout << "End CreateTexture()" << std::endl;
	}

	return osgTextureChunk;
}

ChunkMaterialPtr vtkOsgConverter::CreateMaterial(bool lit, bool hasTexCoords)
{
	if (_verbose)
		std::cout << "Start CreateMaterial()" << std::endl;

	vtkProperty* prop = _actor->GetProperty();
	double* diffuseColor = prop->GetDiffuseColor();
	double* ambientColor = prop->GetAmbientColor();
	double* specularColor = prop->GetSpecularColor();
	double specularPower = prop->GetSpecularPower();

	double diffuse = prop->GetDiffuse();
	double ambient = prop->GetAmbient();
	double specular = prop->GetSpecular();
	double opacity = prop->GetOpacity();

	float pointSize = prop->GetPointSize();
	float lineWidth = prop->GetLineWidth();
	// int lineStipplePattern = prop->GetLineStipplePattern();

	int representation = prop->GetRepresentation();

	if (_verbose)
	{
		std::cout << "    Colors:" << std::endl;
		std::cout << "    diffuse " << diffuse << " * " << diffuseColor[0] << " " <<
		diffuseColor[1] << " " << diffuseColor[2] << std::endl;
		std::cout << "    ambient " << ambient << " * " << ambientColor[0] << " " <<
		ambientColor[1] << " " << ambientColor[2] << std::endl;
		std::cout << "    specular " << specular << " * " << specularColor[0] << " " <<
		specularColor[1] << " " << specularColor[2] << std::endl;
	}

	PolygonChunkPtr polygonChunk = PolygonChunk::create();
	beginEditCP(polygonChunk);
	{
		if (representation == VTK_SURFACE)
		{
			polygonChunk->setFrontMode(GL_FILL);
			polygonChunk->setBackMode(GL_FILL);
		}
		else if (representation == VTK_WIREFRAME)
		{
			polygonChunk->setFrontMode(GL_LINE);
			polygonChunk->setBackMode(GL_LINE);
		}
		else
		{
			polygonChunk->setFrontMode(GL_POINT);
			polygonChunk->setBackMode(GL_POINT);
		}
	} endEditCP(polygonChunk);

	MaterialChunkPtr osgMaterialChunk = MaterialChunk::create();
	beginEditCP(osgMaterialChunk);
	{
		osgMaterialChunk->setDiffuse(Color4f(diffuseColor[0] * diffuse, diffuseColor[1] *
		                                     diffuse, diffuseColor[2] * diffuse, opacity));
		osgMaterialChunk->setSpecular(Color4f(specularColor[0] * specular,
		                                      specularColor[1] * specular,
		                                      specularColor[2] * specular, 1.0));
		osgMaterialChunk->setAmbient(Color4f(ambientColor[0] * ambient, ambientColor[1] *
		                                     ambient, ambientColor[2] * ambient, opacity));
		osgMaterialChunk->setShininess(specularPower);

		//if(opacity < 1.0)
		//{
		// osgMaterialChunk->setColorMaterial(GL_AMBIENT); // HACK: Opacity does not work with GL_AMBIENT_AND_DIFFUSE
		//osgMaterialChunk->setTransparency(1.0f - opacity);
		//}
		//else
		osgMaterialChunk->setColorMaterial(GL_AMBIENT_AND_DIFFUSE);

		// On objects consisting only of points or lines, dont lit
		if(!lit)
			osgMaterialChunk->setLit(false);
	} endEditCP(osgMaterialChunk);

	ChunkMaterialPtr osgChunkMaterial = ChunkMaterial::create();
	beginEditCP(osgChunkMaterial);
	{
		osgChunkMaterial->addChunk(osgMaterialChunk);
		osgChunkMaterial->addChunk(TwoSidedLightingChunk::create());
		osgChunkMaterial->addChunk(polygonChunk);

		if(pointSize > 1.0f)
		{
			PointChunkPtr pointChunk = PointChunk::create();
			pointChunk->setSize(pointSize);
			osgChunkMaterial->addChunk(pointChunk);
		}

		if(lineWidth > 1.0f)
		{
			LineChunkPtr lineChunk = LineChunk::create();
			lineChunk->setWidth(lineWidth);
			osgChunkMaterial->addChunk(lineChunk);
		}

		// TEXTURE
		if (hasTexCoords)
		{
			vtkTexture* vtkTexture = _actor->GetTexture();
			if (vtkTexture)
			{
				TextureChunkPtr osgTextureChunk = NullFC;
				osgTextureChunk = CreateTexture(vtkTexture);

				if(osgTextureChunk != NullFC)
				{
					if (_verbose)
						std::cout << "    Add TextureChunk" << std::endl;
					osgChunkMaterial->addChunk(osgTextureChunk);
				}
			}
		}
	} endEditCP(osgChunkMaterial);

	if (_verbose)
		std::cout << "    End CreateMaterial()" << std::endl;

	return osgChunkMaterial;
}
