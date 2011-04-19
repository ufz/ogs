/////////////////////////////////////////////////////////////////////////////////////////
//                       By Bjoern Zehner 2006/2010                                    //
//			UFZ-Helmholtz Center for Environmental Research Leipzig                    //
//			      Permoser Str. 15, D-04318 Leizig, Germany                            //
/////////////////////////////////////////////////////////////////////////////////////////
// Class for translating VTK Objects into OpenSG representations, derives from vtkActor//
// See Project/file VtkOsgActorTestXXX for an example for the use                      //


#include "vtkOsgActor.h"

#include "vtkMapper.h"
#include "vtkDataSet.h"
#include "vtkDataSetMapper.h"



vtkOsgActor::vtkOsgActor(void) : vtkOpenGLActor(),
	m_pvtkNormals(NULL),
	m_pvtkTexCoords(NULL),
	m_pvtkColors(NULL),
	m_pvtkTexture(NULL),
	m_bTextureHasChanged(false),
	m_iColorType(NOT_GIVEN),
	m_iNormalType(NOT_GIVEN),
	m_bVerbose(false),
	m_iNumPoints(0),
	m_iNumNormals(0),
	m_iNumColors(0),
	m_iNumGLPoints(0),
	m_iNumGLLineStrips(0),
	m_iNumGLPolygons(0),
	m_iNumGLTriStrips(0),
	m_iNumGLPrimitives(0),
	m_posgRoot(NullFC),
	m_posgTransform(NullFC),
	m_posgGeomNode(NullFC),
	m_posgGeometry(NullFC),
	m_posgMaterial(NullFC),
	m_posgMaterialChunk(NullFC),
	m_posgTextureChunk(NullFC),
	m_posgPolygonChunk(NullFC),
	m_posgImage(NullFC),
	m_posgTypes(NullFC),
	m_posgLengths(NullFC),
	m_posgIndices(NullFC),
	m_posgPoints(NullFC),
	m_posgColors(NullFC),
	m_posgNormals(NullFC),
	m_posgTexCoords(NullFC)
{
	TransformPtr tptr;
	m_posgRoot = makeCoredNode<osg::Transform>(&tptr);
	m_posgTransform = tptr;
}

vtkOsgActor::~vtkOsgActor(void)
{
	m_posgRoot = NullFC;

	//Open SG Objects are deleted via the reference counting scheme
	ClearOsg();
}

void vtkOsgActor::InitOpenSG(){
	m_posgGeomNode = Node::create();
	m_posgGeometry = Geometry::create();
	m_posgMaterial = ChunkMaterial::create();
	m_posgMaterialChunk = MaterialChunk::create();
	m_posgTextureChunk = TextureChunk::create();
	m_posgPolygonChunk = PolygonChunk::create();
	m_posgImage = Image::create();
	beginEditCP(m_posgRoot);{
		m_posgRoot->addChild(m_posgGeomNode);
	};endEditCP(m_posgRoot);
	beginEditCP(m_posgGeomNode);{
		m_posgGeomNode->setCore(m_posgGeometry);
	};endEditCP(m_posgGeomNode);
	m_posgTypes = GeoPTypesUI8::create();
	m_posgLengths = GeoPLengthsUI32::create();
	m_posgIndices = GeoIndicesUI32::create();
	m_posgPoints = GeoPositions3f::create();
	m_posgColors = GeoColors3f::create();
	m_posgNormals = GeoNormals3f::create();
	m_posgTexCoords = GeoTexCoords2d::create();
}

vtkOsgActor *vtkOsgActor::New(){
	// First try to create the object from the vtkGraphicsFactory
	vtkObject* ret = vtkObjectFactory::CreateInstance("vtkOsgActor");
  
	if (ret) {
		return static_cast<vtkOsgActor*>(ret);
	}
	return new vtkOsgActor;
}

void vtkOsgActor::PrintSelf(ostream& os, vtkIndent indent){
	vtkOpenGLActor::PrintSelf(os, indent);

	// beside this nothing is done so far ...
}

void vtkOsgActor::Render(vtkRenderer *ren, vtkMapper *mapper){
	vtkOpenGLActor::Render(ren, mapper);
}

void vtkOsgActor::UpdateOsg(){
	if (m_posgGeomNode == NullFC) InitOpenSG();
	if (m_bTextureHasChanged) CreateTexture();

	this->GetMapper()->Update();

	LookForNormals();
	LookForColors();
	LookForTexCoords();
	LookForArraySizes();
	CreateTexture();

	double scaling[3];
	//double translation[3];
	//double rotation[3];

	//this->GetPosition(translation);
	this->GetScale(scaling);

	if (m_bVerbose){
		std::cout << "set scaling: " << scaling[0] << " " << scaling[1] << " " << scaling[2] << std::endl;
	}

	osg::Matrix m;
	m.setIdentity();
	m.setScale(scaling[0], scaling[1], scaling[2]);
	beginEditCP(m_posgTransform);{
		m_posgTransform->setMatrix(m);
	};endEditCP(m_posgTransform);

	NodePtr newNodePtr = this->GetNodePtr();

	if (m_iNormalType == NOT_GIVEN){
		GeometryPtr newGeometryPtr = GeometryPtr::dcast(newNodePtr->getCore());
		if ((newGeometryPtr != NullFC) && (m_iColorType == PER_VERTEX)){
			std::cout << "WARNING: Normals are missing in the vtk layer, calculating normals per vertex!" << std::endl;
			calcVertexNormals(newGeometryPtr);
		} else if ((newGeometryPtr != NullFC) && (m_iColorType == PER_CELL)){
			std::cout << "WARNING: Normals are missing in the vtk layer, calculating normals per face!" << std::endl;
			calcFaceNormals(newGeometryPtr);
		} else if (newGeometryPtr != NullFC){
			std::cout << "WARNING: Normals are missing in the vtk layer, calculating normals per vertex!" << std::endl;
			calcVertexNormals(newGeometryPtr);
		}
	}

	beginEditCP(m_posgRoot);{
		m_posgRoot->addChild(newNodePtr);
	};endEditCP(m_posgRoot);
}

void vtkOsgActor::ClearOsg(){
	//This also decrements the reference count, possibly deleting the objects
	m_posgGeomNode = NullFC;
	m_posgGeometry = NullFC;
	m_posgMaterial = NullFC;
	m_posgMaterialChunk = NullFC;
	m_posgTextureChunk = NullFC;
	m_posgPolygonChunk = NullFC;
	m_posgImage = NullFC;
	m_posgTypes = NullFC;
	m_posgLengths = NullFC;
	m_posgIndices = NullFC;
	m_posgPoints = NullFC;
	m_posgColors = NullFC;
	m_posgNormals = NullFC;
	m_posgTexCoords = NullFC;
	m_bTextureHasChanged = true;
}

void vtkOsgActor::SetVerbose(bool value){
	m_bVerbose = value;
}

void vtkOsgActor::SetTexture(vtkTexture *vtkTex){
	m_pvtkTexture = vtkTex;
	m_bTextureHasChanged = true;
	vtkOpenGLActor::SetTexture(vtkTex);
}

NodePtr vtkOsgActor::GetOsgRoot(){
	return m_posgRoot;
}

void vtkOsgActor::LookForNormals(){
	vtkPolyData *pPolyData = NULL;
	if (dynamic_cast<vtkPolyDataMapper*>(this->GetMapper())){
		pPolyData = (vtkPolyData*) this->GetMapper()->GetInput();
		if (m_bVerbose){
			std::cerr << "Using vtkPolyDataMapper directly" << std::endl;
		}
	} else if (dynamic_cast<vtkDataSetMapper*>(this->GetMapper())){
		vtkDataSetMapper *dataSetMapper = (vtkDataSetMapper*) this->GetMapper();
		pPolyData = (vtkPolyData*) dataSetMapper->GetPolyDataMapper()->GetInput();
		if (m_bVerbose){
			std::cerr << "Using vtkPolyDataMapper via the vtkDataSetMapper" << std::endl;
		}
	}

	m_iNormalType = NOT_GIVEN;
	if (pPolyData == NULL) return;

	if (this->GetProperty()->GetInterpolation() == VTK_FLAT){
		m_pvtkNormals = pPolyData->GetCellData()->GetNormals();
		if (m_pvtkNormals != NULL) m_iNormalType = PER_CELL;
	}else{
		m_pvtkNormals = pPolyData->GetPointData()->GetNormals();
		if (m_pvtkNormals != NULL) m_iNormalType = PER_VERTEX;
	}
	if (m_bVerbose){
		if (m_iNormalType != NOT_GIVEN){
			std::cerr << "	number of normals: " << m_pvtkNormals->GetNumberOfTuples() << std::endl;
			std::cerr << "	normals are given: ";
			std::cerr << ((m_iNormalType == PER_VERTEX) ? "per vertex" : "per cell") << std::endl;
		}else{
			std::cerr << "	no normals are given" << std::endl;
		}
	}
}

void vtkOsgActor::LookForColors(){
	vtkPolyData *pPolyData = NULL;
	vtkPolyDataMapper *pPolyDataMapper = NULL;
	if (dynamic_cast<vtkPolyDataMapper*>(this->GetMapper())){
		pPolyDataMapper = (vtkPolyDataMapper*) this->GetMapper();
		pPolyData = (vtkPolyData*) pPolyDataMapper->GetInput();
		if (m_bVerbose){
			std::cerr << "Using vtkPolyDataMapper directly" << std::endl;
		}
	} else if (dynamic_cast<vtkDataSetMapper*>(this->GetMapper())){
		vtkDataSetMapper *dataSetMapper = (vtkDataSetMapper*) this->GetMapper();
		pPolyDataMapper = dataSetMapper->GetPolyDataMapper();
		pPolyData = (vtkPolyData*) pPolyDataMapper->GetInput();
		if (m_bVerbose){
			std::cerr << "Using vtkPolyDataMapper via the vtkDataSetMapper" << std::endl;
		}
	}
	m_iColorType = NOT_GIVEN;
	if (pPolyData == NULL) return;
	if (pPolyDataMapper == NULL) return;

	//m_pvtkColors = pPolyDataMapper->MapScalars(1.0);
	//if (pPolyDataMapper->GetScalarVisibility() && (m_pvtkColors != NULL)){
	if (pPolyDataMapper->GetScalarVisibility()){
		int iScalarMode = pPolyDataMapper->GetScalarMode();
		m_pvtkColors = pPolyDataMapper->MapScalars(1.0);
		if (m_pvtkColors == NULL){
			m_iColorType = NOT_GIVEN;
		} else if (iScalarMode == VTK_SCALAR_MODE_USE_CELL_DATA){
			m_iColorType = PER_CELL;
		} else if (iScalarMode == VTK_SCALAR_MODE_USE_POINT_DATA){
			m_iColorType = PER_VERTEX;
		} else if (iScalarMode == VTK_SCALAR_MODE_USE_CELL_FIELD_DATA){
			std::cerr << "WARNING: Can not process colours with scalar mode using cell field data!" << std::endl;
			m_iColorType = NOT_GIVEN;
		} else if (iScalarMode == VTK_SCALAR_MODE_USE_POINT_FIELD_DATA){
			std::cerr << "WARNING: Can not process colours with scalar mode using point field data!" << std::endl;
			m_iColorType = NOT_GIVEN;
		} else if (iScalarMode == VTK_SCALAR_MODE_DEFAULT){
			//Bummer, we do not know what it is. may be we can make a guess
			int numColors = m_pvtkColors->GetNumberOfTuples();
			vtkPoints *points = pPolyData->GetPoints();
			int numPoints = points->GetNumberOfPoints();
			int numCells = pPolyData->GetVerts()->GetNumberOfCells();
			numCells += pPolyData->GetLines()->GetNumberOfCells();
			numCells += pPolyData->GetPolys()->GetNumberOfCells();
			numCells += pPolyData->GetStrips()->GetNumberOfCells();
			if (numColors == 0){
				m_iColorType = NOT_GIVEN;
			} else if (numColors == numPoints){
				m_iColorType = PER_VERTEX;
			} else if (numColors == numCells){
				m_iColorType = PER_CELL;
			} else {
				m_iColorType = NOT_GIVEN;
			}
		}
	}
	if (m_bVerbose){	
		if (m_iColorType != NOT_GIVEN){
			std::cerr << "  number of colors: " << m_pvtkColors->GetNumberOfTuples() << std::endl;
			std::cerr << "	colors are given: ";
			std::cerr << ((m_iColorType == PER_VERTEX) ? "per vertex" : "per cell") << std::endl;
		}else{
			std::cerr << "	no colors are given" << std::endl;
		}
	}
}

void vtkOsgActor::LookForTexCoords(){
	vtkPolyData *pPolyData = NULL;
	if (dynamic_cast<vtkPolyDataMapper*>(this->GetMapper())){
		pPolyData = (vtkPolyData*) this->GetMapper()->GetInput();
		if (m_bVerbose){
			std::cerr << "Using vtkPolyDataMapper directly" << std::endl;
		}
	} else if (dynamic_cast<vtkDataSetMapper*>(this->GetMapper())){
		vtkDataSetMapper *dataSetMapper = (vtkDataSetMapper*) this->GetMapper();
		pPolyData = (vtkPolyData*) dataSetMapper->GetPolyDataMapper()->GetInput();
		if (m_bVerbose){
			std::cerr << "Using vtkPolyDataMapper via the vtkDataSetMapper" << std::endl;
		}
	}
	if (pPolyData == NULL) return;

	vtkPointData *pointData = pPolyData->GetPointData();
	if (pointData == NULL) return;
	m_pvtkTexCoords = pointData->GetTCoords();
	if (m_bVerbose){
		if (m_pvtkTexCoords != NULL){
			std::cerr << "	found tex-coords: " << m_pvtkTexCoords->GetNumberOfTuples() << std::endl;
		} else {
			std::cerr << "	no tex-coords where given" << std::endl;
		}
	}
}

void vtkOsgActor::LookForArraySizes(){
	vtkPolyData *pPolyData = NULL;
	if (dynamic_cast<vtkPolyDataMapper*>(this->GetMapper())){
		pPolyData = (vtkPolyData*) this->GetMapper()->GetInput();
		if (m_bVerbose){
			std::cerr << "Using vtkPolyDataMapper directly" << std::endl;
		}
	} else if (dynamic_cast<vtkDataSetMapper*>(this->GetMapper())){
		vtkDataSetMapper *dataSetMapper = (vtkDataSetMapper*) this->GetMapper();
		pPolyData = (vtkPolyData*) dataSetMapper->GetPolyDataMapper()->GetInput();
		if (m_bVerbose){
			std::cerr << "Using vtkPolyDataMapper via the vtkDataSetMapper" << std::endl;
		}
	}
	if (pPolyData == NULL) return;

	m_iNumPoints = pPolyData->GetNumberOfPoints();
	m_iNumGLPoints = pPolyData->GetVerts()->GetNumberOfCells();
	m_iNumGLLineStrips = pPolyData->GetLines()->GetNumberOfCells();
	m_iNumGLPolygons = pPolyData->GetPolys()->GetNumberOfCells();
	m_iNumGLTriStrips = pPolyData->GetStrips()->GetNumberOfCells();
	m_iNumGLPrimitives = m_iNumGLPoints + m_iNumGLLineStrips + m_iNumGLPolygons + m_iNumGLTriStrips; 


	if (m_bVerbose){
		std::cerr << "	number of vertices: " << m_iNumPoints << std::endl;
		std::cerr << "	number of GL_POINTS: " << m_iNumGLPoints << std::endl;
		std::cerr << "	number of GL_LINE_STRIPS: " << m_iNumGLLineStrips << std::endl;
		std::cerr << "	number of GL_POLYGON's: " << m_iNumGLPolygons << std::endl;
		std::cerr << "	number of GL_TRIANGLE_STRIPS: " << m_iNumGLTriStrips << std::endl;
		std::cerr << "	number of primitives: " << m_iNumGLPrimitives << std::endl;
	}
}

void vtkOsgActor::CreateTexture(){
	if (m_bVerbose){
		std::cerr << "Calling CreateTexture()" << std::endl;
	}
	if (! m_bTextureHasChanged){
		if (m_bVerbose){
			std::cerr << "		... nothing to do" << std::endl;
			std::cerr << "		End CreateTexture()" << std::endl;
		}
		//can we check if the actual image has been updated, even if the texture is the same?
		return;
	} else if (m_pvtkTexture == NULL){
		if (m_bVerbose){
			std::cerr << "		... texture is (still ?) NULL" << std::endl;
			std::cerr << "		EndCreateTexture()" << std::endl;
		}
			//the texture has changed but it is NULL now. We should remove the texture from the material
		return;
	}
	m_bTextureHasChanged = false;

	vtkImageData *imgData = m_pvtkTexture->GetInput();
	int iImgDims[3];
	imgData->GetDimensions(iImgDims);

	vtkPointData *imgPointData = imgData->GetPointData();
	vtkCellData *imgCellData = imgData->GetCellData();

	vtkDataArray *data = NULL;

	if (imgPointData != NULL){
		if (NULL != imgPointData->GetScalars()){
			data = imgPointData->GetScalars();
			if (m_bVerbose) std::cout << "		found texture data in point data" << std::endl;
		}
	}

	if (imgCellData != NULL){
		if (NULL != imgCellData->GetScalars()){
			data = imgCellData->GetScalars();
			if (m_bVerbose) std::cout << "		found texture data in cell data" << std::endl;
		}
	}

	if (data == NULL){
		std::cerr << "		could not load texture data" << std::endl;
		return;
	}
	
	int iImgComps = data->GetNumberOfComponents();
	int iImgPixels = data->GetNumberOfTuples();
	if (iImgPixels != (iImgDims[0] * iImgDims[1] * iImgDims[2])){
		std::cout << "Number of pixels in data array seems to be wrong!" << std::endl;
		return;
	}

	UInt8 *newImageData = new UInt8[iImgDims[0] * iImgDims[1] * iImgDims[2] * iImgComps];
	vtkUnsignedCharArray *oldImageUChar = NULL;
	oldImageUChar = dynamic_cast<vtkUnsignedCharArray*>(data);
	int ucharCounter = 0;
	if (oldImageUChar != NULL){
		for (int i=0; i<iImgPixels; i++){
			unsigned char pixel[4];
			oldImageUChar->GetTupleValue(i, pixel);
			for (int j=0; j<iImgComps; j++){
				newImageData[ucharCounter + j] = pixel[j];
			}
			ucharCounter += iImgComps;
		}
	} else {
		std::cout << "Pixel data come in unsupported vtk type" << std::endl;
	}
	
	beginEditCP(m_posgImage);{
		m_posgImage->setWidth(iImgDims[0]);
		m_posgImage->setHeight(iImgDims[1]);
		m_posgImage->setDepth(iImgDims[2]);
		m_posgImage->setDataType(Image::OSG_UINT8_IMAGEDATA);
		if (iImgComps == 1){
			m_posgImage->setPixelFormat(Image::OSG_L_PF);
		} else if (iImgComps == 3){
			m_posgImage->setPixelFormat(Image::OSG_RGB_PF);
		} else if (iImgComps == 4){
			m_posgImage->setPixelFormat(Image::OSG_RGBA_PF);
		} else {
			std::cout << "Unsupported image type!" << std::endl;
			delete [] newImageData;
			return;
		}
		m_posgImage->setData(newImageData);
	};endEditCP(m_posgImage);

	beginEditCP(m_posgTextureChunk);{
		m_posgTextureChunk->setImage(m_posgImage.get());
	};endEditCP(m_posgTextureChunk);

	if (m_bVerbose){
		std::cerr << "		Loading image with " << iImgDims[0] << " x " << iImgDims[1] << " x " << iImgDims[2] << "pixels." << std::endl;
		std::cerr << "		components: " << iImgComps << std::endl;
		std::cerr << "		End CreateTexture()" << std::endl;
	}
}

ChunkMaterialPtr vtkOsgActor::CreateMaterial(){
	if (m_bVerbose){
		std::cerr << "Start CreateMaterial()" << std::endl;
	}
	vtkProperty *prop = this->GetProperty();
	double *diffuseColor = prop->GetDiffuseColor();
	double *ambientColor = prop->GetAmbientColor();
	double *specularColor = prop->GetSpecularColor();
	double specularPower = prop->GetSpecularPower();

	double diffuse = prop->GetDiffuse();
	double ambient = prop->GetAmbient();
	double specular = prop->GetSpecular();

	//float opacity = prop->GetOpacity();
	int representation = prop->GetRepresentation();

	if (m_bVerbose){
		std::cerr << "		Colors:" << std::endl;
		std::cerr << "		diffuse " << diffuse << " * " << diffuseColor[0] << " " << diffuseColor[1] << " " << diffuseColor[2] << std::endl;
		std::cerr << "		ambient " << ambient << " * " << ambientColor[0] << " " << ambientColor[1] << " " << ambientColor[2] << std::endl;
		std::cerr << "		specular " << specular << " * " << specularColor[0] << " " << specularColor[1] << " " << specularColor[2] << std::endl;
	}

	beginEditCP(m_posgPolygonChunk);{
		if (representation == VTK_SURFACE){
			m_posgPolygonChunk->setFrontMode(GL_FILL);
			m_posgPolygonChunk->setBackMode(GL_FILL);
		}else if (representation == VTK_WIREFRAME){
			m_posgPolygonChunk->setFrontMode(GL_LINE);
			m_posgPolygonChunk->setBackMode(GL_LINE);
		}else{
			m_posgPolygonChunk->setFrontMode(GL_POINT);
			m_posgPolygonChunk->setBackMode(GL_POINT);
		}
	};endEditCP(m_posgPolygonChunk);

	beginEditCP(m_posgMaterialChunk);{
		m_posgMaterialChunk->setDiffuse(Color4f(diffuseColor[0]*diffuse, diffuseColor[1]*diffuse, diffuseColor[2]*diffuse, 1.0));
		m_posgMaterialChunk->setSpecular(Color4f(specularColor[0]*specular, specularColor[1]*specular, specularColor[2]*specular, 1.0));
		m_posgMaterialChunk->setAmbient(Color4f(ambientColor[0]*ambient, ambientColor[1]*ambient, ambientColor[2]*ambient, 1.0));
		m_posgMaterialChunk->setShininess(specularPower);
		m_posgMaterialChunk->setColorMaterial(GL_AMBIENT_AND_DIFFUSE);
		//m_posgMaterialChunk->setTransparency(1.0f - opacity); // 1-opacity ?
	};endEditCP(m_posgMaterialChunk);

	beginEditCP(m_posgMaterial);{
		m_posgMaterial->addChunk(m_posgMaterialChunk);
		m_posgMaterial->addChunk(TwoSidedLightingChunk::create());
		m_posgMaterial->addChunk(m_posgPolygonChunk);

		if (m_pvtkTexCoords != NULL){
			if (m_bVerbose){
				std::cerr << "		Add TextureChunk" << std::endl;
			}
			m_posgMaterial->addChunk(m_posgTextureChunk);
		} else {
			if (m_bVerbose){
				std::cerr << "		Not adding TextureChunk" << std::endl;
			}
		}
	}endEditCP(m_posgMaterial);
	return m_posgMaterial;
	if (m_bVerbose){
		std::cerr << "		End CreateMaterial()" << std::endl;
	}
}

NodePtr vtkOsgActor::ProcessGeometryNormalsAndColorsPerVertex(){
	if (m_bVerbose){
		std::cerr << "Start ProcessGeometryNormalsAndColorsPerVertex()" << std::endl;
	}

	beginEditCP(m_posgTypes);{
		m_posgTypes->clear();
	};endEditCP(m_posgTypes);

	beginEditCP(m_posgLengths);{
		m_posgLengths->clear();
	};endEditCP(m_posgLengths);

	beginEditCP(m_posgIndices);{
		m_posgIndices->clear();
	};endEditCP(m_posgIndices);

	beginEditCP(m_posgPoints);{
		m_posgPoints->clear();
	};endEditCP(m_posgPoints);

	beginEditCP(m_posgColors);{
		m_posgColors->clear();
	};endEditCP(m_posgColors);

	beginEditCP(m_posgNormals);{
		m_posgNormals->clear();
	};endEditCP(m_posgNormals);

	beginEditCP(m_posgTexCoords);{
		m_posgTexCoords->clear();
	};endEditCP(m_posgTexCoords);

	int iNumPoints = 0;
	int iNumNormals = 0;
	int iNumColors = 0;
	int i;

	vtkPolyData *pPolyData = NULL;
	if (dynamic_cast<vtkPolyDataMapper*>(this->GetMapper())){
		pPolyData = (vtkPolyData*) this->GetMapper()->GetInput();
		if (m_bVerbose){
			std::cerr << "		Using vtkPolyDataMapper directly" << std::endl;
		}
	} else if (dynamic_cast<vtkDataSetMapper*>(this->GetMapper())){
		vtkDataSetMapper *dataSetMapper = (vtkDataSetMapper*) this->GetMapper();
		pPolyData = (vtkPolyData*) dataSetMapper->GetPolyDataMapper()->GetInput();
		if (m_bVerbose){
			std::cerr << "		Using vtkPolyDataMapper via the vtkDataSetMapper" << std::endl;
		}
	}
	if (pPolyData == NULL) return NullFC;

	//getting the vertices:
	beginEditCP(m_posgPoints);{
		iNumPoints = pPolyData->GetNumberOfPoints();
		for (i=0; i<iNumPoints; i++){
			double *aVertex = pPolyData->GetPoint(i);
			m_posgPoints->addValue(Vec3f(aVertex[0], aVertex[1], aVertex[2]));
		}
	}endEditCP(m_posgPoints);

	//possibly getting the normals
	if (m_iNormalType == PER_VERTEX){
		iNumNormals = m_pvtkNormals->GetNumberOfTuples();
		beginEditCP(m_posgNormals);{
			double *aNormal;
			for (i=0; i<iNumNormals; i++){
				aNormal = m_pvtkNormals->GetTuple(i);
				m_posgNormals->addValue(Vec3f(aNormal[0], aNormal[1], aNormal[2]));
			}
		}endEditCP(m_posgNormals);
		if (iNumNormals != iNumPoints){
			std::cerr << "WARNING: CVtkActorToOpenSG::ProcessGeometryNormalsAndColorsPerVertex() number of normals" << std::endl;
			std::cerr << "should equal the number of vertices (points)!" << std::endl << std::endl;
		}
	}

	//possibly getting the colors
	if (m_iColorType == PER_VERTEX){
		iNumColors = m_pvtkColors->GetNumberOfTuples();
		beginEditCP(m_posgColors);{
			unsigned char aColor[4];
			for (i=0; i<iNumColors; i++){
				m_pvtkColors->GetTupleValue(i, aColor);
				float r = ((float) aColor[0]) / 255.0f;
				float g = ((float) aColor[1]) / 255.0f;
				float b = ((float) aColor[2]) / 255.0f;
				m_posgColors->addValue(Color3f(r, g, b));
			}
		}endEditCP(m_posgColors);
		if (iNumColors != iNumPoints){
			std::cerr << "WARNING: CVtkActorToOpenSG::ProcessGeometryNormalsAndColorsPerVertex() number of colors" << std::endl;
			std::cerr << "should equal the number of vertices (points)!" << std::endl << std::endl;
		}
	}

	//possibly getting the texture coordinates. These are alwary per vertex
	if (m_pvtkTexCoords != NULL){
		int numTuples = m_pvtkTexCoords->GetNumberOfTuples();
		for (i=0; i<numTuples; i++){
			double texCoords[3];
			m_pvtkTexCoords->GetTuple(i, texCoords);
			m_posgTexCoords->addValue(Vec2f(texCoords[0], texCoords[1]));
		}
	}

	//getting the cells
	beginEditCP(m_posgTypes);
	beginEditCP(m_posgLengths);
	beginEditCP(m_posgIndices);{
		vtkCellArray *pCells;
		vtkIdType npts, *pts;
		int prim;

		prim = 0;
		pCells = pPolyData->GetVerts();
		if (pCells->GetNumberOfCells() > 0){
			for (pCells->InitTraversal(); pCells->GetNextCell(npts, pts); prim++){
				m_posgLengths->addValue(npts);
				m_posgTypes->addValue(GL_POINTS);
				for (i=0; i<npts; i++){
					m_posgIndices->addValue(pts[i]);
				}
			}
		}

		prim = 0;
		pCells = pPolyData->GetLines();
		if (pCells->GetNumberOfCells() > 0){
			for (pCells->InitTraversal(); pCells->GetNextCell(npts, pts); prim++){
				m_posgLengths->addValue(npts);
				m_posgTypes->addValue(GL_LINE_STRIP);
				for (i=0; i<npts; i++){
					m_posgIndices->addValue(pts[i]);
				}
			}
		}

		prim = 0;
		pCells = pPolyData->GetPolys();
		if (pCells->GetNumberOfCells() > 0){
			for (pCells->InitTraversal(); pCells->GetNextCell(npts, pts); prim++){
				m_posgLengths->addValue(npts);
				m_posgTypes->addValue(GL_POLYGON);
				for (i=0; i<npts; i++){
					m_posgIndices->addValue(pts[i]);
				}
			}
		}

		prim = 0;
		pCells = pPolyData->GetStrips();
		if (pCells->GetNumberOfCells() > 0){
			for (pCells->InitTraversal(); pCells->GetNextCell(npts, pts); prim++){
				m_posgLengths->addValue(npts);
				m_posgTypes->addValue(GL_TRIANGLE_STRIP);
				for (i=0; i<npts; i++){
					m_posgIndices->addValue(pts[i]);
				}
			}
		}
	}endEditCP(m_posgIndices);
	endEditCP(m_posgLengths);
	endEditCP(m_posgTypes);

	ChunkMaterialPtr material = CreateMaterial();
	beginEditCP(m_posgGeometry);{
		m_posgGeometry->setPositions(m_posgPoints);
		m_posgGeometry->setTypes(m_posgTypes);
		m_posgGeometry->setLengths(m_posgLengths);
		m_posgGeometry->setIndices(m_posgIndices);
		m_posgGeometry->setMaterial(material);

		if (m_iNormalType == PER_VERTEX) m_posgGeometry->setNormals(m_posgNormals);
		if (m_iColorType == PER_VERTEX) m_posgGeometry->setColors(m_posgColors);
		if (m_posgTexCoords->getSize() > 0) m_posgGeometry->setTexCoords(m_posgTexCoords);
	};endEditCP(m_posgGeometry);
	if (m_bVerbose){
		std::cerr << "		Start ProcessGeometryNormalsAndColorsPerVertex()" << std::endl;
	}
	return m_posgGeomNode;
}

NodePtr vtkOsgActor::ProcessGeometryNonIndexedCopyAttributes(int gl_primitive_type){
	if (m_bVerbose){
		std::cout << "Start ProcessGeometryNonIndexedCopyAttributes(int gl_primitive_type)" << std::endl;
	}

	beginEditCP(m_posgTypes);{
		m_posgTypes->clear();
	};endEditCP(m_posgTypes);

	beginEditCP(m_posgLengths);{
		m_posgLengths->clear();
	};endEditCP(m_posgLengths);

	beginEditCP(m_posgIndices);{
		m_posgIndices->clear();
	};endEditCP(m_posgIndices);

	beginEditCP(m_posgPoints);{
		m_posgPoints->clear();
	};endEditCP(m_posgPoints);

	beginEditCP(m_posgColors);{
		m_posgColors->clear();
	};endEditCP(m_posgColors);

	beginEditCP(m_posgNormals);{
		m_posgNormals->clear();
	};endEditCP(m_posgNormals);

	beginEditCP(m_posgTexCoords);{
		m_posgTexCoords->clear();
	};endEditCP(m_posgTexCoords);

	vtkPolyData *pPolyData = NULL;
	if (dynamic_cast<vtkPolyDataMapper*>(this->GetMapper())){
		pPolyData = (vtkPolyData*) this->GetMapper()->GetInput();
		if (m_bVerbose){
			std::cerr << "		Using vtkPolyDataMapper directly" << std::endl;
		}
	} else if (dynamic_cast<vtkDataSetMapper*>(this->GetMapper())){
		vtkDataSetMapper *dataSetMapper = (vtkDataSetMapper*) this->GetMapper();
		pPolyData = (vtkPolyData*) dataSetMapper->GetPolyDataMapper()->GetInput();
		if (m_bVerbose){
			std::cerr << "		Using vtkPolyDataMapper via the vtkDataSetMapper" << std::endl;
		}
	}
	if (pPolyData == NULL) return NullFC;

	vtkCellArray *pCells;
	if (gl_primitive_type == GL_POINTS){
		pCells = pPolyData->GetVerts();
	} else if (gl_primitive_type == GL_LINE_STRIP){
		pCells = pPolyData->GetLines();
	} else if (gl_primitive_type == GL_POLYGON){
		pCells = pPolyData->GetPolys();
	} else if (gl_primitive_type == GL_TRIANGLE_STRIP){
		pCells = pPolyData->GetStrips();
	} else {
		std::cerr << "CVtkActorToOpenSG::ProcessGeometryNonIndexedCopyAttributes(int gl_primitive_type)" << std::endl;
		std::cerr << "	was called with non implemented gl_primitive_type!" << std::endl;
		return NullFC;
	}

	beginEditCP(m_posgTypes);
	beginEditCP(m_posgLengths);
	beginEditCP(m_posgPoints);
	beginEditCP(m_posgColors);
	beginEditCP(m_posgNormals);{
		int prim = 0;
		if (pCells->GetNumberOfCells() > 0){
			vtkIdType npts, *pts;
			for (pCells->InitTraversal(); pCells->GetNextCell(npts, pts); prim++){
				m_posgLengths->addValue(npts);
				m_posgTypes->addValue(GL_POLYGON);
				for (int i=0; i<npts; i++){
					double *aVertex;
					double *aNormal;
					unsigned char aColor[4];

					aVertex = pPolyData->GetPoint(pts[i]);
					m_posgPoints->addValue(Vec3f(aVertex[0], aVertex[1], aVertex[2]));

					if (m_iNormalType == PER_VERTEX){
						aNormal = m_pvtkNormals->GetTuple(pts[i]);
						m_posgNormals->addValue(Vec3f(aNormal[0], aNormal[1], aNormal[2]));
					} else if (m_iNormalType == PER_CELL){
						aNormal = m_pvtkNormals->GetTuple(prim);
						m_posgNormals->addValue(Vec3f(aNormal[0], aNormal[1], aNormal[2]));
					}

					if (m_iColorType == PER_VERTEX){
						m_pvtkColors->GetTupleValue(pts[i], aColor);
						float r = ((float) aColor[0]) / 255.0f;
						float g = ((float) aColor[1]) / 255.0f;
						float b = ((float) aColor[2]) / 255.0f;
						m_posgColors->addValue(Color3f(r, g, b));
					} else if (m_iColorType == PER_CELL){
						m_pvtkColors->GetTupleValue(prim, aColor);
						float r = ((float) aColor[0]) / 255.0f;
						float g = ((float) aColor[1]) / 255.0f;
						float b = ((float) aColor[2]) / 255.0f;
						m_posgColors->addValue(Color3f(r, g, b));
					}
				}
			}
		}
	};endEditCP(m_posgTypes);
	endEditCP(m_posgLengths);
	endEditCP(m_posgPoints);
	endEditCP(m_posgColors);
	endEditCP(m_posgNormals);

	//possibly getting the texture coordinates. These are always per vertex
	vtkPoints *points = pPolyData->GetPoints();
	if ((m_pvtkTexCoords != NULL) && (points != NULL)){
		int numPoints = points->GetNumberOfPoints();
		int numTexCoords = m_pvtkTexCoords->GetNumberOfTuples();
		if (numPoints == numTexCoords){
			beginEditCP(m_posgTexCoords);{
				int numTuples = m_pvtkTexCoords->GetNumberOfTuples();
				for (int i=0; i<numTuples; i++){
					double texCoords[3];
					m_pvtkTexCoords->GetTuple(i, texCoords);
					m_posgTexCoords->addValue(Vec2f(texCoords[0], texCoords[1]));
				}
			};endEditCP(m_posgTexCoords);
		}
	}

	ChunkMaterialPtr material = CreateMaterial();
	//GeometryPtr geo = Geometry::create();
	beginEditCP(m_posgGeometry);{
		m_posgGeometry->setPositions(m_posgPoints);
		m_posgGeometry->setTypes(m_posgTypes);
		m_posgGeometry->setLengths(m_posgLengths);
		m_posgGeometry->setMaterial(material);

		if (m_iNormalType != NOT_GIVEN) m_posgGeometry->setNormals(m_posgNormals);
		if (m_iColorType != NOT_GIVEN) m_posgGeometry->setColors(m_posgColors);
		if (m_posgTexCoords->getSize() > 0) m_posgGeometry->setTexCoords(m_posgTexCoords);
		//geo->setMaterial(getDefaultMaterial());
	}endEditCP(m_posgGeometry);

	if (m_bVerbose){
		std::cout << "		End ProcessGeometryNonIndexedCopyAttributes(int gl_primitive_type)" << std::endl;
	}
	return m_posgGeomNode;
}

NodePtr vtkOsgActor::GetNodePtr(){
	if (m_bVerbose){
		std::cerr << "Calling GetNodePtr()" << std::endl;
	}

	this->GetMapper()->Update();

	//Rendering with OpenSG simple indexed geometry
	if (((m_iNormalType == PER_VERTEX) || (m_iNormalType == NOT_GIVEN))  &&
		((m_iColorType == PER_VERTEX) || (m_iColorType == NOT_GIVEN))){
			return this->ProcessGeometryNormalsAndColorsPerVertex();
	}else{
		//Rendering with OpenSG non indexed geometry by copying a lot of attribute data
		if (m_iNumGLPolygons > 0){
			if (m_iNumGLPolygons != m_iNumGLPrimitives){
				std::cerr << "WARNING: vtkActor contains different kind of primitives" << std::endl;
			}
			return this->ProcessGeometryNonIndexedCopyAttributes(GL_POLYGON);
		} else if (m_iNumGLLineStrips > 0){
			if (m_iNumGLLineStrips != m_iNumGLPrimitives){
				std::cerr << "WARNING: vtkActor contains different kind of primitives" << std::endl;
			}
			return this->ProcessGeometryNonIndexedCopyAttributes(GL_LINE_STRIP);
		} else if (m_iNumGLTriStrips > 0){
			if (m_iNumGLTriStrips != m_iNumGLPrimitives){
				std::cerr << "WARNING: vtkActor contains different kind of primitives" << std::endl;
			}
			return this->ProcessGeometryNonIndexedCopyAttributes(GL_TRIANGLE_STRIP);
		} else if (m_iNumGLPoints > 0){
			if (m_iNumGLPoints != m_iNumGLPrimitives){
				std::cerr << "WARNING: vtkActor contains different kind of primitives" << std::endl;
			}
			return this->ProcessGeometryNonIndexedCopyAttributes(GL_POINTS);
		} else {
			return NullFC;
		}
	}
	return NullFC;
}
