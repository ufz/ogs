/**
 * \file VtkOGSFilter.cpp
 * 14/04/2010 KR Initial implementation
 *
 */


#include "VtkOGSFilter.h"

/* VTK Includes */
#include <vtkSmartPointer.h>
#include <vtkPolyDataAlgorithm.h>
#include <vtkImageAlgorithm.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridAlgorithm.h>
#include <vtkSphereSource.h>
#include <vtkCubeSource.h>
#include <vtkConeSource.h>
#include <vtkGlyph3D.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkThreshold.h>
#include <vtkTubeFilter.h>
#include <vtkPointData.h>
#include <vtkTransform.h>
#include <vtkImageMapToColors.h>
#include <vtkLookupTable.h>
#include <vtkCleanPolyData.h>
//test
#include <vtkPointData.h>
#include <vtkCellData.h>

/* Qt Includes */
#include <QImage>
#include <QPointF>
#include <QFileDialog>
#include <QFileInfo>

/* OGS Includes */
#include "OGSFilterInfo.h"
#include "OGSRaster.h"
#include "VtkColorByHeightFilter.h"
#include "VtkGeoImageSource.h"
#include "VtkMeshSource.h"
#include "VtkTextureOnSurfaceFilter.h"
#include "VtkImageDataToLinePolyDataFilter.h"
#include "VtkApplyColorTableFilter.h"

//test
#include <vtkDoubleArray.h>

VtkOGSFilter::VtkOGSFilter(vtkAlgorithm* algorithm,
		VtkOGSFilter::OGSVisFilter filter,
		TreeItem* parentItem, vtkPointSet* input,
		const QList<QVariant> data) 
: VtkVisPipelineItem(algorithm, parentItem, input, data)
{
	_algorithm = this->apply(algorithm, filter);
	//_parameterList["dummy"] = 0.0;
}


VtkOGSFilter::~VtkOGSFilter() 
{

}


std::vector<OGSFilterInfo> VtkOGSFilter::getAvailableFilters() 
{
	std::vector<OGSFilterInfo> availableFilters;
	availableFilters.push_back(OGSFilterInfo("Add Colormap to Image",		VtkOGSFilter::COLORMAPTOIMAGEFILTER,	OGSFilterInfo::IMAGEDATA));
	availableFilters.push_back(OGSFilterInfo("Apply Texture to Grid",	    VtkOGSFilter::TEXTOGRIDFILTER,		    OGSFilterInfo::UNSTRUCTUREDGRID));
	availableFilters.push_back(OGSFilterInfo("Apply Texture to Surface",    VtkOGSFilter::TEXTOSURFACEFILTER,		OGSFilterInfo::POLYDATA));
	availableFilters.push_back(OGSFilterInfo("Build Mesh from Image",		VtkOGSFilter::MESHFROMIMAGEFILTER,		OGSFilterInfo::IMAGEDATA));
	availableFilters.push_back(OGSFilterInfo("Convert Lines to Cylinders",  VtkOGSFilter::LINETOCYLINDERFILTER,     OGSFilterInfo::POLYDATA));
	availableFilters.push_back(OGSFilterInfo("Convert Points to Glyphs",    VtkOGSFilter::POINTTOGLYPHFILTER,       OGSFilterInfo::POLYDATA));
	availableFilters.push_back(OGSFilterInfo("Elevation-based Colour",      VtkOGSFilter::COLORBYHEIGHTFILTER_GRID, OGSFilterInfo::UNSTRUCTUREDGRID));
	availableFilters.push_back(OGSFilterInfo("Elevation-based Colour",      VtkOGSFilter::COLORBYHEIGHTFILTER_POLY, OGSFilterInfo::POLYDATA));
	availableFilters.push_back(OGSFilterInfo("Generate Surface",            VtkOGSFilter::SURFACEFILTER,            OGSFilterInfo::UNSTRUCTUREDGRID));
	availableFilters.push_back(OGSFilterInfo("Show Material Groups",	    VtkOGSFilter::MATGROUPFILTER,		    OGSFilterInfo::UNSTRUCTUREDGRID));
	availableFilters.push_back(OGSFilterInfo("Specify Thresholds",			VtkOGSFilter::THRESHOLDINGFILTER,		OGSFilterInfo::UNSTRUCTUREDGRID));
	availableFilters.push_back(OGSFilterInfo("Convert Image to Cylinders",	VtkOGSFilter::IMAGETOCYLINDERSFILTER, OGSFilterInfo::IMAGEDATA));
	
	return availableFilters;
}

vtkAlgorithm* VtkOGSFilter::apply(vtkAlgorithm* input, OGSVisFilter filter)
{
	switch (filter)
	{
	case POINTTOGLYPHFILTER:
		return static_cast<OGSPoint2GlyphFilter*>(this)->applyFilter(static_cast<vtkPolyDataAlgorithm*>(input));
	case LINETOCYLINDERFILTER:
		return static_cast<OGSLine2CylinderFilter*>(this)->applyFilter(static_cast<vtkPolyDataAlgorithm*>(input));
	case SURFACEFILTER:
		return static_cast<OGSSurfaceFilter*>(this)->applyFilter(static_cast<vtkUnstructuredGridAlgorithm*>(input));
	case COLORBYHEIGHTFILTER_GRID:
		return static_cast<OGSColorByHeightFilter*>(this)->applyFilter(static_cast<vtkUnstructuredGridAlgorithm*>(input));
	case COLORBYHEIGHTFILTER_POLY:
		return static_cast<OGSColorByHeightFilter*>(this)->applyFilter(static_cast<vtkPolyDataAlgorithm*>(input));
	case MATGROUPFILTER:
		return static_cast<OGSMaterialGroupFilter*>(this)->applyFilter(static_cast<vtkUnstructuredGridAlgorithm*>(input));
	case TEXTOGRIDFILTER:
		return static_cast<OGSTextureToSurfaceFilter*>(this)->applyFilter(static_cast<vtkUnstructuredGridAlgorithm*>(input));
	case TEXTOSURFACEFILTER:
		return static_cast<OGSTextureToSurfaceFilter*>(this)->applyFilter(static_cast<vtkPolyDataAlgorithm*>(input));
	case THRESHOLDINGFILTER:
		return static_cast<OGSThresholdFilter*>(this)->applyFilter(static_cast<vtkUnstructuredGridAlgorithm*>(input));
	case MESHFROMIMAGEFILTER:
		return static_cast<OGSMeshFromImageFilter*>(this)->applyFilter(static_cast<vtkImageAlgorithm*>(input));
	case COLORMAPTOIMAGEFILTER:
		return static_cast<OGSColormapToImageFilter*>(this)->applyFilter(static_cast<vtkImageAlgorithm*>(input));
	case IMAGETOCYLINDERSFILTER:
		return static_cast<OGSImageToCylindersFilter*>(this)->applyFilter(static_cast<vtkImageAlgorithm*>(input));
	default:
		return input;
	}
}

std::map<std::string, double> VtkOGSFilter::getParameters()
{
	return _parameterList;
}

void VtkOGSFilter::setParameter( std::string parameterName, double value )
{
	_parameterList[parameterName]=value;
}

vtkPolyDataAlgorithm* OGSPoint2GlyphFilter::applyFilter(vtkPolyDataAlgorithm* algorithm)
{
	double radius = _parameterList["Radius"];
	vtkSmartPointer<vtkSphereSource> sphere = vtkSmartPointer<vtkSphereSource>::New();
		sphere->SetRadius(250);
		sphere->SetPhiResolution(10);
		sphere->SetThetaResolution(10);
		sphere->SetReleaseDataFlag(1);

	vtkGlyph3D* glyphs = vtkGlyph3D::New();
		glyphs->ScalingOn();
		glyphs->SetScaleModeToScaleByScalar();
		glyphs->SetScaleModeToDataScalingOff();
		glyphs->SetSource(sphere->GetOutput());
		glyphs->SetInput(algorithm->GetOutput());

	return glyphs;
}


vtkPolyDataAlgorithm* OGSLine2CylinderFilter::applyFilter(vtkPolyDataAlgorithm* algorithm)
{
	double radius = this->_parameterList["Radius"];

	// collapse coincident points
	vtkSmartPointer<vtkCleanPolyData> mergePoints = vtkSmartPointer<vtkCleanPolyData>::New();
		mergePoints->SetInputConnection(0, algorithm->GetOutputPort(0));
		mergePoints->SetTolerance(0.0);
		mergePoints->ConvertLinesToPointsOn();

/*
	vtkSmartPointer<vtkUnsignedCharArray> colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
		colors->SetNumberOfComponents(3);
		colors->SetName("StratColors");

	int nPoints = algorithm->GetOutput()->GetNumberOfPoints();
	int nLines = algorithm->GetOutput()->GetNumberOfCells();

	vtkUnsignedCharArray* orgColors = vtkUnsignedCharArray::SafeDownCast(algorithm->GetOutput()->GetCellData()->GetScalars());
	
	int start = nPoints-nLines;
	for (int i=644; i<nPoints; i++)
	{
		unsigned char c[3];
		orgColors->GetTupleValue(i, c);
		colors->InsertNextTupleValue(c);
	}
*/
	vtkTubeFilter* tubes = vtkTubeFilter::New();
		tubes->SetInputConnection(0, mergePoints->GetOutputPort(0));
		//tubes->GetOutput()->GetCellData()->RemoveArray("StationColors");
		//tubes->GetOutput()->GetCellData()->AddArray(colors);
		//tubes->GetOutput()->GetCellData()->SetActiveScalars("StratColors");
		tubes->SetInputArrayToProcess(1,0,0,vtkDataObject::FIELD_ASSOCIATION_CELLS,"StratColors");
		tubes->SetRadius(150/*2500*/);
		tubes->SetNumberOfSides(10);
	return tubes;
}

vtkUnstructuredGridAlgorithm* OGSMaterialGroupFilter::applyFilter(vtkUnstructuredGridAlgorithm* algorithm)
{
	vtkSmartPointer<VtkMeshSource> meshsource = static_cast<VtkMeshSource*>(algorithm);
		meshsource->setColorsFromMaterials();
		meshsource->SetScalarVisibility(true);

	return meshsource;
}

vtkPolyDataAlgorithm* OGSSurfaceFilter::applyFilter(vtkUnstructuredGridAlgorithm* algorithm)
{
	vtkDataSetSurfaceFilter* surfaceFilter = vtkDataSetSurfaceFilter::New();
		surfaceFilter->SetInputConnection(0, algorithm->GetOutputPort(0));	
		surfaceFilter->Update();
	return surfaceFilter;
}

vtkPolyDataAlgorithm* OGSColorByHeightFilter::applyFilter(vtkUnstructuredGridAlgorithm* algorithm)
{
	VtkMeshSource* meshsource = static_cast<VtkMeshSource*>(algorithm);
		meshsource->SetScalarVisibility(true);
	vtkDataSetSurfaceFilter* surfaceFilter = vtkDataSetSurfaceFilter::New();
		surfaceFilter->SetInputConnection(0, meshsource->GetOutputPort(0));	
		surfaceFilter->Update();
	return this->applyFilter(static_cast<vtkPolyDataAlgorithm*>(surfaceFilter));
}

vtkPolyDataAlgorithm* OGSColorByHeightFilter::applyFilter(vtkPolyDataAlgorithm* algorithm)
{
	VtkColorByHeightFilter* heightFilter = VtkColorByHeightFilter::New();
		heightFilter->SetInputConnection(0, algorithm->GetOutputPort(0));
		heightFilter->GetColorLookupTable()->SetTableRange(-35, 800); // default min- and max-height (default: blue to red)
		//heightFilter->GetColorLookupTable()->setInterpolationType(ColorLookupTable::EXPONENTIAL);
	unsigned char a[4] = { 0, 0, 255, 255 };
	unsigned char b[4] = { 0, 255, 0, 255 };
	unsigned char c[4] = { 255, 255, 0, 255 };
	unsigned char d[4] = { 255, 0, 0, 255 };
	unsigned char e[4] = { 255, 255, 255, 255 };
		heightFilter->GetColorLookupTable()->setColor(0.0, a);
		heightFilter->GetColorLookupTable()->setColor(0.2, b); // green at about 150m
		heightFilter->GetColorLookupTable()->setColor(0.6, c); // yellow at about 450m and changing to red from then on
		heightFilter->GetColorLookupTable()->setColor(1.0, d);
		//heightFilter->GetColorLookupTable()->setColor(1.0, e);
		heightFilter->Update(); 

	return heightFilter;
}

vtkPolyDataAlgorithm* OGSTextureToSurfaceFilter::applyFilter(vtkUnstructuredGridAlgorithm* algorithm)
{
	vtkDataSetSurfaceFilter* surfaceFilter = vtkDataSetSurfaceFilter::New();
		surfaceFilter->SetInputConnection(0, algorithm->GetOutputPort(0));	
		surfaceFilter->Update();
	return this->applyFilter(static_cast<vtkPolyDataAlgorithm*>(surfaceFilter));
}

vtkPolyDataAlgorithm* OGSTextureToSurfaceFilter::applyFilter(vtkPolyDataAlgorithm* algorithm)
{
	QWidget* parent = 0;
	QString fileName = QFileDialog::getOpenFileName(parent, "Select raster file to apply as texture", "","Raster files (*.asc *.bmp *.jpg *.png *.tif);;");
	QFileInfo fi(fileName);

	if ((fi.suffix().toLower() == "asc") || (fi.suffix().toLower() == "tif") || (fi.suffix().toLower() == "png") || 
		(fi.suffix().toLower() == "jpg") || (fi.suffix().toLower() == "bmp"))
	{
		QImage img;
		QPointF origin;
		double scalingFactor=0;

		OGSRaster::loadImage(fileName, img, origin, scalingFactor);
		std::pair<float, float> org(origin.x(), origin.y()); 
		VtkTextureOnSurfaceFilter* surface = VtkTextureOnSurfaceFilter::New();
			surface->SetRaster(img, org, scalingFactor);
			surface->SetInputConnection(0, algorithm->GetOutputPort(0));
			surface->Update();
		return surface;
		
	}
	return NULL;
}

vtkUnstructuredGridAlgorithm* OGSThresholdFilter::applyFilter(vtkUnstructuredGridAlgorithm* algorithm)
{
	size_t min = static_cast<size_t>(_parameterList["Minimum"]);
	size_t max = static_cast<size_t>(_parameterList["Maximum"]);
	vtkThreshold* threshold = vtkThreshold::New();
		//threshold->SetAttributeModeToUseCellData();
		//threshold->SetComponentModeToUseSelected();
		threshold->SetInputConnection(0, algorithm->GetOutputPort(0));
	//	threshold->ThresholdByLower(1);
		threshold->ThresholdBetween(4,4);
		threshold->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, static_cast<VtkMeshSource*>(algorithm)->GetMaterialArrayName());
		threshold->Update();
	
	return threshold;
}

vtkUnstructuredGridAlgorithm* OGSMeshFromImageFilter::applyFilter(vtkImageAlgorithm* algorithm)
{
	vtkSmartPointer<VtkGeoImageSource> imageSource = vtkSmartPointer<VtkGeoImageSource>(VtkGeoImageSource::SafeDownCast(algorithm));
	vtkSmartPointer<vtkImageData> image = vtkSmartPointer<vtkImageData>(imageSource->GetOutput());

	Mesh_Group::CFEMesh* mesh = GridAdapter::convertImgToMesh(image, imageSource->getOrigin(), imageSource->getSpacing());
	GridAdapter* grid = new GridAdapter(mesh);
	VtkMeshSource* meshSource = VtkMeshSource::New();
	meshSource->SetGrid(grid);
	delete mesh;
	return meshSource;
}

vtkImageAlgorithm* OGSColormapToImageFilter::applyFilter(vtkImageAlgorithm* algorithm)
{
	vtkSmartPointer<vtkLookupTable> colormap = vtkSmartPointer<vtkLookupTable>::New();
		colormap->SetTableRange(0, 100);
		colormap->SetHueRange(0.0, 0.666);
		colormap->ForceBuild();

	vtkImageMapToColors* map = vtkImageMapToColors::New();
		map->SetInputConnection(0, algorithm->GetOutputPort(0));
		map->SetLookupTable(colormap);
		map->SetPassAlphaToOutput(1);
		map->Update();

	return map;
}

OGSImageToCylindersFilter::OGSImageToCylindersFilter(
	vtkAlgorithm* algorithm, VtkOGSFilter::OGSVisFilter filter,
	TreeItem* parentItem, vtkPointSet* input,
	const QList<QVariant> data /*= QList<QVariant>()*/)
	: VtkOGSFilter(algorithm, filter, parentItem, input, data)
{
}

vtkPolyDataAlgorithm* OGSImageToCylindersFilter::applyFilter(vtkImageAlgorithm* algorithm)
{
	vtkSmartPointer<VtkImageDataToLinePolyDataFilter> lineFilter = vtkSmartPointer<VtkImageDataToLinePolyDataFilter>::New();
	lineFilter->SetInputConnection(algorithm->GetOutputPort());
	lineFilter->SetLengthScaleFactor(30); // needs to be a parameter
	lineFilter->Update();

	// Build colour map
	vtkSmartPointer<vtkLookupTable> colormap = vtkSmartPointer<vtkLookupTable>::New();
		colormap->SetTableRange(20, 200);
		//colormap->SetTableRange(0, 100);
		//colormap->SetTableRange(1, 255 * lineFilter->GetLengthScaleFactor());
		colormap->SetHueRange(0.0, 0.666);
		colormap->ForceBuild();


	VtkApplyColorTableFilter* ctf = VtkApplyColorTableFilter::New();
		ctf->SetInputConnection(0, lineFilter->GetOutputPort());
		ctf->SetColorLookupTable(colormap);
		ctf->Update();

/*
	vtkSmartPointer<vtkUnsignedCharArray> colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
		colors->SetNumberOfComponents(3);
		colors->SetName("Colors");

	// Read the active scalar array of the linefilter and use it to define colours for the resulting tubes
	double color[4];
	vtkSmartPointer<vtkUnsignedCharArray> scalars = vtkUnsignedCharArray::SafeDownCast(lineFilter->GetOutput()->GetPointData()->GetScalars());
	int nPoints = lineFilter->GetOutput()->GetNumberOfPoints();
	for (int i=1; i<nPoints; i+=2)
	{
		double* value = scalars->GetTuple(i);
		colormap->GetColor(value[0], color);
		unsigned char tubecolor[3] = { color[0]*255, color[1]*255, color[2]*255 };
		// set colour for start- and end-point of the line
		colors->InsertNextTupleValue(tubecolor); 
		colors->InsertNextTupleValue(tubecolor);
	}

	lineFilter->GetOutput()->GetPointData()->AddArray(colors);
	lineFilter->GetOutput()->GetPointData()->SetActiveScalars("Colors");
	lineFilter->Update();
*/

	vtkTubeFilter* tubeFilter = vtkTubeFilter::New();
	tubeFilter->SetInputConnection(ctf->GetOutputPort());
	tubeFilter->CappingOn();
	tubeFilter->SetNumberOfSides(6);
	tubeFilter->SetRadius(lineFilter->GetImageSpacing() * 0.25);

	return tubeFilter;

}
