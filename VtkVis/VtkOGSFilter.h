/**
 * \file VtkOGSFilter.h
 * 14/04/2010 KR Initial implementation
 *
 */

#ifndef VTKOGSFILTER_H
#define VTKOGSFILTER_H

#include <map>
#include <string>
#include <vector>
#include <cstddef>

#include "VtkVisPipelineItem.h"

class OGSFilterInfo;
class vtkImageAlgorithm;
class vtkPolyDataAlgorithm;
class vtkUnstructuredGridAlgorithm;
class vtkImageAlgorithm;

class OGSFilterInfo;

/**
 * \brief Base class providing access functions for use of VTK filters in OGS.
 *
 * Specific filters are derived from this class and must implement the applyFilter()-function. Existing filter should
 * registered in the OGSVisFilter Enum and the getAvailableFilters()-method to be available in the GUI filter dialog.
 */
class VtkOGSFilter : public VtkVisPipelineItem
{
public:
	/// Registered filter aliases
	enum OGSVisFilter
	{
		INVALIDFILTER = 0,
		POINTTOGLYPHFILTER,       
		LINETOCYLINDERFILTER,     
		SURFACEFILTER,            
		COLORBYHEIGHTFILTER_GRID,
		COLORBYHEIGHTFILTER_POLY, 
		MATGROUPFILTER,
		TEXTOGRIDFILTER,          
		TEXTOSURFACEFILTER,       
		THRESHOLDINGFILTER,
		MESHFROMIMAGEFILTER,
		COLORMAPTOIMAGEFILTER,
		IMAGETOCYLINDERSFILTER
	};

	/// Constructor
	VtkOGSFilter(vtkAlgorithm* algorithm,
		VtkOGSFilter::OGSVisFilter filter,
		TreeItem* parentItem, vtkPointSet* input,
		const QList<QVariant> data = QList<QVariant>());

	~VtkOGSFilter();

	/// Implementation of the filter or filter pipeline (must be implemented in derived classes!)
	virtual vtkAlgorithm* applyFilter(vtkAlgorithm* input) { return NULL; };

	/// Get all parameters from the parameter list.
	std::map<std::string, double> getParameters();

	/// Reset the value of the parameter parameterName.
	void setParameter(std::string parameterName, double value);

	/// Returns all registered filters as VtkFilterItem-Objects.
	static std::vector<OGSFilterInfo> getAvailableFilters();

protected:
/*
	int RequestData(vtkInformation* request, 
		            vtkInformationVector** inputVector, 
					vtkInformationVector* outputVector);
*/
	std::map<std::string, double> _parameterList;

private:
	/// Access function for applying filters via the GUI
	vtkAlgorithm* apply(vtkAlgorithm* input, VtkOGSFilter::OGSVisFilter filter);

	//std::vector<VtkFilterItem> _availableFilters;

};

class OGSColorByHeightFilter : public VtkOGSFilter
{
public:
	OGSColorByHeightFilter(vtkAlgorithm* algorithm,
		VtkOGSFilter::OGSVisFilter filter,
		TreeItem* parentItem, vtkPointSet* input,
		const QList<QVariant> data = QList<QVariant>()) 
	: VtkOGSFilter(algorithm, filter, parentItem, input, data) {};

	/**
	 * \brief Elevation filter for colouring objects based on their height. 
	 * 
	 * This function uses the ColorLookupTable class and allows customisaton of the applied transfer function.
	 */
	vtkPolyDataAlgorithm* applyFilter(vtkPolyDataAlgorithm* algorithm);

	/**
	 * \brief Elevation filter for colouring objects based on their height. 
	 * 
	 * This function uses the ColorLookupTable class and allows customisaton of the applied transfer function.
	 */
	vtkPolyDataAlgorithm* applyFilter(vtkUnstructuredGridAlgorithm* algorithm);
};

class OGSMaterialGroupFilter : public VtkOGSFilter
{
public:
	OGSMaterialGroupFilter(vtkAlgorithm* algorithm,
		VtkOGSFilter::OGSVisFilter filter,
		TreeItem* parentItem, vtkPointSet* input,
		const QList<QVariant> data = QList<QVariant>()) 
		: VtkOGSFilter(algorithm, filter, parentItem, input, data) {};

	/// Colorises the mesh based on its material groups.
	vtkUnstructuredGridAlgorithm* applyFilter(vtkUnstructuredGridAlgorithm* algorithm);
};

class OGSPoint2GlyphFilter : public VtkOGSFilter
{
public:
	OGSPoint2GlyphFilter(vtkAlgorithm* algorithm,
		VtkOGSFilter::OGSVisFilter filter,
		TreeItem* parentItem, vtkPointSet* input,
		const QList<QVariant> data = QList<QVariant>()) 
	: VtkOGSFilter(algorithm, filter, parentItem, input, data)
	{ _parameterList.insert(std::pair<std::string, double>("Radius", 150)); };

	/// Filter for generating spheres from points.
	vtkPolyDataAlgorithm* applyFilter(vtkPolyDataAlgorithm* algorithm);
};

class OGSLine2CylinderFilter : public VtkOGSFilter
{
public:
	OGSLine2CylinderFilter(vtkAlgorithm* algorithm,
		VtkOGSFilter::OGSVisFilter filter,
		TreeItem* parentItem, vtkPointSet* input,
		const QList<QVariant> data = QList<QVariant>()) 
	: VtkOGSFilter(algorithm, filter, parentItem, input, data) 
	{ _parameterList.insert(std::pair<std::string, double>("Radius", 150)); };

	/// Filter for generating spheres from points.
	vtkPolyDataAlgorithm* applyFilter(vtkPolyDataAlgorithm* algorithm);
};

class OGSSurfaceFilter : public VtkOGSFilter
{
public:
	OGSSurfaceFilter(vtkAlgorithm* algorithm,
		VtkOGSFilter::OGSVisFilter filter,
		TreeItem* parentItem, vtkPointSet* input,
		const QList<QVariant> data = QList<QVariant>()) 
	: VtkOGSFilter(algorithm, filter, parentItem, input, data) {};

	/// Filter for generating a surface (i.e. a polydata-object) from a mesh (i.e. a grid-object)
	/// This is a base-filter for various other filters that require the input to be of type PolyData.
	vtkPolyDataAlgorithm* applyFilter(vtkUnstructuredGridAlgorithm* algorithm);
};

class OGSTextureToSurfaceFilter : public VtkOGSFilter
{
public:
	OGSTextureToSurfaceFilter(vtkAlgorithm* algorithm,
		VtkOGSFilter::OGSVisFilter filter,
		TreeItem* parentItem, vtkPointSet* input,
		const QList<QVariant> data = QList<QVariant>()) 
	: VtkOGSFilter(algorithm, filter, parentItem, input, data) {};

	/// Maps a texture onto a grid.
	vtkPolyDataAlgorithm* applyFilter(vtkUnstructuredGridAlgorithm* algorithm);

	/// Maps a texture onto a surface.
	vtkPolyDataAlgorithm* applyFilter(vtkPolyDataAlgorithm* algorithm);
};

class OGSThresholdFilter : public VtkOGSFilter
{
public:
	OGSThresholdFilter(vtkAlgorithm* algorithm,
		VtkOGSFilter::OGSVisFilter filter,
		TreeItem* parentItem, vtkPointSet* input,
		const QList<QVariant> data = QList<QVariant>()) 
	: VtkOGSFilter(algorithm, filter, parentItem, input, data)
	{		
		_parameterList.insert(std::pair<std::string, double>("Minimum", 1));
		_parameterList.insert(std::pair<std::string, double>("Maximum", 2));
	};

	/// Filtering data based on given thresholds.
	vtkUnstructuredGridAlgorithm* applyFilter(vtkUnstructuredGridAlgorithm* algorithm);
};

class OGSMeshFromImageFilter : public VtkOGSFilter
{
public:
	OGSMeshFromImageFilter(vtkAlgorithm* algorithm,
		VtkOGSFilter::OGSVisFilter filter,
		TreeItem* parentItem, vtkPointSet* input,
		const QList<QVariant> data = QList<QVariant>()) 
	: VtkOGSFilter(algorithm, filter, parentItem, input, data){};

	/// Filtering data based on given thresholds. 
	vtkUnstructuredGridAlgorithm* applyFilter(vtkImageAlgorithm* algorithm);
};

class OGSColormapToImageFilter : public VtkOGSFilter
{
public:
	OGSColormapToImageFilter(vtkAlgorithm* algorithm,
		VtkOGSFilter::OGSVisFilter filter,
		TreeItem* parentItem, vtkPointSet* input,
		const QList<QVariant> data = QList<QVariant>()) 
	: VtkOGSFilter(algorithm, filter, parentItem, input, data){};

	/// Filtering data based on given thresholds.
	vtkImageAlgorithm* applyFilter(vtkImageAlgorithm* algorithm);
};

class OGSImageToCylindersFilter : public VtkOGSFilter
{
public:
	OGSImageToCylindersFilter(
		vtkAlgorithm* algorithm, VtkOGSFilter::OGSVisFilter filter,
		TreeItem* parentItem, vtkPointSet* input,
		const QList<QVariant> data = QList<QVariant>());

	vtkPolyDataAlgorithm* applyFilter(vtkImageAlgorithm* algorithm);
};

#endif // VTKOGSFILTER_H
