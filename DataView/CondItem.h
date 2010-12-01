/**
 * \file CondItem.h
 * 20/10/2010 KR Initial implementation
 */

#ifndef CONDITEM_H
#define CONDITEM_H

#include "TreeItem.h"
#include "VtkPointsSource.h"


/**
 * \brief A TreeItem containing conditions of a FEM as well as the vtk object representing these conditions.
 * \sa TreeItem
 */
class CondItem : public TreeItem
{

public:
	/// Constructor, automatically generates VTK object
	CondItem(const QList<QVariant> &data, TreeItem *parent, const std::vector<CSourceTerm*> *sourceterms) : TreeItem(data, parent)
	{
		(void)sourceterms;
		_vtkSource = VtkPointsSource::New();
	};

	~CondItem()	{ _vtkSource->Delete();	};

	/// Returns the VTK object.
	VtkPointsSource* vtkSource() const { return _vtkSource; };	


private:
	VtkPointsSource* _vtkSource;
};

#endif //CONDITEM_H
