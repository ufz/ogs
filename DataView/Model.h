/**
 * \file Model.h
 * 24/9/2009 LB Initial implementation
 * 05/05/2010 KR 2d graphic functionality removed and various layout changes
 */
#ifndef MODEL_H
#define MODEL_H

// ** INCLUDES **
#include <QVector>
#include <QItemSelection>

#include "TreeModel.h"

class ModelItem;
class QModelIndex;
class vtkAlgorithm;

/**
 * The Model is a base class for model implementation representing
 *  data. Data is stored in _data .
 */
class Model : public TreeModel
{
	Q_OBJECT

public:
	/// Constructor
	Model(QString name, QObject* parent = 0);

	/// Virtual empty destructor because this is a polymorphic base class
	virtual ~Model() {};

	/// Returns the number of rows
	//virtual int rowCount(const QModelIndex& parent = QModelIndex()) const = 0;

	/**
	 * Returns flags that the items are editable by default. Overwrite this
	 * function if want non editable items.
	 */
	Qt::ItemFlags flags(const QModelIndex& index) const;

	/// Removes count rows starting at row. Eventually this must be extended
	/// in a derived class to delete wrapped data, see PntsModel::removeRows()
	/// for an example.
	virtual bool removeRows(int row, int count, const QModelIndex & parent = QModelIndex() );

	/// Returns a ModelItem from an QModelIndex.
	//virtual ModelItem* itemFromIndex(const QModelIndex& index) const;

	/// Returns QModelIndexes from ModelItem.
	//QModelIndex indexFromItem(const ModelItem* item) const;

	/// Returns the dependent submodels
	QVector<Model*> subModels() const;

	/// Returns the name of the model.
	QString name() const { return _name; }

	/// Returns the Vtk polydata source object
	vtkAlgorithm* vtkSource() const { return _vtkSource; }

public slots:
	/// Deletes all data.
	void clearData();

	/// Deletes the selected data.
	void clearSelectedData();

	/// Reloads all data
	virtual void updateData();


protected:
	/// On a table model the model items are stored here.
	//QVector<ModelItem*> _data;

	/// Models which are dependent and created/managed by this model.
	QVector<Model*> _subModels;

	/// The actually selected items. This is used for removing the actually
	/// selected rows when calling clearSelectedData().
	QItemSelection _selectedItems;

	/// The name of the model. Usually this is the filename where the data
	/// comes from.
	QString _name;

	/// The Vtk data source object. This is the starting point for a Vtk data
	/// visualization pipeline.
	vtkAlgorithm* _vtkSource;

signals:

	
	/// Is emitted when the selection of the model is cleared due to removing items
	void selectionCleared();

};

#endif // MODEL_H
