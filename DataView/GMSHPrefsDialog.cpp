/**
 * \file GMSHPresDialog.cpp
 * 2010/01/21 KR Initial implementation
 */

#include "GeoObjects.h"
#include "GMSHPrefsDialog.h"
#include <QStringList>
#include <QStringListModel>


GMSHPrefsDialog::GMSHPrefsDialog(const GEOLIB::GEOObjects* geoObjects, QDialog* parent) 
: QDialog(parent), _allGeo(new QStringListModel), _selGeo(new QStringListModel)
{
	setupUi(this);	

	// default parameters
	this->param1->setText("0");
	this->param2->setText("0");
	this->param3->setText("0");
	this->param4->setText("0");

	std::vector<std::string> geoNames;
	geoObjects->getGeometryNames(geoNames);
	size_t nGeoObjects(geoNames.size());
	QStringList list;

	for (size_t i=0; i<nGeoObjects; i++)
	{
		list.append(QString::fromStdString(geoNames[i]));
	}

	if (list.empty())
	{
		this->selectGeoButton->setDisabled(true);
		this->deselectGeoButton->setDisabled(true);
		list.append("(No geometry available.)");
	}
	_allGeo->setStringList(list);
	this->allGeoView->setModel(_allGeo);
	this->selectedGeoView->setModel(_selGeo);
	this->radioAdaptive->toggle();	// default is adaptive meshing
	this->on_radioAdaptive_toggled(true);
}

GMSHPrefsDialog::~GMSHPrefsDialog()
{
	delete _allGeo;
	delete _selGeo;
}

void GMSHPrefsDialog::on_selectGeoButton_pressed()
{
	QModelIndexList selected = this->allGeoView->selectionModel()->selectedIndexes();
	QStringList list = _selGeo->stringList();

	for (QModelIndexList::iterator it = selected.begin(); it != selected.end(); ++it)
	{
		list.append(it->data().toString());

		_allGeo->removeRow(it->row());
	}
	_selGeo->setStringList(list);
}

void GMSHPrefsDialog::on_deselectGeoButton_pressed()
{
	QModelIndexList selected = this->selectedGeoView->selectionModel()->selectedIndexes();
	QStringList list = _allGeo->stringList();

	for (QModelIndexList::iterator it = selected.begin(); it != selected.end(); ++it)
	{
		list.append(it->data().toString());

		_selGeo->removeRow(it->row());
	}
	_allGeo->setStringList(list);
}

void GMSHPrefsDialog::on_radioAdaptive_toggled(bool isTrue)
{
	if (isTrue) // meshing set to adaptive
	{
		this->param1->setEnabled(true);
		this->param2->setEnabled(true);
		this->param3->setEnabled(true);
		this->param4->setEnabled(false);
	}
	else // meshing set to homogeneous
	{
		this->param1->setEnabled(false);
		this->param2->setEnabled(false);
		this->param3->setEnabled(false);
		this->param4->setEnabled(true);
	}
}

void GMSHPrefsDialog::accept()
{
	std::vector<std::string> selectedObjects = this->getSelectedObjects(_selGeo->stringList());
	double val1(-1), val2(-1), val3(-1), val4(-1);
	
	if (this->radioAdaptive->isChecked())
	{
		val1 = strtod(param1->text().toStdString().c_str(), 0);
		val2 = strtod(param2->text().toStdString().c_str(), 0);
		val3 = strtod(param3->text().toStdString().c_str(), 0);
	}
	else
	{
		val4 = strtod(param4->text().toStdString().c_str(), 0);
	}

	emit requestMeshing(selectedObjects, val1, val2, val3, val4);
	this->done(QDialog::Accepted);
}

void GMSHPrefsDialog::reject()
{
	this->done(QDialog::Rejected);
}

std::vector<std::string> GMSHPrefsDialog::getSelectedObjects(QStringList list)
{
	std::vector<std::string> indexList;
	for (QStringList::iterator it = list.begin(); it != list.end(); ++it)
	{
		indexList.push_back(it->toStdString());
	}
	return indexList;
}