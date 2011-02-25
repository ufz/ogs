/**
 * \file LineEditDialog.cpp
 * 2010/12/09 KR Initial implementation
 */

#include "LineEditDialog.h"
#include <QStringList>
#include <QStringListModel>


LineEditDialog::LineEditDialog(const GEOLIB::PolylineVec &ply_vec, QDialog* parent)
: QDialog(parent), _allPly(new QStringListModel), _selPly(new QStringListModel), _geoName(ply_vec.getName())
{
	setupUi(this);

	size_t nPly(ply_vec.size());
	QStringList list;
	for (size_t i=0; i<nPly; i++)
	{
		std::string ply_name("");
		ply_vec.getNameOfElementByID(i, ply_name);
		list.append("Line " + QString::number(i) + "  " + QString::fromStdString(ply_name));
	}
	_allPly->setStringList(list);

	this->allPlyView->setModel(_allPly);
	this->selectedPlyView->setModel(_selPly);
}

LineEditDialog::~LineEditDialog()
{
	delete _allPly;
	delete _selPly;
}

void LineEditDialog::on_selectPlyButton_pressed()
{
	QModelIndexList selected = this->allPlyView->selectionModel()->selectedIndexes();
	QStringList list = _selPly->stringList();

	for (QModelIndexList::iterator it = selected.begin(); it != selected.end(); ++it)
	{
		list.append(it->data().toString());

		_allPly->removeRow(it->row());
	}
	_selPly->setStringList(list);
}

void LineEditDialog::on_deselectPlyButton_pressed()
{
	QModelIndexList selected = this->selectedPlyView->selectionModel()->selectedIndexes();
	QStringList list = _allPly->stringList();

	for (QModelIndexList::iterator it = selected.begin(); it != selected.end(); ++it)
	{
		list.append(it->data().toString());

		_selPly->removeRow(it->row());
	}
	_allPly->setStringList(list);
}

void LineEditDialog::accept()
{
	std::vector<size_t> selectedIndeces = this->getSelectedIndeces(_selPly->stringList());
	std::string prox_string = this->proximityEdit->text().toStdString();
	double prox = (prox_string.empty()) ? 0.0 : strtod( prox_string.c_str(), 0 );
	emit connectPolylines(_geoName, selectedIndeces, prox, this->closePlyCheckBox->isChecked(), this->createSfcCheckBox->isChecked());
	this->done(QDialog::Accepted);
}

void LineEditDialog::reject()
{
	this->done(QDialog::Rejected);
}

std::vector<size_t> LineEditDialog::getSelectedIndeces(QStringList list)
{
	std::vector<size_t> indexList;
	for (QStringList::iterator it = list.begin(); it != list.end(); ++it)
	{
		QString s = it->mid(5, it->indexOf("  ")-5);
		indexList.push_back(atoi(s.toStdString().c_str()));
	}
	return indexList;
}
