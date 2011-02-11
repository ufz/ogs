/**
 * \file LineEditDialog.h
 * 2010/12/09 KR Initial implementation
 */

#ifndef LINEEDITDIALOG_H
#define LINEEDITDIALOG_H

#include <QtGui/QMainWindow>
#include "ui_LineEdit.h"

#include "PolylineVec.h"

class QStringListModel;

/**
 * \brief A dialog window for manipulation of polylines.
 * Currently included functionality is the concatenation of polylines 
 * as well as creating polygons or surfaces from polylines.
 */
class LineEditDialog : public QDialog, private Ui_LineEdit
{
	Q_OBJECT

public:
	LineEditDialog(const GEOLIB::PolylineVec &ply_vec, QDialog* parent = 0);
	~LineEditDialog(void);



private:
	std::vector<size_t> getSelectedIndeces(QStringList list);

	QStringListModel* _allPly;
	QStringListModel* _selPly;
	std::string _geoName;


private slots:
	/// Instructions when polylines are selected.
	void on_selectPlyButton_pressed();

	/// Instructions when polylines are deselected.
	void on_deselectPlyButton_pressed();

	/// Instructions if the OK-Button has been pressed.
	void accept();

	/// Instructions if the Cancel-Button has been pressed.
	void reject();

signals:
	void connectPolylines(const std::string&, std::vector<size_t>, bool, bool);	
	void triangulateSurface(const GEOLIB::Polyline);

};

#endif //LINEEDITDIALOG_H
