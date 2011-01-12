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
 * \brief A dialog window for editing meshes in various ways
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


private slots:
	void on_selectPlyButton_pressed();

	void on_deselectPlyButton_pressed();

	/// Instructions if the OK-Button has been pressed.
	void accept();

	/// Instructions if the Cancel-Button has been pressed.
	void reject();

signals:
	void connectPolylines(std::vector<size_t>, bool, bool);	

};

#endif //LINEEDITDIALOG_H
