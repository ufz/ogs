/**
 * \file MshEditDialog.h
 * 2010/11/09 KR Initial implementation
 */

#ifndef MSHEDITDIALOG_H
#define MSHEDITDIALOG_H

#include "ui_MshEdit.h"
#include <QDialog>

#include "MshLayerMapper.h"

class QPushButton;
class QCheckBox;

namespace MeshLib
{
class Mesh;
}

/**
 * \brief A dialog window for editing meshes in various ways
 */
class MshEditDialog : public QDialog, private Ui_MshEdit
{
	Q_OBJECT

public:
	MshEditDialog(const MeshLib::Mesh* mesh, QDialog* parent = 0);
	~MshEditDialog(void);

private:
	const MeshLib::Mesh* _msh;
	QVector<QLabel*> _labels;
	QMap<QPushButton*, QLineEdit*> _fileButtonMap;
	QVector<QLineEdit*> _edits;
	QVector<QPushButton*> _buttons;
	QCheckBox* _noDataDeleteBox;

private slots:
	void getFileName();

	/// Instructions if the OK-Button has been pressed.
	void accept();

	/// Instructions if the Cancel-Button has been pressed.
	void reject();

signals:
	void mshEditFinished(MeshLib::Mesh*, std::string&);
};

#endif //MSHEDITDIALOG_H
