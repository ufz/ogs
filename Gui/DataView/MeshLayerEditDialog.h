/**
 * \file
 * \author Karsten Rink
 * \date   2010-11-09
 * \brief  Definition of the MeshLayerEditDialog class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MESHLAYEREDITDIALOG_H
#define MESHLAYEREDITDIALOG_H

#include "ui_MeshLayerEdit.h"
#include <QDialog>

#include "MshLayerMapper.h"

class QPushButton;
class QCheckBox;
class QGroupBox;
class QGridLayout;
class QVBoxLayout;
class QRadioButton;

namespace MeshLib
{
class Mesh;
}

/**
 * \brief A dialog window for editing meshes in various ways
 */
class MeshLayerEditDialog : public QDialog, private Ui_MeshLayerEdit
{
	Q_OBJECT

public:
	MeshLayerEditDialog(const MeshLib::Mesh* mesh, QDialog* parent = 0);
	~MeshLayerEditDialog(void);

private:
	const MeshLib::Mesh* _msh;
	QVector<QLabel*> _labels;
	QMap<QPushButton*, QLineEdit*> _fileButtonMap;
	QVector<QLineEdit*> _edits;
	QVector<QPushButton*> _buttons;
	QCheckBox* _noDataDeleteBox;

	QLabel* _nLayerLabel;
	QLabel* _nLayerExplanation;
	QLineEdit* _layerEdit;
	QPushButton* _nextButton;
	QGroupBox* _layerBox;
	QGroupBox* _radioButtonBox;
	QGridLayout* _layerSelectionLayout;
	QVBoxLayout* _radiobuttonLayout;
	QRadioButton* _selectButton1;
	QRadioButton* _selectButton2;
	unsigned _n_layers;
	bool _use_rasters;

private slots:
	void getFileName();

	void nextButtonPressed();

	void createWithRasters();

	void createStatic();

	/// Instructions if the OK-Button has been pressed.
	void accept();

	/// Instructions if the Cancel-Button has been pressed.
	void reject();

signals:
	void mshEditFinished(MeshLib::Mesh*);
};

#endif //MESHLAYEREDITDIALOG_H
