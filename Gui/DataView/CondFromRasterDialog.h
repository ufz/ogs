/**
 * \file CondFromRasterDialog.h
 * 2012/01/04 KR Initial implementation
 */

#ifndef CONDFROMRASTERDIALOG_H
#define CONDFROMRASTERDIALOG_H

#include "ui_CondFromRaster.h"
#include <QDialog>

#include "ProjectData.h"

namespace {
	class Mesh;
}

class StrictDoubleValidator;

/**
 * \brief A dialog window for creating DIRECT boundary conditions from raster files
 */
class CondFromRasterDialog : public QDialog, private Ui_CondFromRaster
{
	Q_OBJECT

public:
	CondFromRasterDialog(const std::vector<MeshLib::Mesh*> &msh_vec, QDialog* parent = 0);
	~CondFromRasterDialog(void);

private:
	const std::vector<MeshLib::Mesh*> _msh_vec;
	StrictDoubleValidator* _scale_validator;

private slots:
	void on_integrateButton_toggled(bool isSelected);
	void on_selectButton_pressed();

	/// Instructions if the OK-Button has been pressed.
	void accept();

	/// Instructions if the Cancel-Button has been pressed.
	void reject();

signals:
	void directNodesWritten(std::string);
	void transmitDisValues(std::vector< std::pair<size_t,double> >);
};

#endif //CONDFROMRASTERDIALOG_H
