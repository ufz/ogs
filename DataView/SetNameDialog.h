/**
 * \file SetNameDialog.h
 * 2011/10/26 KR Initial implementation
 */

#ifndef SETNAMEDIALOG_H
#define SETNAMEDIALOG_H

#include "GeoType.h"

#include <QDialog>

class QDialogButtonBox;
class QLabel;
class QLineEdit;
class QVBoxLayout;

/**
 * \brief Small dialog for setting a name for an object.
 */
class SetNameDialog : public QDialog
{
	Q_OBJECT

public:
	/// Constructor
	SetNameDialog(const std::string &parent_name, 
				  const std::string &object_type_name, 
				  size_t id, 
				  const std::string &old_name, 
				  QDialog* parent = 0);
	~SetNameDialog();

	QDialogButtonBox* _buttonBox; /// The buttons used in this dialog.

private:
	/// Constructs a dialog window
	void setupDialog(const std::string &old_name);

	QLabel* _txt_label;
	QLineEdit* _new_name;
	QVBoxLayout* _layout;

	std::string _parent_name;
	std::string _object_type_name;
	size_t _id;

private slots:
	/// Instructions if the OK-Button has been pressed.
	void accept();

	/// Instructions if the Cancel-Button has been pressed.
	void reject();

signals:
	void requestNameChange(const std::string&, const GEOLIB::GEOTYPE, size_t, std::string);
};

#endif //SETNAMEDIALOG_H
