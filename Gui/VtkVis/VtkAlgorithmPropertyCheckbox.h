/**
 * \file VtkAlgorithmPropertyCheckbox.h
 * 20/10/2010 LB Initial implementation
 */

#ifndef VTKALGORITHMPROPERTIESCHECKBOX_H
#define VTKALGORITHMPROPERTIESCHECKBOX_H

#include <QCheckBox>

class VtkAlgorithmProperties;

/// @brief This checkbox sets a user property on the given VtkAlgorithmProperties
/// object automatically.
class VtkAlgorithmPropertyCheckbox : public QCheckBox
{
	Q_OBJECT

public:
	/// @brief Constructor.
	/// @param value The initial check state.
	/// @param name The name of the user property to set.
	/// @param algProps The VtkAlgorithmProperties object.
	/// @param parent The parent widget.
	VtkAlgorithmPropertyCheckbox(const bool value, const QString& name,
	                             VtkAlgorithmProperties* algProps, QWidget* parent = 0);

	/// @brief Destructor.
	virtual ~VtkAlgorithmPropertyCheckbox();

private:
	const QString _name;
	VtkAlgorithmProperties* _algProps;

private slots:
	/// @brief This slots is automatically called when the checkbox state changed.
	void setNewValue(int state);
};

#endif // VTKALGORITHMPROPERTIESCHECKBOX_H
