/**
 * \file VtkAlgorithmPropertyLineEdit.h
 * 18/10/2010 LB Initial implementation
 */

#ifndef VTKALGORITHMPROPERTYLINEEDIT_H
#define VTKALGORITHMPROPERTYLINEEDIT_H

#include <QLineEdit>
#include <QVariant>

class VtkAlgorithmProperties;
class QString;

/// @brief This QLineEdit sets a user property on the given VtkAlgorithmProperties
/// object automatically.
class VtkAlgorithmPropertyLineEdit : public QLineEdit
{
	Q_OBJECT
	
public:
	/// @brief Constructor.
	/// @param value The initial value of the text.
	/// @param name The name of the user property to set.
	/// @param type The type of the property.
	/// @param algProps The VtkAlgorithmProperties object.
	/// @param parent The parent widget.
	VtkAlgorithmPropertyLineEdit(const QString& contents, const QString& name,
		QVariant::Type type, VtkAlgorithmProperties* algProps, QWidget* parent = 0);
	virtual ~VtkAlgorithmPropertyLineEdit();

private:
	const QString _name;
	VtkAlgorithmProperties* _algProps;
	QVariant::Type _type;
	
private slots:
	/// @brief This slots is automatically called when the text changed.
	void setNewValue();
};

#endif // VTKALGORITHMPROPERTYLINEEDIT_H
