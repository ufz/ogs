/**
 * \file VtkAlgoritmProperties.h
 * 24/03/2010 KR Initial implementation
 *
 */


#ifndef VTKALGORITHMPROPERTIES_H
#define VTKALGORITHMPROPERTIES_H

// ** INCLUDES **
#include <QObject>
#include <vtkProperty.h>
#include <vtkTexture.h>
#include <QString>
#include <QMap>
#include <QVariant>
#include <QList>

#include "VtkColorLookupTable.h"

#define ogsUserPropertyMacro(name,type) \
virtual void Set##name (type _arg) \
{ \
	vtkDebugMacro(<< this->GetClassName() << " (" << this << "): setting " #name " to " << _arg); \
	if (this->name != _arg)  \
	{ \
		this->name = _arg; \
		this->Modified(); \
		(*(this->_algorithmUserProperties))[QString(#name)] = QVariant(_arg); \
	} \
} \
\
type name;
// End of ogsUserPropertyMacro

#define ogsUserVec2PropertyMacro(name,type) \
virtual void Set##name (type _arg1, type _arg2) \
{ \
	vtkDebugMacro(<< this->GetClassName() << " (" << this << "): setting " << #name " to (" << _arg1 << "," << _arg2 << ")"); \
	if ((this->name[0] != _arg1)||(this->name[1] != _arg2)) \
	{ \
		this->name[0] = _arg1; \
		this->name[1] = _arg2; \
		this->Modified(); \
		QList<QVariant> list; \
		list.push_back(QVariant(_arg1)); \
		list.push_back(QVariant(_arg2)); \
		(*(this->_algorithmUserVectorProperties))[QString(#name)] = list; \
	} \
} \
\
virtual void Set##name (type _arg[2]) \
{ \
	this->Set##name (_arg[0], _arg[1]);\
} \
\
type name[2];
// End of ogsUserVec2PropertyMacro

#define ogsUserVec3PropertyMacro(name,type) \
virtual void Set##name (type _arg1, type _arg2, type _arg3) \
{ \
	vtkDebugMacro(<< this->GetClassName() << " (" << this << "): setting " << #name " to (" << _arg1 << "," << _arg2 << "," << _arg3 << ")"); \
	if ((this->name[0] != _arg1)||(this->name[1] != _arg2)||(this->name[2] != _arg3)) \
	{ \
		this->name[0] = _arg1; \
		this->name[1] = _arg2; \
		this->name[2] = _arg3; \
		this->Modified(); \
		QList<QVariant> list; \
		list.push_back(QVariant(_arg1)); \
		list.push_back(QVariant(_arg2)); \
		list.push_back(QVariant(_arg3)); \
		(*(this->_algorithmUserVectorProperties))[QString(#name)] = list; \
	} \
} \
\
virtual void Set##name (type _arg[3]) \
{ \
	this->Set##name (_arg[0], _arg[1], _arg[2]);\
} \
\
type name[3];
// End of ogsUserVec3PropertyMacro

#define ogsUserVec4PropertyMacro(name,type) \
virtual void Set##name (type _arg1, type _arg2, type _arg3, type _arg4) \
{ \
	vtkDebugMacro(<< this->GetClassName() << " (" << this << "): setting " << #name " to (" << _arg1 << "," << _arg2 << "," << _arg3 << "," << _arg4 << ")"); \
	if ((this->name[0] != _arg1)||(this->name[1] != _arg2)||(this->name[2] != _arg3)||(this->name[3] != _arg4)) \
	{ \
		this->name[0] = _arg1; \
		this->name[1] = _arg2; \
		this->name[2] = _arg3; \
		this->name[3] = _arg4; \
		this->Modified(); \
		QList<QVariant> list; \
		list.push_back(QVariant(_arg1)); \
		list.push_back(QVariant(_arg2)); \
		list.push_back(QVariant(_arg3)); \
		list.push_back(QVariant(_arg4)); \
		(*(this->_algorithmUserVectorProperties))[QString(#name)] = list; \
	} \
} \
\
virtual void Set##name (type _arg[4]) \
{ \
	this->Set##name (_arg[0], _arg[1], _arg[2], _arg[3]);\
} \
\
type name[4];
// End of ogsUserVec4PropertyMacro

/**
 * \brief Contains properties for the visualization of objects as VtkVisPipelineItems
 */
class VtkAlgorithmProperties : public QObject
{
	Q_OBJECT

public:
	/// Constructor (sets default values)
	VtkAlgorithmProperties(QObject* parent = NULL) 
		: QObject(parent)
	{ 
		_property = vtkProperty::New(); 
		_texture  = NULL;
		_lut      = NULL;
		_scalarVisibility = true;
		_algorithmUserProperties = new QMap<QString, QVariant>;
		_algorithmUserVectorProperties = new QMap<QString, QList<QVariant> >;
	}

	virtual ~VtkAlgorithmProperties() 
	{
		_property->Delete();
		if (_texture != NULL) _texture->Delete();
		if (_lut != NULL) _lut->Delete();
		delete _algorithmUserProperties;
		delete _algorithmUserVectorProperties;
	};
	
	/// @brief Returns the vtk properties
	vtkProperty* GetProperties() const { return _property; };
	
	/// @brief Returns a texture (if one has been assigned).
	vtkTexture* GetTexture() { return _texture; };
	/// @brief Sets a texture for the VtkVisPipelineItem.
	void SetTexture(vtkTexture* t) { _texture = t; };

	/// @brief Returns the colour lookup table (if one has been assigned).
	VtkColorLookupTable* GetLookupTable() { return _lut; };
	/// @brief Sets a colour lookup table for the VtkVisPipelineItem.
	void SetLookUpTable(VtkColorLookupTable* lut) { _lut = lut; };
	
	/// Loads a predefined color lookup table from a file.
	void SetLookUpTable(const std::string &filename)
	{ 
		VtkColorLookupTable* colorLookupTable = VtkColorLookupTable::New();
		colorLookupTable->readFromFile(filename);
		colorLookupTable->setInterpolationType(VtkColorLookupTable::NONE);
		colorLookupTable->Build();
		SetLookUpTable(colorLookupTable);
	};

	/// @brief Returns the scalar visibility.
	bool GetScalarVisibility() const { return _scalarVisibility; }
	/// @brief Sets the scalar visibility.
	void SetScalarVisibility(bool on)
	{
		_scalarVisibility = on;
		emit ScalarVisibilityChanged(on);
	}
	
	/// @brief Returns the name. This is set to the file path if it is a source algorithm.
	QString GetName() const { return _name; }
	/// @brief Sets the name.
	void SetName(QString name) { _name = name; }
	
	/// @brief Returns a map of user properties.
	QMap<QString, QVariant>* GetAlgorithmUserProperties() const
	{
		return _algorithmUserProperties;
	}

	/// @brief Returns a map of vector user properties.
	QMap<QString, QList<QVariant> >* GetAlgorithmUserVectorProperties() const
	{
		return _algorithmUserVectorProperties;
	}

	/// @brief Sets a user property. This should be implemented by subclasses.
	virtual void SetUserProperty(QString name, QVariant value)
	{
		(*_algorithmUserProperties)[name] = value;
	}

	/// @brief Returns the value of a user property.
	QVariant GetUserProperty(QString name) const
	{
		if (this->_algorithmUserProperties->contains(name))
			return this->_algorithmUserProperties->value(name);
		else
		{
			std::cout << "Not a valid property: " << name.toStdString() << std::endl;
			return QVariant();
		}
	}

	/// @brief Sets a vector user property. This should be implemented by subclasses.
	virtual void SetUserVectorProperty(QString name, QList<QVariant> values)
	{
		(*_algorithmUserVectorProperties)[name] = values;
	}

	/// @brief Returns a list of values of a vector user property.
	QList<QVariant> GetUserVectorProperty(QString name) const
	{
		if (this->_algorithmUserVectorProperties->contains(name))
			return this->_algorithmUserVectorProperties->value(name);
		else
		{
			std::cout << "Not a valid property: " << name.toStdString() << std::endl;
			return QList<QVariant>();
		}
	}
	
protected:

	// Properties set on vtkActor
	vtkProperty* _property;
	vtkTexture* _texture;

	// Properties set on vtkMapper
	bool _scalarVisibility;
	VtkColorLookupTable* _lut;
	
	// Properties used in the GUI
	QString _name;
	
	QMap<QString, QVariant>* _algorithmUserProperties;
	QMap<QString, QList<QVariant> >* _algorithmUserVectorProperties;

signals:
	void ScalarVisibilityChanged(bool on);

};

#endif // VTKALGORITHMPROPERTIES_H
