/*
 * TemplateVec.h
 *
 *  Created on: Feb 26, 2010
 *      Author: TF
 */

#ifndef TEMPLATEVEC_H_
#define TEMPLATEVEC_H_

namespace GeoLib {

/**
 * \ingroup GeoLib
 *
 * \brief The class TemplateVec takes a unique name and manages
 * a std::vector of pointers to data elements of type T.
 *
 * Instances are classes PolylineVec and SurfaceVec.
 * */
template <class T> class TemplateVec
{
public:
	/**
	 * Constructor of class TemlateVec.
	 * @param name unique name of the project the elements belonging to.
	 * In order to access the data elements a unique name is required.
	 * @param data_vec vector of data elements
	 * @param elem_name_map Names of data elements can be given by a
	 * std::map<std::string, size_t>. Here the std::string is the name
	 * of the element and the value for size_t stands for an index in
	 * the data_vec.

	 */
	TemplateVec (const std::string &name, std::vector<T*> *data_vec, std::map<std::string, size_t>* elem_name_map = NULL) :
		_name(name), _data_vec(data_vec), _name_id_map (elem_name_map)
	{}

	/**
	 * destructor, deletes all data elements
	 */
	virtual ~TemplateVec ()
	{
		for (size_t k(0); k<size(); k++) delete (*_data_vec)[k];
		delete _data_vec;
		delete _name_id_map;
	}

	/** sets the name of the vector of geometric objects
	 * the data elements belonging to
	 * \param n the name as standard string */
	void setName (const std::string & n) { _name = n; }
	/**
	 * the name, the data element belonging to
	 * @return the name of the object
	 */
	std::string getName () const { return _name; }

	/**
	 * @return the number of data elements
	 */
	size_t size () const { return _data_vec->size(); }

	/**
	 * get a pointer to a standard vector containing the data elements
	 * @return the data elements
	 */
	const std::vector<T*>* getVector () const { return _data_vec; }

	/**
	 * search the vector of names for the ID of the geometric element with the given name
	 * @param name the name of the geometric element
	 * @param id the id of the geometric element
	 * @return
	 */
	bool getElementIDByName (const std::string& name, size_t &id) const
	{
		std::map<std::string,size_t>::const_iterator it (_name_id_map->find (name));

		if (it != _name_id_map->end()) {
			id = it->second;
			return true;
		} else return false;
	}

	const T* getElementByName (const std::string& name) const
	{
		size_t id;
		bool ret (getElementIDByName (name, id));
		if (ret) {
			return (*_data_vec)[id];
		} else {
			return NULL;
		}
	}

	/**
	 * The method returns true if there is a name associated
	 * with the given id, else method returns false.
	 * @param id the id
	 * @param element_name if a name associated with the id
	 * is found name is assigned to element_name
	 * @return if there is name associated with the given id:
	 * true, else false
	 */
	bool getNameOfElementByID (size_t id, std::string& element_name) const
	{
		if (! _name_id_map) return false;
		// search in map for id
		std::map<std::string,size_t>::const_iterator it (_name_id_map->begin());
		while (it != _name_id_map->end()) {
			if (it->second == id) {
				element_name = it->first;
				return true;
			}
			it++;
		}
		return false;
	}

	void setNameOfElementByID (size_t id, std::string& element_name)
	{
		if (! _name_id_map) return;
		_name_id_map->insert( std::pair<std::string, size_t>(element_name, id) );
	}

	/**
	 * The method returns true if the given element of type T
	 * can be found and the element has a name, else method returns false.
	 * @param data the data element, one wants to know the name
	 * @param name the name of the data element (if the data element is
	 * found and the data element has a name)
	 * @return if element is found and has a name: true, else false
	 */
	bool getNameOfElement (const T* data, std::string& name) const
	{
		for (size_t k(0); k<_data_vec->size(); k++) {
			if ((*_data_vec)[k] == data) {
				return getNameOfElementByID (k, name);
			}
		}
		return false;
	}

	void push_back (T* data_element, std::string const * const name = NULL)
	{
		_data_vec->push_back (data_element);
		if (name == NULL) return;
		if (! name->empty()) {
			if (_name_id_map == NULL) {
				_name_id_map = new std::map <std::string, size_t>;
			}
			_name_id_map->insert (std::pair<std::string,size_t>(*name, _data_vec->size()-1));
		}
	}


private:
	/** copy constructor doesn't have an implementation */
	// compiler does not create a (possible unwanted) copy constructor
	TemplateVec (const TemplateVec &);
	/** assignment operator doesn't have an implementation */
	// this way the compiler does not create a (possible unwanted) assignment operator
	TemplateVec& operator= (const TemplateVec& rhs);

	/** the name of the object */
	std::string _name;

	/**
	 * pointer to a vector of data elements
	 */
	std::vector <T*> *_data_vec;
	/**
	 * store names associated with the element ids
	 */
	std::map<std::string, size_t>* _name_id_map;
};

} // end namespace GeoLib

#endif /* TEMPLATEVEC_H_ */
