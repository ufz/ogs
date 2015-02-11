# Get all files matching .ExternalData_<algo>_<hash>
file(GLOB_RECURSE FILES "" ${ExternalData_SOURCE_ROOT}/.ExternalData_${ExternalData_LINK_CONTENT}_*)
foreach(HASH_FILE ${FILES})
	string(REGEX MATCH [^_]+$ HASH ${HASH_FILE})
	message("Copying ${HASH_FILE} to ${ExternalData_OBJECT_STORES}/${HASH}")
	file(RENAME ${HASH_FILE} ${ExternalData_OBJECT_STORES}/${ExternalData_LINK_CONTENT}/${HASH})
endforeach()
