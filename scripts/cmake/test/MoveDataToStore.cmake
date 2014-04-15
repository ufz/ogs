# Get all files matching .ExternalData_<algo>_<hash>
FILE(GLOB_RECURSE FILES "" ${ExternalData_SOURCE_ROOT}/.ExternalData_${ExternalData_LINK_CONTENT}_*)
FOREACH(HASH_FILE ${FILES})
	STRING(REGEX MATCH [^_]+$ HASH ${HASH_FILE})
	MESSAGE("Copying ${HASH_FILE} to ${ExternalData_OBJECT_STORES}/${HASH}")
	FILE(RENAME ${HASH_FILE} ${ExternalData_OBJECT_STORES}/${ExternalData_LINK_CONTENT}/${HASH})
ENDFOREACH()
