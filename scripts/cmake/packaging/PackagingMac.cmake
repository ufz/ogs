SET(CPACK_GENERATOR DragNDrop ZIP)
SET(CPACK_DMG_FORMAT "UDBZ")

# See http://stackoverflow.com/a/16662169/80480 how to create the DS_Store file.
SET(CPACK_DMG_BACKGROUND_IMAGE ${CMAKE_SOURCE_DIR}/Documentation/OpenGeoSys-Logo.png)
SET(CPACK_DMG_DS_STORE ${CMAKE_SOURCE_DIR}/scripts/packaging/.DS_Store)
