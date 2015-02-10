set(CPACK_GENERATOR DragNDrop ZIP)
set(CPACK_DMG_FORMAT "UDBZ")

# See http://stackoverflow.com/a/16662169/80480 how to create the DS_Store file.
set(CPACK_DMG_BACKGROUND_IMAGE ${CMAKE_SOURCE_DIR}/Documentation/OpenGeoSys-Logo.png)
set(CPACK_DMG_DS_STORE ${CMAKE_SOURCE_DIR}/scripts/packaging/.DS_Store)
