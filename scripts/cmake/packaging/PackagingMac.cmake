option(OGS_PACK_DMG "package targets creates a .dmg disk image instead of .tar.gz" FALSE)
if(OGS_PACK_DMG)
    if(NOT OGS_BUILD_GUI)
        message(WARNING "OGS_PACK_DMG requires OGS_BUILD_GUI=ON!")
    endif()
    set(CPACK_GENERATOR DragNDrop)
endif()
set(CPACK_DMG_FORMAT "UDBZ")

# See http://stackoverflow.com/a/16662169/80480 how to create the DS_Store file.
set(CPACK_DMG_BACKGROUND_IMAGE ${CMAKE_SOURCE_DIR}/Documentation/OpenGeoSys-Logo.png)
set(CPACK_DMG_DS_STORE ${CMAKE_SOURCE_DIR}/scripts/packaging/.DS_Store)
