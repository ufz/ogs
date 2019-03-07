# Source files
set(SOURCES
    OGSError.cpp
    QNonScalableGraphicsTextItem.cpp
    RecentFiles.cpp
    TreeItem.cpp
    TreeModel.cpp
    ColorPickerPushButton.cpp
    TreeModelIterator.cpp
    CheckboxDelegate.cpp
    QValueTooltipSlider.cpp
)

# Header files
set(HEADERS
    LastSavedFileDirectory.h
    OGSError.h
    QNonScalableGraphicsTextItem.h
    TreeItem.h
    TreeModelIterator.h
    StrictDoubleValidator.h
    RecentFiles.h
    TreeModel.h
    ColorPickerPushButton.h
    CheckboxDelegate.h
    QValueTooltipSlider.h
)

# Put moc files in a project folder
source_group("UI Files" REGULAR_EXPRESSION "\\w*\\.ui")
source_group("Moc Files" REGULAR_EXPRESSION "moc_.*")

# Create the library
add_library(QtBase
    ${SOURCES}
    ${HEADERS}
)
if(BUILD_SHARED_LIBS)
    install(TARGETS QtBase LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}}
                           RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()

target_link_libraries(QtBase Qt5::Widgets)

set_property(TARGET QtBase PROPERTY FOLDER "DataExplorer")

if(OGS_USE_PCH)
    # Workaround for error on MSVC with BUILD_SHARED_LIBS=ON
    # 1>exports.def : warning LNK4022: cannot find unique match for symbol '__'
    # 1>exports.def : warning LNK4002: __@@_PchSym_@00@UfhvihUqvmprmhUxlwvUlthGUyfrowRivovzhvUzkkorxzgrlmhUurovrlUzkkorxzgrlmhurovrlOwriUivovzhvUzhxrrizhgvirmgviuzxvOlyq@4B2008FD98C1DD4 defined in ApplicationsFileIO.dir\Release\AsciiRasterInterface.obj
    # 1>exports.def : warning LNK4002: __@@_PchSym_@00@UfhvihUqvmprmhUxlwvUlthGUyfrowRivovzhvUzkkorxzgrlmhUwzgzvckoliviUyzhvUjgyzhvOwriUivovzhvUlthviiliOlyq@4B2008FD98C1DD4 defined in QtBase.dir\Release\OGSError.obj
    # 1>exports.def : error LNK2001: unresolved external symbol __
    if(NOT MSVC)
        cotire(QtBase)
    endif()
endif()
