# Source files
set(SOURCES
    OGSError.cpp
    QNonScalableGraphicsTextItem.cpp
    RecentFiles.cpp
    TreeItem.cpp
    TreeModel.cpp
    ColorPickerPushButton.cpp
    TreeModelIterator.cpp
    # modeltest.cpp # Not needed
    CheckboxDelegate.cpp
    QValueTooltipSlider.cpp
    Color.cpp
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
    # modeltest.h # Not needed
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

target_link_libraries(QtBase
    Qt4::QtGui
    logog
    GeoLib)

set_property(TARGET QtBase PROPERTY FOLDER "DataExplorer")
