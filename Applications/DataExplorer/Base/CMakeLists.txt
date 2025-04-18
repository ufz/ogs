# Source files
set(SOURCES
    Utils.cpp
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
    Utils.h
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
ogs_add_library(QtBase ${SOURCES} ${HEADERS})

target_link_libraries(QtBase Qt5::Widgets)

set_property(TARGET QtBase PROPERTY FOLDER "DataExplorer")
