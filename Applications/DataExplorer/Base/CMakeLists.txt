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
)

# Moc headers
set(MOC_HEADERS
	RecentFiles.h
	TreeModel.h
	ColorPickerPushButton.h
	# modeltest.h # Not needed
	CheckboxDelegate.h
	QValueTooltipSlider.h
)

# Header files
set(HEADERS
	LastSavedFileDirectory.h
	OGSError.h
	QNonScalableGraphicsTextItem.h
	TreeItem.h
	TreeModelIterator.h
	StrictDoubleValidator.h
)

# Put moc files in a project folder
source_group("UI Files" REGULAR_EXPRESSION "\\w*\\.ui")
source_group("Moc Files" REGULAR_EXPRESSION "moc_.*")

# Run Qts meta object compiler moc on header files
qt4_wrap_cpp(MOC_SOURCES ${MOC_HEADERS})

# Create the library
add_library(QtBase STATIC
	${SOURCES}
	${HEADERS}
	${MOC_HEADERS}
	${MOC_SOURCES}
)

target_link_libraries(QtBase PUBLIC Qt4::QtGui)

set_property(TARGET QtBase PROPERTY FOLDER "DataExplorer")
