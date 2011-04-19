# Source files
set( SOURCES
	OGSError.cpp
	QNonScalableGraphicsTextItem.cpp
	RecentFiles.cpp
	TreeItem.cpp
	TreeModel.cpp
	ColorPickerPushButton.cpp
	TreeModelIterator.cpp
	modeltest.cpp
	CheckboxDelegate.cpp
	QValueTooltipSlider.cpp
)

# Moc headers
set( MOC_HEADERS
	RecentFiles.h
	TreeModel.h
	ColorPickerPushButton.h
	modeltest.h
	CheckboxDelegate.h
	QValueTooltipSlider.h
)

# Header files
set( HEADERS
	OGSError.h
	QNonScalableGraphicsTextItem.h
	TreeItem.h
	TreeModelIterator.h
	StrictDoubleValidator.h
)

# UI files
set ( UIS
)

# Put moc files in a project folder
source_group("UI Files" REGULAR_EXPRESSION "\\w*\\.ui")
source_group("Moc Files" REGULAR_EXPRESSION "moc_.*")

# Run Qts user interface compiler uic on .ui files
qt4_wrap_ui( UI_HEADERS ${UIS} )

# Run Qts meta object compiler moc on header files
qt4_wrap_cpp( MOC_SOURCES ${MOC_HEADERS} )

# Create the library
add_library( QtBase STATIC
	${SOURCES}
	${HEADERS}
	${MOC_HEADERS}
	${MOC_SOURCES}
	${UIS}
)

target_link_libraries( QtBase ${QT_LIBRARIES} )
