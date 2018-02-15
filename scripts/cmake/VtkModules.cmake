# CLI modules
set(VTK_MODULES
    vtkIOXML
    CACHE INTERNAL "Required VTK / ParaView modules"
)

if(OGS_USE_MPI)
    set(VTK_MODULES ${VTK_MODULES} vtkIOParallelXML vtkParallelMPI)
endif()

if(OGS_BUILD_GUI)
    set(VTK_MODULES ${VTK_MODULES}
        vtkRenderingCore
        vtknetcdfcpp
        vtkIOLegacy
        vtkIOImage
        vtkRenderingAnnotation
        vtkFiltersExtraction
        vtkFiltersGeometry
        vtkFiltersTexture
        vtkFiltersModeling
        vtkFiltersSources
        vtkImagingCore
        vtkInteractionWidgets
        vtkInteractionStyle
        vtkIOExport
        vtkRenderingFreeType
    )
endif()

if(OGS_INSITU)
    set(VTK_MODULES ${VTK_MODULES} vtkPVPythonCatalyst)
endif()
