

IF(BUILD_WITH_ITK)
  IF(BUILD_WITH_VTK)
  
    FIND_PACKAGE(ITK)
    IF(ITK_FOUND)
      INCLUDE(${ITK_USE_FILE})
    ELSE(ITK_FOUND)
      MESSAGE(FATAL_ERROR
              "Cannot build CDC Display VTK without ITK.  Please set ITK_DIR.")
    ENDIF(ITK_FOUND)
    
    
    FIND_PACKAGE(VTK)
    IF (VTK_FOUND)
      INCLUDE (${VTK_USE_FILE})
    ELSE(VTK_FOUND)
      MESSAGE(FATAL_ERROR
              "Cannot build CDC Display VTK without VTK.  Please set VTK_DIR.")
    ENDIF(VTK_FOUND)

    SET( displayVTK_sources
         cdcl_extract_data.h       cdcl_extract_data.txx
         vtkTensorGlyphScaled.h    vtkTensorGlyphScaled.cxx
    )
        
    AUX_SOURCE_DIRECTORY(Templates  displayVTK_sources)
         
    ADD_LIBRARY( displayVTK ${displayVTK_sources} )         
    
    TARGET_LINK_LIBRARIES(
      displayVTK

      ITKBasicFilters
      ITKCommon
      ITKIO
#      ITKFEM
      ITKNumerics
      ITKMetaIO
      ITKStatistics

      vtkRendering
      vtkGraphics
      vtkHybrid
      vtkImaging
      vtkIO
      vtkFiltering
      vtkCommon
      vtkVolumeRendering

      )
    
    
  ENDIF(BUILD_WITH_VTK)
ENDIF(BUILD_WITH_ITK)


