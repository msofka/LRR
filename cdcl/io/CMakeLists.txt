

IF(BUILD_WITH_ITK)
  IF(BUILD_WITH_VTK)
    FIND_PACKAGE(ITK)
    IF(ITK_FOUND)
      INCLUDE(${ITK_USE_FILE})
    ELSE(ITK_FOUND)
      MESSAGE(FATAL_ERROR
              "Cannot build CDC GUI without ITK.  Please set ITK_DIR.")
    ENDIF(ITK_FOUND)
       
    FIND_PACKAGE(VTK)
    IF (VTK_FOUND)
        INCLUDE (${VTK_USE_FILE})
    ELSE(VTK_FOUND)
      MESSAGE(FATAL_ERROR
              "Cannot build CDC GUI without VTK.  Please set VTK_DIR.")
    ENDIF (VTK_FOUND)
     
    SET( cdcl_io_sources
#      vtkPolyDataToFeaturesFilter.h            vtkPolyDataToFeaturesFilter.cxx
      vtkPolyDataToFeaturesWithShapeFilter.h   vtkPolyDataToFeaturesWithShapeFilter.cxx
      vtkPolyDataToFeaturesICPFilter.h         vtkPolyDataToFeaturesICPFilter.cxx
      itkImageSlicesWithGeometryFilter.h       itkImageSlicesWithGeometryFilter.txx
    )  
    ADD_LIBRARY( cdcl_io ${cdcl_io_sources} )  
    TARGET_LINK_LIBRARIES(cdcl_io
                          ITKNumerics ITKIO
                          vtkIO vtkFiltering vtkCommon cdcl
    )       
            
  ENDIF(BUILD_WITH_VTK)
ENDIF(BUILD_WITH_ITK)

