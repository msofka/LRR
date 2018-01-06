
IF(BUILD_WITH_ITK)

#  ADD_EXECUTABLE( compute_shape_context3d
#    compute_shape_context3d.cxx
#  )
#  TARGET_LINK_LIBRARIES( compute_shape_context3d rrlitk cdcl vnl vnl_algo cdcl_io )
 
  ADD_EXECUTABLE( compute_distance_map
    compute_distance_map.cxx
  )
  TARGET_LINK_LIBRARIES( compute_distance_map cdcl ${VNL_LIB} cdcl_io ITKStatistics ITKCommon ITKIO vtkIO )

#IF( BUILD_THESE )
#  ADD_EXECUTABLE( compute_keypoints3d
#    compute_keypoints3d.cxx
#  )
#  TARGET_LINK_LIBRARIES( compute_keypoints3d rrlitk cdcl ${VNL_LIB} cdcl_io )
#ENDIF( BUILD_THESE )
 
  # include this because of ITK/vtk converting filters
  INCLUDE_DIRECTORIES( ${lrr_SOURCE_DIR}/cdcl/displayVTK )
  
  # for rsdl
  INCLUDE_DIRECTORIES( ${VXL_RPI_INCLUDE_DIR} )
 
  ADD_EXECUTABLE( location_registration
    location_registration.cxx
  )
  TARGET_LINK_LIBRARIES( location_registration ${VNL_LIB} rrlitk cdcl_io cdcl displayVTK libsvm )
  ADD_EXECUTABLE( indexing_one_descriptor
    indexing_one_descriptor.cxx
  )
  TARGET_LINK_LIBRARIES( indexing_one_descriptor rrlitk cdcl ${VNL_LIB} displayVTK ITKStatistics )
    
  ADD_EXECUTABLE( feature_repeatability
    feature_repeatability.cxx
  )
  TARGET_LINK_LIBRARIES( feature_repeatability cdcl vnl vnl_algo cdcl_io ) 

  ADD_EXECUTABLE( indexing_shape_context3dITK
    indexing_shape_context3dITK.cxx
  )
  TARGET_LINK_LIBRARIES( indexing_shape_context3dITK rrlitk cdcl ${VNL_LIB} )

  ADD_EXECUTABLE( compare_descriptors
    compare_descriptors.cxx
  )
  TARGET_LINK_LIBRARIES( compare_descriptors rrlitk cdcl ${VNL_LIB} )
  

  ADD_EXECUTABLE(NoduleRegistration2 NoduleRegistration2.cxx )
  TARGET_LINK_LIBRARIES(NoduleRegistration2 ITKNumerics ITKIO displayVTK ITKStatistics)  
  
  ADD_EXECUTABLE(MultiResImageRegistration2 MultiResImageRegistration2.cxx )
  TARGET_LINK_LIBRARIES(MultiResImageRegistration2 ITKNumerics ITKIO ITKStatistics)  
    
ENDIF(BUILD_WITH_ITK)



