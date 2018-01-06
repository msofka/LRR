IF(vxl_SOURCE_DIR)
  SET(VXL_MUL_INCLUDE_DIR ${MUL_INCLUDE_DIR})
ENDIF(vxl_SOURCE_DIR)

INCLUDE_DIRECTORIES( ${VXL_MUL_INCLUDE_DIR} )
INCLUDE_DIRECTORIES( ${VXL_RPL_INCLUDE_DIR} )
INCLUDE_DIRECTORIES( ${lrr_SOURCE_DIR} )

# Setup ITK
IF(BUILD_WITH_ITK)
  INCLUDE(${ITK_USE_FILE})
ENDIF(BUILD_WITH_ITK)

IF(BUILD_WITH_ITK)
 
  INCLUDE_DIRECTORIES( ${lrr_SOURCE_DIR}/cdcl/displayVTK )
  INCLUDE_DIRECTORIES( ${lrr_SOURCE_DIR}/cdcl )
  
  SET( rrlitk_sources								   
#    rrl_features_with_covariances.h              rrl_features_with_covariances.txx
#    rrl_shape_context3d.h                        rrl_shape_context3d.cxx								   
									   
    rrl_estimation_ICP_matching_all.h            rrl_estimation_ICP_matching_all.txx
    rrl_estimation_symmetric_ICP_matching_all.h  rrl_estimation_symmetric_ICP_matching_all.txx
    itkLocationRegistration.h                    itkLocationRegistration.cxx
  )
  
  #INCLUDE_DIRECTORIES( ${lrr_SOURCE_DIR}/r3p )
   
  AUX_SOURCE_DIRECTORY( Templates_itk rrlitk_sources )
  AUX_SOURCE_DIRECTORY( Templates rrlitk_sources )
  ADD_LIBRARY( rrlitk ${rrlitk_sources} )
  TARGET_LINK_LIBRARIES( rrlitk cdcl_io cdcl ITKCommon ITKIO vtkIO vtkFiltering vtkCommon libsvm )
ENDIF(BUILD_WITH_ITK)

IF( BUILD_EXAMPLES )
  SUBDIRS( examples )
ENDIF( BUILD_EXAMPLES )


