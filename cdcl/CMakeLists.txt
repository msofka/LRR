
INCLUDE_DIRECTORIES( ${VXL_RPL_INCLUDE_DIR} )
INCLUDE_DIRECTORIES( ${lrr_SOURCE_DIR}/cdcl/displayVTK )
INCLUDE_DIRECTORIES( ${lrr_SOURCE_DIR}/cdcl )
INCLUDE_DIRECTORIES( ${lrr_SOURCE_DIR} )

SET( estimation_sources
  cdcl_estimation_abs.h
  cdcl_estimation.h                            cdcl_estimation.txx
  cdcl_estimation_ICP.h                        cdcl_estimation_ICP.txx
  cdcl_estimation_ICP_matching_all.h           cdcl_estimation_ICP_matching_all.txx
  cdcl_estimation_symmetric_ICP_matching_all.h cdcl_estimation_symmetric_ICP_matching_all.txx  
  cdcl_estimation_symmetric.h                  cdcl_estimation_symmetric.txx
  cdcl_obj_fun_par.h                           cdcl_obj_fun_par.txx
  cdcl_obj_fun_par_sym.h                       cdcl_obj_fun_par_sym.txx
  cdcl_obj_fun_cov.h                           cdcl_obj_fun_cov.txx
  cdcl_feature.h
  cdcl_feature_with_shape.h
  cdcl_feature_ICP.h
  cdcl_keypoint.h                              
  cdcl_match.h                 
  cdcl_macros.h                
  cdcl_trans.h                                 cdcl_trans.txx
  cdcl_trans_similarity2d.h                    cdcl_trans_similarity2d.cxx
  cdcl_trans_affine.h                          cdcl_trans_affine.txx
  cdcl_trans_rigid3d.h                         cdcl_trans_rigid3d.cxx
  cdcl_utils.h                                 cdcl_utils.txx
  )

IF(BUILD_WITH_VTK)
  FIND_PACKAGE(VTK QUIET)
  IF( EXISTS ${VTK_USE_FILE} )
    INCLUDE(${VTK_USE_FILE})
    ADD_DEFINITIONS( -DVTK_FOUND )
    # VTK dependent code
    SET(estimation_sources ${estimation_sources}
      cdcl_utils_VTK.h                  cdcl_utils_VTK.txx
    )
    AUX_SOURCE_DIRECTORY(TemplatesVTK estimation_sources)
  ENDIF( EXISTS ${VTK_USE_FILE} ) 
ENDIF(BUILD_WITH_VTK)

AUX_SOURCE_DIRECTORY(Templates estimation_sources)
AUX_SOURCE_DIRECTORY(TemplatesVTK estimation_sources)

# dependencies: 
# rsdl for kd_tree
ADD_LIBRARY( cdcl ${estimation_sources} )
IF( EXISTS ${VTK_USE_FILE} )
  TARGET_LINK_LIBRARIES( cdcl rrel rsdl ${VNL_LIB} vul vtkIO vtkCommon )
ELSE( EXISTS ${VTK_USE_FILE} )
  TARGET_LINK_LIBRARIES( cdcl rrel rsdl ${VNL_LIB} vul )
ENDIF( EXISTS ${VTK_USE_FILE} )

SUBDIRS(display)

IF( BUILD_EXAMPLES )
  SUBDIRS( examples )
ENDIF( BUILD_EXAMPLES )

SUBDIRS(displayVTK)

SUBDIRS(io)

