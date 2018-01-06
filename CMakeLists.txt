
PROJECT(lrr)

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

if( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )
endif( COMMAND cmake_policy )

SET (LIBRARY_OUTPUT_PATH ${lrr_BINARY_DIR}/lib CACHE PATH "Single output directory for building all libraries.")

FIND_PACKAGE(VXL)
IF(VXL_FOUND)
  # include the standard stuff, such as Dart test
  SET(VXL_PROVIDE_STANDARD_OPTIONS 1)
  INCLUDE(${VXL_CMAKE_DIR}/UseVXL.cmake)

  SET (LIBRARY_OUTPUT_PATH ${lrr_BINARY_DIR}/lib/ CACHE PATH "Single output directory for building all libraries.")

  OPTION( BUILD_WITH_ITK "Build with ITK (required in rkpl and eic)" ON)
  OPTION( BUILD_WITH_VTK "Build with VTK (required in cdcl)" ON)
  OPTION( BUILD_CDC "Build Covariance Driven Correspondences" ON)
  OPTION( BUILD_LRR "Build Location Registration and Recognition" ON)  

  # If needed, find an ITK built with system VXL.
  IF(BUILD_WITH_ITK)
    FIND_PACKAGE(ITK)
    IF(ITK_FOUND)
      #IF(NOT ITK_USE_SYSTEM_VXL)
      #  MESSAGE(FATAL_ERROR "ITK must be built with ITK_USE_SYSTEM_VXL set ON.")
      #ENDIF(NOT ITK_USE_SYSTEM_VXL)
      
      # only include this when needed
      # INCLUDE(${ITK_USE_FILE})
	  #message( fatal_error ${ITK_USE_FILE} )
	  #SET( VNL_LIB itkvnl itkvnl_algo )
	  SET( VNL_LIB vnl vnl_algo )
	  
	  #GET_TARGET_PROPERTY(LIB_NAME itkvnl LOCATION)
	  #GET_TARGET_PROPERTY(vnl_prefix itkvnl PREFIX)
	  #GET_TARGET_PROPERTY(vnl_suffix itkvnl SUFFIX)
	  #SET(NEW_LIB_NAME ${vnl_prefix}vnl${vnl_suffix})

	  #ADD_CUSTOM_COMMAND(
	  #  TARGET itkvnl
	  #  POST_BUILD
	  #  COMMAND ${CMAKE_COMMAND} -E copy ${LIB_NAME} ${NEW_LIB_NAME}
	  #)
	  #message( fatal_error "${CMAKE_COMMAND} -E copy ${LIB_NAME} ${NEW_LIB_NAME}" )

 
      #CONFIGURE_FILE( ${ITK_DIR}/bin/itkvnl.lib ${ITK_DIR}/vnl.lib COPY_ONLY IMMEDIATE)
      #CONFIGURE_FILE( ${ITK_DIR}/bin/itkvnl_algo.lib ${ITK_DIR}/vnl_algo.lib COPY_ONLY IMMEDIATE)

    ELSE(ITK_FOUND)
      MESSAGE(FATAL_ERROR "ITK is required.  Set ITK_DIR.")
	  SET( VNL_LIB vnl vnl_algo )
    ENDIF(ITK_FOUND)
  ENDIF(BUILD_WITH_ITK)
  
  IF(BUILD_WITH_VTK)
    # VTK is added to the top of the hierarchy to make sure the link
    # directories it needs are provided to all libraries that link to
    # rpisrc libraries that chain to VTK libraries.
    FIND_PACKAGE(VTK QUIET)
    IF( EXISTS ${VTK_USE_FILE} )
      INCLUDE(${VTK_USE_FILE})
      ADD_DEFINITIONS( -DVTK_FOUND )
    ENDIF( EXISTS ${VTK_USE_FILE} )
  ENDIF(BUILD_WITH_VTK)

  INCLUDE_DIRECTORIES( ${lrr_SOURCE_DIR} )
  
  INCLUDE(${ITK_DIR}/UseITK.cmake)
  INCLUDE_DIRECTORIES( ${ITK_INCLUDE_DIRS} )
  
  # Covariance Driven Correspondences
  IF(BUILD_CDC)
    SUBDIRS( cdcl )
  ENDIF(BUILD_CDC)

  # Location Registration and Recognition
  IF(BUILD_LRR)
    SUBDIRS( rrl )
    SUBDIRS( libsvm )
    SUBDIRS( features )
  ENDIF(BUILD_LRR)

ELSE(VXL_FOUND)
  IF( VXL_BINARY_PATH )
    SET( VXL_DIR ${VXL_BINARY_PATH} )
    MESSAGE( "Use preferred VXL_DIR in place of VXL_BINARY_PATH" )
  ENDIF( VXL_BINARY_PATH )
ENDIF(VXL_FOUND)

