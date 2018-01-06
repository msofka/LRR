

FIND_PACKAGE(ITK REQUIRED)
INCLUDE(${ITK_USE_FILE})


FIND_PACKAGE(VTK REQUIRED)
INCLUDE (${VTK_USE_FILE})

INCLUDE_DIRECTORIES( ${VXL_RPL_INCLUDE_DIR}/rsdl )

SET( FeatureExtraction_SRCS
     itkFeatureImageFilter.h                  itkFeatureImageFilter.txx
     itkMeshSpatialFilterClean.h              itkMeshSpatialFilterClean.txx
     vtkITKPointSetToPolyDataFilter.h         vtkITKPointSetToPolyDataFilter.cxx
     ExtractFeatures.cxx
  )

ADD_EXECUTABLE(ExtractFeatures ${FeatureExtraction_SRCS} )
TARGET_LINK_LIBRARIES(ExtractFeatures
                      ITKNumerics ITKIO
                      vtkIO vtkFiltering vtkCommon
)

SET( KeypointExtraction_SRCS
     itkFeatureImageFilter.h                  itkFeatureImageFilter.txx
     itkMeshSpatialFilterClean.h              itkMeshSpatialFilterClean.txx
     vtkITKPointSetToPolyDataFilter.h         vtkITKPointSetToPolyDataFilter.cxx
     ExtractKeypoints.cxx
  )

ADD_EXECUTABLE(ExtractKeypoints ${KeypointExtraction_SRCS} )
TARGET_LINK_LIBRARIES(ExtractKeypoints
                      ITKNumerics ITKIO
                      vtkIO vtkFiltering vtkCommon
)


ADD_EXECUTABLE(WatershedSegmentation1 WatershedSegmentation1.cxx )
TARGET_LINK_LIBRARIES(WatershedSegmentation1 ITKAlgorithms ITKIO
)

ADD_EXECUTABLE(CropDicom CropDicom.cxx )
TARGET_LINK_LIBRARIES(CropDicom ITKNumerics ITKIO)

#SET( ComputeKeypoints_SRCS
#  #vtkFeaturesToPolyDataFilter.h            vtkFeaturesToPolyDataFilter.cxx
#  ComputeKeypoints.cxx
#  )

#ADD_EXECUTABLE(ComputeKeypoints ${ComputeKeypoints_SRCS} )
#TARGET_LINK_LIBRARIES(ComputeKeypoints
#                      ITKNumerics ITKIO
#                      vtkIO vtkFiltering vtkCommon
#)

SET( ComputeDescriptors_SRCS
  itkDescriptorMeshFilter.h                   itkDescriptorMeshFilter.txx
  vtkITKDescriptorPointSetToPolyDataFilter.h  vtkITKDescriptorPointSetToPolyDataFilter.cxx
  ComputeDescriptors.cxx
  )
  
ADD_EXECUTABLE(ComputeDescriptors ${ComputeDescriptors_SRCS} )
TARGET_LINK_LIBRARIES(ComputeDescriptors
                      ITKNumerics ITKIO
                      vtkIO vtkFiltering vtkCommon
)


ADD_EXECUTABLE(ExtractSlice ExtractSlice.cxx )
TARGET_LINK_LIBRARIES(ExtractSlice ITKIO
)




