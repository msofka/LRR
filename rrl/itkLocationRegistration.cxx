#ifndef __itkLocationRegistration_cxx
#define __itkLocationRegistration_cxx

#include <vul/vul_file.h>
#include <vul/vul_arg.h>

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <fstream>

#include "itkImageSeriesReader.h"
#include "itkImageSeriesWriter.h"
#include "itkNumericSeriesFileNames.h"
#include "itkRegularExpressionSeriesFileNames.h"

#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"

#include "itkRawImageIO.h"

#include "itkKdTree.h"
#include "itkKdTreeGenerator.h"
#include "itkListSample.h"

#include "itkImage.h"
#include "itkWarpImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkIdentityTransform.h"
#include "itkImageRegion.h"
#include "itkVectorLinearInterpolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"

#include "itkBSplineDeformableTransform.h"
#include "itkTransformFileReader.h"
#include "itkTransformFileWriter.h"
#include "itkTransformFactory.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkCenteredAffineTransform.h"
#include "itkResampleImageFilter.h"
#include "itkRegionOfInterestImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkCheckerBoardImageFilter.h"
#include "itkTileImageFilter.h"
#include "itkMaximumProjectionImageFilter.h"
#include "itkImageSeriesWriter.h"
#include "itkCastImageFilter.h"
#include "itkVTKImageIO.h"
#include "itkTimeProbe.h"

#include "itkDanielssonDistanceMapImageFilter.h"

#include "itkLightProcessObject.h"

#include "itkRescaleIntensityImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkMaximumProjectionImageFilter.h"
#include "itkResampleImageFilter.h"
#include "itkTileImageFilter.h"
#include "itkPNGImageIO.h"

#include <rrl/itkLocationRegistration.h>

#include <cdcl/cdcl_keypoint.h>
#include <cdcl/cdcl_estimation_ICP.h>
#include <cdcl/cdcl_estimation_ICP_matching_all.h>
#include <cdcl/cdcl_estimation_symmetric_ICP_matching_all.h>
#include <rrl/rrl_estimation_symmetric_ICP_matching_all.h>
#include <cdcl/cdcl_estimation_transfer.h>
#include <cdcl/cdcl_trans_affine.h>
#include <cdcl/cdcl_trans_rigid3d.h>
#include <cdcl/cdcl_trans.h>
#include <cdcl/cdcl_feature_with_shape.h>
#include <cdcl/cdcl_feature_ICP.h>
#include <cdcl/cdcl_match.h>
#include <cdcl/cdcl_utils.h>
#include <cdcl/cdcl_utils_VTK.h>

#include <cdcl/displayVTK/cdcl_extract_data.h>
#include <cdcl/io/vtkPolyDataToFeaturesWithShapeFilter.h>
#include <cdcl/io/vtkPolyDataToFeaturesICPFilter.h>

#include "vtkPolyData.h"
#include "vtkPolyDataWriter.h"
#include "vtkSmartPointer.h"

#include <cdcl/io/itkImageSlicesWithGeometryFilter.h>

#include <vnl/vnl_inverse.h>
#include <vnl/vnl_cross.h>
#include <vnl/vnl_trace.h>
#include <vnl/vnl_vector_fixed.txx>

VNL_VECTOR_FIXED_INSTANTIATE(double,12);

#include <libsvm/svm_interface.h>


//#include <rrl/rrl_cdcl_registration.h>

//:
// \file
// \brief  Register volume location using rrl_cdcl_registration object.
// \author Michal Sofka
// \date   Oct 2007



namespace {

template <class ImageType>
void WriteROIImage( typename ImageType::Pointer const &  Image,
                    itk::Point< double, 3 > const &  Point,
                    std::string  FileName )
{
  typedef typename itk::ImageRegion<3>  RegionOfInterestType;
  typedef itk::Point< double, 3 >  PointType;
  
  RegionOfInterestType  m_ROI;

  RegionOfInterestType::SizeType  size;
  RegionOfInterestType::IndexType  index;
  Image->TransformPhysicalPointToIndex( Point, index );

  std::cout << "Index of Point: " << index << std::endl;

  typename ImageType::SpacingType  spacing = Image->GetSpacing();

  size[0] = typename ImageType::SizeValueType( 100 / spacing[0] );
  size[1] = typename ImageType::SizeValueType( 100 / spacing[1] );
  size[2] = typename ImageType::SizeValueType( 100 / spacing[2] );

  index[0] = typename ImageType::IndexValueType( index[0] - size[0] / 2 );
  index[1] = typename ImageType::IndexValueType( index[1] - size[1] / 2 );
  index[2] = typename ImageType::IndexValueType( index[2] - size[2] / 2 );

  m_ROI.SetSize( size );
  m_ROI.SetIndex( index );

  bool canCrop = m_ROI.Crop( Image->GetLargestPossibleRegion() );
  if( !canCrop || !Image->GetLargestPossibleRegion().IsInside( index ) ) {
    std::cout << "Error: point outside the image." << std::endl;
    return;
  }


  // GenerateSlices
  //
  typedef itk::Image<signed short, 2>                           SliceImageType;
  typedef itk::ExtractImageFilter< ImageType, SliceImageType >  ExtractFilterType;
  typedef itk::ExtractImageFilter< ImageType, ImageType >       ExtractFilter3DType;
  typedef itk::MaximumProjectionImageFilter< ImageType, SliceImageType >  MaximumProjectionFilterType;

  // for some reason, having one filter and just changing the region didn't work here
  //
  // XY plane
  //
  SliceImageType::Pointer  movingSliceXY;
  {
  //ExtractFilter3DAdaptorType::Pointer  extractFilterXY = ExtractFilter3DAdaptorType::New();
  //ExtractFilter3DAdaptorType::InputImageRegionType  extractionRegion;

  //extractFilterXY->SetInput( m_MovingMappedImage );
  //extractFilterXY->SetInput( m_PointsImage );

  //// using m_ROI -- since m_MovingMappedImage is in the fixed coordinate system
  //ExtractFilter3DAdaptorType::InputImageRegionType::SizeType  size = m_ROI.GetSize();
  //ExtractFilter3DAdaptorType::InputImageRegionType::IndexType  index = m_ROI.GetIndex();
  //// add half size to get the slice in the middle of the ROI
  //index[2] += size[2] / 2 - 3;
  //// extracting Z slice
  //size[2] = 6;
  typename ExtractFilter3DType::Pointer  extractFilterXY = ExtractFilter3DType::New();
  typename ExtractFilter3DType::InputImageRegionType  extractionRegion;

  extractFilterXY->SetInput( Image );

  // using m_ROI -- since m_MovingMappedImage is in the fixed coordinate system
  typename ExtractFilter3DType::InputImageRegionType::SizeType  size = m_ROI.GetSize();
  typename ExtractFilter3DType::InputImageRegionType::IndexType  index = m_ROI.GetIndex();
  // add half size to get the slice in the middle of the ROI
  index[2] += size[2] / 2 - 3;
  // extracting Z slice
  size[2] = 6;

  extractionRegion.SetSize( size );
  extractionRegion.SetIndex( index );
  extractionRegion.Crop( Image->GetLargestPossibleRegion() );

  extractFilterXY->SetExtractionRegion( extractionRegion );

  typename MaximumProjectionFilterType::Pointer  maximumProjectionFilter = MaximumProjectionFilterType::New();
  maximumProjectionFilter->SetInput( extractFilterXY->GetOutput() );
  maximumProjectionFilter->Update();
  movingSliceXY = maximumProjectionFilter->GetOutput();
  movingSliceXY->DisconnectPipeline();
  }

  // factor by which the slice in XZ or YZ plane will be resampled along Z direction
  double resamplingFactor = spacing[2] / spacing[0];

  // YZ plane
  //
  SliceImageType::Pointer  movingSliceYZ;
  {
  typename ExtractFilterType::Pointer  extractFilterYZ = ExtractFilterType::New();
  typename ExtractFilterType::InputImageRegionType  extractionRegion;

  extractFilterYZ->SetInput( Image );

  // using m_ROI -- since m_MovingMappedImage is in the fixed coordinate system
  typename ExtractFilterType::InputImageRegionType::SizeType  size = m_ROI.GetSize();
  typename ExtractFilterType::InputImageRegionType::IndexType  index = m_ROI.GetIndex();
  // add half size to get the slice in the middle of the ROI
  index[0] += size[0] / 2;
  // extracting X slice
  size[0] = 0;

  extractionRegion.SetSize( size );
  extractionRegion.SetIndex( index );
  extractionRegion.Crop( Image->GetLargestPossibleRegion() );

  extractFilterYZ->SetExtractionRegion( extractionRegion );
  extractFilterYZ->Update();

  SliceImageType::Pointer  extractedYZ = extractFilterYZ->GetOutput();
  extractedYZ->DisconnectPipeline();

  // The resolution in the Z direction is small, so resample the image
  typedef itk::ResampleImageFilter< SliceImageType, SliceImageType >  ResamplerType;
  ResamplerType::Pointer movingMappedResampler = ResamplerType::New();
  movingMappedResampler->SetInput( extractedYZ );

  movingMappedResampler->SetOutputOrigin( extractedYZ->GetOrigin() );
  movingMappedResampler->SetOutputDirection( extractedYZ->GetDirection() );

  // Trick: Need to change the size of the index as well
  typename SliceImageType::IndexType    movingMappedIndex   = extractedYZ->GetLargestPossibleRegion().GetIndex();
  movingMappedIndex[1] = typename SliceImageType::IndexValueType( movingMappedIndex[1] * resamplingFactor );
  movingMappedResampler->SetOutputStartIndex( movingMappedIndex);

  typename SliceImageType::SizeType     movingMappedSize    = extractedYZ->GetLargestPossibleRegion().GetSize();
  movingMappedSize[1] = typename SliceImageType::SizeValueType( movingMappedSize[1] * resamplingFactor );
  movingMappedResampler->SetSize( movingMappedSize );

  SliceImageType::SpacingType  movingMappedSpacing = extractedYZ->GetSpacing();
  movingMappedSpacing[1] = movingMappedSpacing[1] / resamplingFactor;
  movingMappedResampler->SetOutputSpacing( movingMappedSpacing);

  movingMappedResampler->Update();

  movingSliceYZ = movingMappedResampler->GetOutput();
  movingSliceYZ->DisconnectPipeline();
  }


  // XZ plane
  //
  SliceImageType::Pointer  movingSliceXZ;
  {
  typename ExtractFilterType::Pointer  extractFilterXZ = ExtractFilterType::New();
  typename ExtractFilterType::InputImageRegionType  extractionRegion;

  extractFilterXZ->SetInput( Image );

  // using m_ROI -- since m_MovingMappedImage is in the fixed coordinate system
  typename ExtractFilterType::InputImageRegionType::SizeType  size = m_ROI.GetSize();
  typename ExtractFilterType::InputImageRegionType::IndexType  index = m_ROI.GetIndex();
  // add half size to get the slice in the middle of the ROI
  index[1] += size[1] / 2;
  // extracting Y slice
  size[1] = 0;

  extractionRegion.SetSize( size );
  extractionRegion.SetIndex( index );
  extractionRegion.Crop( Image->GetLargestPossibleRegion() );

  extractFilterXZ->SetExtractionRegion( extractionRegion );
  extractFilterXZ->Update();

  SliceImageType::Pointer  extractedXZ = extractFilterXZ->GetOutput();
  extractedXZ->DisconnectPipeline();

  // The resolution in the Z direction is small, so resample the image
  typedef itk::ResampleImageFilter< SliceImageType, SliceImageType >  ResamplerType;
  ResamplerType::Pointer movingMappedResampler = ResamplerType::New();
  movingMappedResampler->SetInput( extractedXZ );

  movingMappedResampler->SetOutputOrigin( extractedXZ->GetOrigin() );
  movingMappedResampler->SetOutputDirection( extractedXZ->GetDirection() );

  // Trick: Need to change the size of the index as well
  typename SliceImageType::IndexType    movingMappedIndex   = extractedXZ->GetLargestPossibleRegion().GetIndex();
  movingMappedIndex[1] = typename SliceImageType::IndexValueType( movingMappedIndex[1] * resamplingFactor );
  movingMappedResampler->SetOutputStartIndex( movingMappedIndex);

  typename SliceImageType::SizeType     movingMappedSize    = extractedXZ->GetLargestPossibleRegion().GetSize();
  movingMappedSize[1] = typename SliceImageType::SizeValueType( movingMappedSize[1] * resamplingFactor );
  movingMappedResampler->SetSize( movingMappedSize );

  typename SliceImageType::SpacingType  movingMappedSpacing = extractedXZ->GetSpacing();
  movingMappedSpacing[1] = movingMappedSpacing[1] / resamplingFactor;
  movingMappedResampler->SetOutputSpacing( movingMappedSpacing);

  movingMappedResampler->Update();

  movingSliceXZ = movingMappedResampler->GetOutput();
  movingSliceXZ->DisconnectPipeline();
  }


  // the trick is to tile 2d images into a 3d volume, so that it can be passed into vtk display in the base class
  typedef itk::TileImageFilter< SliceImageType, SliceImageType >  TileFilterType;
  //typedef itk::TileImageFilter< typename cdcl_extract_data_callback3d< dim, dof >::OutputImageType, typename cdcl_extract_data_callback3d< dim, dof >::OutputImageType >  TileFilterType;
  //typedef itk::TileImageFilter< typename cdcl_extract_data_callback3d< dim, dof >::OutputImageType, SliceImageType >  TileFilterType;
  TileFilterType::Pointer  tileFilter = TileFilterType::New();
  tileFilter->SetInput( 0, movingSliceXY );
  tileFilter->SetInput( 1, movingSliceYZ );
  tileFilter->SetInput( 2, movingSliceXZ );
  TileFilterType::LayoutArrayType  layout;
  layout[0] = 3;
  layout[1] = 1;
  tileFilter->SetLayout( layout );
  tileFilter->Update();


  SliceImageType::Pointer                                     m_SlicesImage;
  m_SlicesImage = tileFilter->GetOutput();



  // Write the checkerboard
  //
  typedef SliceImageType                  WriteImageType;
  typedef itk::Image<vxl_byte, 2>                                            ByteImageType;
  typedef itk::ImageFileWriter<ByteImageType>                               WriterType;

  typedef itk::RescaleIntensityImageFilter<WriteImageType, ByteImageType>  RescalerType;
  //typedef itk::RescaleIntensityImageFilter<SliceImageType, WriteImageType>  RescalerType;
  //  Rescaler to set the final image to 0..255
  RescalerType::Pointer rescaler = RescalerType::New();
  rescaler->SetOutputMinimum( 0 );
  rescaler->SetOutputMaximum( 255 );
  rescaler->SetInput( m_SlicesImage );

  // need to set IO because the filename might include periods from decimals
  // and the factory does not know which writer to use (based on extension)
  typedef itk::PNGImageIO  ImageIOType;
  ImageIOType::Pointer  io = ImageIOType::New();

  //  Writer
  WriterType::Pointer writer = WriterType::New();
  //writer->SetImageIO( io );
  writer->SetInput( rescaler->GetOutput() );

  writer->SetFileName( FileName.c_str() );
  writer->Update();

}



}



namespace {

template <class T, class Ran>
void rrel_util_median_and_scale( Ran first, Ran last,
                                 T& median, T& scale,
                                 int dof )
{
  long count = long(last-first); // VC6 & 7 has broken iterator_traits
  assert( count > 0 );
  assert( count > dof );

  //// ignore these values in median computation
  //T  notComputed = -1.0;
  //long notComputedCount = 0;

  //for ( Ran i=first; i!=last; ++i ) {
  //  if( *i == notComputed ) ++notComputedCount;
  //}

  Ran loc = first + count/2;// + notComputedCount;
  std::nth_element( first, loc, last );
  median = *loc;
  for ( Ran i=first; i!=last; ++i ) {
    *i = vnl_math_abs(*i-median);
  }
  ++loc;
  std::nth_element( first, loc, last );
  scale = T(1.4826 * (1 + 5.0/(count-dof)) * *loc);
}

}


namespace itk {

LocationRegistration::LocationRegistration()
{
  m_FixedImageDir = "";
  m_MovingImageDir = "";
  m_TransformFile = "";
  m_Prefix = "";
  m_SegmentationSuffix = "";

  locRegDir = "";
  matchesDirPath = "";

  m_FixedImage = 0;
  m_MovingImage = 0;

  kdTreeMovingDescriptors = 0;
  kdTreeMovingKeypoints = 0;

  m_FinalTransform = 0;

  estimation = 0;

  m_ROISize[0] = 50;
  m_ROISize[1] = 50;
  m_ROISize[2] = 50;
}



void LocationRegistration::GenerateSlices( 
                           vtkSmartPointer< vtkPolyData > const &  polyDataMovingMapped,
                           vtkSmartPointer< vtkPolyData > const &  polyDataFixed,
                           ImageType::Pointer const &  m_MovingMappedImage,
                           ImageType::Pointer const & m_FixedImage,
                           itk::ImageRegion< 3 > const &  m_ROI,
                           OutputImageType::Pointer &  m_SlicesImage )
{
  // GenerateSlices
  //
  typedef itk::ImageSlicesWithGeometryFilter< ImageType, OutputImageType >  ImageSlicesFilter;

  // Moving Images
  //
  // XY plane
  //
  OutputImageType::Pointer  movingSliceXY;

  ImageSlicesFilter::Pointer  movingSlicesFilterXY = ImageSlicesFilter::New();
  movingSlicesFilterXY->SetInput( m_MovingMappedImage );
  movingSlicesFilterXY->SetROI( m_ROI );
  movingSlicesFilterXY->SetPolyData( polyDataMovingMapped );
  movingSlicesFilterXY->SetSliceOrientation( ImageSlicesFilter::SLICE_ORIENTATION_XY );
  movingSlicesFilterXY->Update();
  movingSliceXY = movingSlicesFilterXY->GetOutput();
  movingSliceXY->DisconnectPipeline();

  // XZ plane
  //
  OutputImageType::Pointer  movingSliceXZ;

  ImageSlicesFilter::Pointer  movingSlicesFilterXZ = ImageSlicesFilter::New();
  movingSlicesFilterXZ->SetInput( m_MovingMappedImage );
  movingSlicesFilterXZ->SetROI( m_ROI );
  movingSlicesFilterXZ->SetPolyData( polyDataMovingMapped );
  movingSlicesFilterXZ->SetSliceOrientation( ImageSlicesFilter::SLICE_ORIENTATION_XZ );
  movingSlicesFilterXZ->Update();
  movingSliceXZ = movingSlicesFilterXZ->GetOutput();
  movingSliceXZ->DisconnectPipeline();

  // YZ plane
  //
  OutputImageType::Pointer  movingSliceYZ;

  ImageSlicesFilter::Pointer  movingSlicesFilterYZ = ImageSlicesFilter::New();
  movingSlicesFilterYZ->SetInput( m_MovingMappedImage );
  movingSlicesFilterYZ->SetROI( m_ROI );
  movingSlicesFilterYZ->SetPolyData( polyDataMovingMapped );
  movingSlicesFilterYZ->SetSliceOrientation( ImageSlicesFilter::SLICE_ORIENTATION_YZ );
  movingSlicesFilterYZ->Update();
  movingSliceYZ = movingSlicesFilterYZ->GetOutput();
  movingSliceYZ->DisconnectPipeline();


  // Fixed Images
  //
  // XY plane
  //
  OutputImageType::Pointer  fixedSliceXY;

  ImageSlicesFilter::Pointer  fixedSlicesFilterXY = ImageSlicesFilter::New();
  fixedSlicesFilterXY->SetInput( m_FixedImage );
  fixedSlicesFilterXY->SetROI( m_ROI );
  fixedSlicesFilterXY->SetPolyData( polyDataFixed );
  fixedSlicesFilterXY->SetSliceOrientation( ImageSlicesFilter::SLICE_ORIENTATION_XY );
  fixedSlicesFilterXY->Update();
  fixedSliceXY = fixedSlicesFilterXY->GetOutput();
  fixedSliceXY->DisconnectPipeline();

  // XZ plane
  //
  OutputImageType::Pointer  fixedSliceXZ;

  ImageSlicesFilter::Pointer  fixedSlicesFilterXZ = ImageSlicesFilter::New();
  fixedSlicesFilterXZ->SetInput( m_FixedImage );
  fixedSlicesFilterXZ->SetROI( m_ROI );
  fixedSlicesFilterXZ->SetPolyData( polyDataFixed );
  fixedSlicesFilterXZ->SetSliceOrientation( ImageSlicesFilter::SLICE_ORIENTATION_XZ );
  fixedSlicesFilterXZ->Update();
  fixedSliceXZ = fixedSlicesFilterXZ->GetOutput();
  fixedSliceXZ->DisconnectPipeline();

  // YZ plane
  //
  OutputImageType::Pointer  fixedSliceYZ;

  ImageSlicesFilter::Pointer  fixedSlicesFilterYZ = ImageSlicesFilter::New();
  fixedSlicesFilterYZ->SetInput( m_FixedImage );
  fixedSlicesFilterYZ->SetROI( m_ROI );
  fixedSlicesFilterYZ->SetPolyData( polyDataFixed );
  fixedSlicesFilterYZ->SetSliceOrientation( ImageSlicesFilter::SLICE_ORIENTATION_YZ );
  fixedSlicesFilterYZ->Update();
  fixedSliceYZ = fixedSlicesFilterYZ->GetOutput();
  fixedSliceYZ->DisconnectPipeline();


  typedef itk::CheckerBoardImageFilter<OutputImageType>  CheckerBoardType;
  CheckerBoardType::Pointer  CheckerboardFilterXY = CheckerBoardType::New();
  CheckerboardFilterXY->SetInput1( movingSliceXY );
  CheckerboardFilterXY->SetInput2( fixedSliceXY );
  CheckerboardFilterXY->Update();
  OutputImageType::Pointer  checkerXY = CheckerboardFilterXY->GetOutput();
  checkerXY->DisconnectPipeline();

  CheckerBoardType::Pointer  CheckerboardFilterXZ = CheckerBoardType::New();
  CheckerboardFilterXZ->SetInput1( movingSliceXZ );
  CheckerboardFilterXZ->SetInput2( fixedSliceXZ );
  CheckerboardFilterXZ->Update();
  OutputImageType::Pointer  checkerXZ = CheckerboardFilterXZ->GetOutput();
  checkerXZ->DisconnectPipeline();

  CheckerBoardType::Pointer  CheckerboardFilterYZ = CheckerBoardType::New();
  CheckerboardFilterYZ->SetInput1( movingSliceYZ );
  CheckerboardFilterYZ->SetInput2( fixedSliceYZ );
  CheckerboardFilterYZ->Update();
  OutputImageType::Pointer  checkerYZ = CheckerboardFilterYZ->GetOutput();
  checkerYZ->DisconnectPipeline();

  // the trick is to tile 2d images into a 3d volume, so that it can be passed into vtk display in the base class
  typedef itk::TileImageFilter< OutputImageType, OutputImageType >  TileFilterType;
  TileFilterType::Pointer  tileFilter = TileFilterType::New();
  tileFilter->SetInput( 0, movingSliceXY );
  tileFilter->SetInput( 1, fixedSliceXY );
  tileFilter->SetInput( 2, checkerXY );
  tileFilter->SetInput( 3, movingSliceXZ );
  tileFilter->SetInput( 4, fixedSliceXZ );
  tileFilter->SetInput( 5, checkerXZ );
  tileFilter->SetInput( 6, movingSliceYZ );
  tileFilter->SetInput( 7, fixedSliceYZ );
  tileFilter->SetInput( 8, checkerYZ );
  //tileFilter->SetInput( 0, movingSliceXY );
  //tileFilter->SetInput( 1, movingSliceXZ );
  //tileFilter->SetInput( 2, movingSliceYZ );
  //tileFilter->SetInput( 3, fixedSliceXY );
  //tileFilter->SetInput( 4, fixedSliceXZ );
  //tileFilter->SetInput( 5, fixedSliceYZ );
  //tileFilter->SetInput( 6, checkerXY );
  //tileFilter->SetInput( 7, checkerXZ );
  //tileFilter->SetInput( 8, checkerYZ );
  TileFilterType::LayoutArrayType  layout;
  layout[0] = 3;
  layout[1] = 3;
  tileFilter->SetLayout( layout );
  tileFilter->Update();


  m_SlicesImage = tileFilter->GetOutput();

}


void LocationRegistration::ReadImages()
{
  // Try using image reader first (for mhd file type, for example).
  //
  typedef itk::ImageFileReader< ImageType > FixedImageReaderType;
  FixedImageReaderType::Pointer fixedImageReader = FixedImageReaderType::New();
  fixedImageReader->SetFileName( m_FixedImageDir );

  fixedImageReader->Update();

  m_FixedImage = fixedImageReader->GetOutput();


  typedef itk::ImageFileReader< ImageType > MovingImageReaderType;
  MovingImageReaderType::Pointer movingImageReader = MovingImageReaderType::New();
  movingImageReader->SetFileName( m_MovingImageDir );

  movingImageReader->Update();

  m_MovingImage = movingImageReader->GetOutput();


  if( m_MovingImage.IsNull() || m_FixedImage.IsNull() ) {
    // Read images
    //
    typedef   itk::ImageSeriesReader< ImageType >  FixedReaderType;
    typedef   itk::ImageSeriesReader< ImageType >  MovingReaderType;


    // Read the input image.
    // USE GDCMImageIO, DICOMImageIO2 is OLD
    typedef itk::GDCMImageIO                        ImageIOType;
    typedef itk::GDCMSeriesFileNames                NamesGeneratorType;
    
    ImageIOType::Pointer io = ImageIOType::New();

    // Get the DICOM filenames from the directory
    NamesGeneratorType::Pointer  fixedNames = NamesGeneratorType::New();
    fixedNames->SetDirectory( m_FixedImageDir.c_str() );
    const FixedReaderType::FileNamesContainer &  fixedFileNames = 
                              fixedNames->GetInputFileNames();

    FixedReaderType::Pointer  fixedImageReader = FixedReaderType::New();
    fixedImageReader->SetFileNames( fixedFileNames );
    fixedImageReader->SetImageIO( io );
    std::cout << fixedNames;

    m_FixedImage  = fixedImageReader->GetOutput();
    fixedImageReader->Update();


    // Get the DICOM filenames from the directory
    NamesGeneratorType::Pointer  movingNames = NamesGeneratorType::New();
    movingNames->SetDirectory( m_MovingImageDir.c_str() );
    const MovingReaderType::FileNamesContainer &  movingFileNames = 
                              movingNames->GetInputFileNames();

    MovingReaderType::Pointer  movingImageReader = MovingReaderType::New();
    movingImageReader->SetFileNames( movingFileNames );
    movingImageReader->SetImageIO( io );
    std::cout << movingNames;

    m_MovingImage = movingImageReader->GetOutput();
    movingImageReader->Update();
  }

  std::cout << "Fixed image: " << m_FixedImage << std::endl;
  std::cout << "Moving image: " << m_MovingImage << std::endl;

  std::cout << "Fixed spacing: " << m_FixedImage->GetSpacing() << std::endl;
  std::cout << "Moving spacing: " << m_MovingImage->GetSpacing() << std::endl;
}


void LocationRegistration::ReadSegmentedImages()
{
  if( m_SegmentationSuffix == "" ) return;

#if 0
  // Read images
  //
  typedef   itk::ImageSeriesReader< SegmentedImageType >  FixedReaderType;
  typedef   itk::ImageSeriesReader< SegmentedImageType >  MovingReaderType;


  // Read the input image.
  // USE GDCMImageIO, DICOMImageIO2 is OLD
  //typedef itk::GDCMImageIO                        ImageIOType;
  //typedef itk::GDCMSeriesFileNames                NamesGeneratorType;
  typedef itk::RegularExpressionSeriesFileNames   NamesGeneratorType;
  const unsigned int subMatch = 0;
  std::string  regularExpression = ".*";
  
  //ImageIOType::Pointer io = ImageIOType::New();

  // Get the DICOM filenames from the directory
  NamesGeneratorType::Pointer  fixedNames = NamesGeneratorType::New();
  std::string dirNameFixed = m_FixedImageDir + m_SegmentationSuffix;
  std::cout << "Reading fixed segmented dir: " << dirNameFixed << std::endl;
  fixedNames->SetDirectory( dirNameFixed.c_str() );
  fixedNames->SetRegularExpression( regularExpression );
  fixedNames->SetSubMatch( subMatch );
  const FixedReaderType::FileNamesContainer &  fixedFileNames = 
                            fixedNames->GetFileNames();

  FixedReaderType::Pointer  fixedImageReader = FixedReaderType::New();
  fixedImageReader->SetFileNames( fixedFileNames );
//  fixedImageReader->SetImageIO( io );
  std::cout << fixedNames;

  m_FixedImageSegmented  = fixedImageReader->GetOutput();
  fixedImageReader->Update();


  // Get the DICOM filenames from the directory
  NamesGeneratorType::Pointer  movingNames = NamesGeneratorType::New();
  std::string dirNameMoving = m_MovingImageDir + m_SegmentationSuffix;
  std::cout << "Reading moving segmented dir: " << dirNameMoving << std::endl;
  movingNames->SetDirectory( dirNameMoving.c_str() );
  movingNames->SetRegularExpression( regularExpression );
  movingNames->SetSubMatch( subMatch );
  const MovingReaderType::FileNamesContainer &  movingFileNames = 
                            movingNames->GetFileNames();

  MovingReaderType::Pointer  movingImageReader = MovingReaderType::New();
  movingImageReader->SetFileNames( movingFileNames );
  //movingImageReader->SetImageIO( io );
  std::cout << movingNames;

  m_MovingImageSegmented = movingImageReader->GetOutput();
  movingImageReader->Update();
#endif

  typedef   itk::ImageFileReader< SegmentedImageType >  FixedReaderType;
  typedef   itk::ImageFileReader< SegmentedImageType >  MovingReaderType;

  //typedef  itk::RawImageIO< SegmentedImageType::PixelType, 3>   RawReaderType;
  //RawReaderType::Pointer  rawReader = RawReaderType::New();

  FixedReaderType::Pointer  fixedImageReader = FixedReaderType::New();
  std::string fileNameFixed = m_FixedImageDir + m_SegmentationSuffix + ".mhd";
  fixedImageReader->SetFileName( fileNameFixed );
  //fixedImageReader->SetImageIO( rawReader );
  m_FixedImageSegmented = fixedImageReader->GetOutput();
  fixedImageReader->Update();

  MovingReaderType::Pointer  movingImageReader = MovingReaderType::New();
  std::string fileNameMoving = m_MovingImageDir + m_SegmentationSuffix + ".mhd";
  movingImageReader->SetFileName( fileNameMoving );
  //movingImageReader->SetImageIO( rawReader );
  m_MovingImageSegmented = movingImageReader->GetOutput();
  movingImageReader->Update();


//save the segmentation as one image rather than series, because dicom can't be used (different image type)
//otherwise origin of Z is lost
//m_FixedImageSegmented->SetOrigin( m_FixedImage->GetOrigin() );
//m_MovingImageSegmented->SetOrigin( m_MovingImage->GetOrigin() );

  std::cout << "Fixed image segmented: " << m_FixedImageSegmented << std::endl;
  std::cout << "Moving image segmented: " << m_MovingImageSegmented << std::endl;

  std::cout << "Fixed Segmented spacing: " << m_FixedImageSegmented->GetSpacing() << std::endl;
  std::cout << "Moving Segmented spacing: " << m_MovingImageSegmented->GetSpacing() << std::endl;
}


int LocationRegistration::ReadVoronoiMaps()
{
  typedef   itk::Image< unsigned int, 3 >  UnsignedIntImageType;
  typedef   itk::ImageFileReader< UnsignedIntImageType >  MapFileReaderType;

  {
  MapFileReaderType::Pointer  mapFileReader = MapFileReaderType::New();

  std::string  fileName = m_FixedImageDir + "voronoi.mhd";
  mapFileReader->SetFileName( fileName );

  try
    {
    mapFileReader->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << "Error while reading the fixed voronoi map" << std::endl;
    std::cerr << excp << std::endl;
    std::cerr << "[FAILED]" << std::endl;
    return EXIT_FAILURE;
    }

  m_FixedVoronoiMap = mapFileReader->GetOutput();
  m_FixedVoronoiMap->DisconnectPipeline();
  }

  {
  MapFileReaderType::Pointer  mapFileReader = MapFileReaderType::New();

  std::string  fileName = m_MovingImageDir + "voronoi.mhd";
  mapFileReader->SetFileName( fileName );

  try
    {
    mapFileReader->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << "Error while reading the moving voronoi map" << std::endl;
    std::cerr << excp << std::endl;
    std::cerr << "[FAILED]" << std::endl;
    return EXIT_FAILURE;
    }

  m_MovingVoronoiMap = mapFileReader->GetOutput();
  m_MovingVoronoiMap->DisconnectPipeline();
  }

  return EXIT_SUCCESS;
}


void LocationRegistration::ReadFeatures()
{
  // Read features
  //

  {
  vtkSmartPointer< vtkXMLPolyDataReader >  featureReader = vtkSmartPointer< vtkXMLPolyDataReader >::New();
  vcl_string  filename = m_MovingImageDir + "_00.vtk";
  featureReader->SetFileName( filename.c_str() );
  featureReader->Update();  // need this update (bug in the converting filter)

  // create the converting filter
  vtkSmartPointer< PolyDataToFeaturesFilterType >  polyDataToFeaturesFilter = vtkSmartPointer< PolyDataToFeaturesFilterType >::New();
  polyDataToFeaturesFilter->SetInput( featureReader->GetOutput() );
  polyDataToFeaturesFilter->Update();
  vtkSmartPointer< FeatureAttributeSet >  featureAttributeSetVTK = polyDataToFeaturesFilter->GetOutput();

  typedef FeatureAttributeSet::FeatureSetType  FeatureSetType;
  FeatureSetType &  features = featureAttributeSetVTK->GetPoints();
  m_FeaturesMovingAll.push_back( features );
  }

  {
  vtkSmartPointer< vtkXMLPolyDataReader >  featureReader = vtkSmartPointer< vtkXMLPolyDataReader >::New();
  vcl_string  filename = m_FixedImageDir + "_00.vtk";
  featureReader->SetFileName( filename.c_str() );
  featureReader->Update();  // need this update (bug in the converting filter)

  // create the converting filter
  vtkSmartPointer< PolyDataToFeaturesFilterType >  polyDataToFeaturesFilter = vtkSmartPointer< PolyDataToFeaturesFilterType >::New();
  polyDataToFeaturesFilter->SetInput( featureReader->GetOutput() );
  polyDataToFeaturesFilter->Update();
  vtkSmartPointer< FeatureAttributeSet >  featureAttributeSetVTK = polyDataToFeaturesFilter->GetOutput();

  typedef FeatureAttributeSet::FeatureSetType  FeatureSetType;
  FeatureSetType &  features = featureAttributeSetVTK->GetPoints();
  m_FeaturesFixedAll.push_back( features );
  }
}


void LocationRegistration::ReadMovingKeypointDescriptors()
{
  // setup indexing
  //
  typedef itk::Array< float >                      DescriptorVectorType;
  typedef itk::Statistics::ListSample< DescriptorVectorType >       DescriptorSampleType;
  typedef itk::Statistics::KdTreeGenerator< DescriptorSampleType >  DescriptorTreeGeneratorType;
  typedef DescriptorTreeGeneratorType::KdTreeType                   DescriptorTreeType;


  DescriptorSampleType::Pointer  sampleMovingDescriptors = DescriptorSampleType::New();
  KeypointSampleType::Pointer    sampleMovingKeypoints = KeypointSampleType::New();

  vcl_vector< cdcl_keypoint< 3 >::sptr >                movingKeypoints;
  // read moving descriptors
  {
    vcl_vector< vnl_vector< float > >                    movingDescriptors;
    vcl_string  movingDescriptorFile = m_MovingImageDir + "desc.vtk";
    cdcl_read_keypoint_descriptors_VTK( movingKeypoints, movingDescriptors, movingDescriptorFile );


    // moving descriptor vectors
    //
    sampleMovingDescriptors->SetMeasurementVectorSize( movingDescriptors[0].size() );

    DescriptorVectorType  mv;
    mv.SetSize( movingDescriptors[0].size() );
    for( unsigned int n = 0 ; n < movingDescriptors.size(); ++n ) {
      for( unsigned int d = 0; d < movingDescriptors[n].size(); ++d ) {
        //mv.SetSize( movingDescriptors[n].size() );
        //mv.SetData( movingDescriptors[n].data_block() );
        mv.SetElement( d, movingDescriptors[n][d] );
      }
      sampleMovingDescriptors->PushBack( mv );
    }
    
    // recommended hack on releasing memory is to swap in empty vector
    vcl_vector< vnl_vector< float > >  empty_vector;
    empty_vector.swap( movingDescriptors );


    // moving keypoint locations
    //
    KeypointVectorType  mv_k;
    for( unsigned int n = 0 ; n < movingKeypoints.size(); ++n ) {
      mv_k[0] = movingKeypoints[n]->location_[0];
      mv_k[1] = movingKeypoints[n]->location_[1];
      mv_k[2] = movingKeypoints[n]->location_[2];
      //mv.Set_vnl_vector( m_FixedKeypoints[n]->location_.as_vector() );
      sampleMovingKeypoints->PushBack( mv_k );
    }

  }

  // moving keypoint descriptor vectors
  DescriptorTreeGeneratorType::Pointer  descriptorTreeGenerator = DescriptorTreeGeneratorType::New();

  descriptorTreeGenerator->SetSample( sampleMovingDescriptors );
  descriptorTreeGenerator->SetBucketSize( 16 );
  descriptorTreeGenerator->Update();

  kdTreeMovingDescriptors = descriptorTreeGenerator->GetOutput();

  vcl_cout << "Kd tree size: " << kdTreeMovingDescriptors->Size() << vcl_endl;


  // moving keypoint locations
  KeypointTreeGeneratorType::Pointer  keypointTreeGenerator = KeypointTreeGeneratorType::New();

  keypointTreeGenerator->SetSample( sampleMovingKeypoints );
  keypointTreeGenerator->SetBucketSize( 16 );
  keypointTreeGenerator->Update();

  kdTreeMovingKeypoints = keypointTreeGenerator->GetOutput();

}


void LocationRegistration::ReadFixedKeypointDescriptors()
{
  // read fixed descriptors
  //
  vcl_string  fixedDescriptorFile = m_FixedImageDir + "desc.vtk";
  cdcl_read_keypoint_descriptors_VTK( fixedKeypoints, fixedDescriptors, fixedDescriptorFile );
}


int LocationRegistration::SetupGroundTruthTransform()
{
  // Setup ground truth transforms  (either BSpline transform or deformation field is read)
  //
  deformationType = NONE;

  if( m_TransformFile != "" ) {
    vcl_string  transformExtension = itksys::SystemTools::GetFilenameLastExtension( m_TransformFile );
    // if warp field supplied
    if( transformExtension == ".vtk" || transformExtension == ".mhd" ) {
      typedef   itk::ImageFileReader< DeformationFieldType >  FieldReaderType;

      std::cout << "Reading deformation field: " << m_TransformFile;

      typedef itk::VTKImageIO                        ImageIOType;
      ImageIOType::Pointer io = ImageIOType::New();

      if( itksys::SystemTools::FileExists( m_TransformFile.c_str() ) && io->CanReadFile( m_TransformFile.c_str() ) ) {
        deformationType = FIELD;

        FieldReaderType::Pointer fieldReader = FieldReaderType::New();
        fieldReader->SetFileName( m_TransformFile );
        fieldReader->SetImageIO( io );
        fieldReader->Update();
        deformationField = fieldReader->GetOutput();

        std::cout << " of size " << deformationField->GetLargestPossibleRegion().GetSize() << std::endl;
      }
      else {
        std::cout << "Warning: The field " << m_TransformFile << " cannot be read." << std::endl;
      }

    }
    else { // assume BSpline transform supplied
      deformationType = BSPLINE;

      // In order to read a transform file, we instantiate a TransformFileReader. 
      // Like the writer, the reader is not templated.
      itk::TransformFileReader::Pointer reader;
      reader = itk::TransformFileReader::New();

      // Some transforms (like the BSpline transform) might not be registered 
      // with the factory so we add them manually. 
      itk::TransformFactory<BSplineTransformType>::RegisterTransform();

      // We then set the name of the file we want to read, and call the
      // Update() function.
      reader->SetFileName( m_TransformFile.c_str() );
    
      try
        {
        reader->Update();
        }
      catch( itk::ExceptionObject & excp )
        {
        std::cerr << "Error while reading the transform file" << std::endl;
        std::cerr << excp << std::endl;
        std::cerr << "[FAILED]" << std::endl;
        return EXIT_FAILURE;
        }

      // The transform reader is not template and therefore it retunrs a list
      // of \doxygen{Transform}. However, the reader instantiate the appropriate
      // transform class when reading the file but it is up to the user to
      // do the approriate cast.
      // To get the output list of transform we use the GetTransformList() function.
      typedef itk::TransformFileReader::TransformListType * TransformListType;
      TransformListType transforms = reader->GetTransformList();
      std::cout << "Number of transforms = " << transforms->size() << std::endl;

      // We then use an STL iterator to go trought the list of transforms. We show here
      // how to do the proper casting of the resulting transform.
      for( itk::TransformFileReader::TransformListType::const_iterator it = transforms->begin(); it != transforms->end(); ++it ) {
        if(!strcmp((*it)->GetNameOfClass(),"BSplineDeformableTransform"))
          {
          bsplineTransform = static_cast<BSplineTransformType*>((*it).GetPointer());
          bsplineTransform->Print(std::cout);
          }
        else
          {
            std::cout << "Do not know how to read: " << (*it)->GetNameOfClass() << std::endl;
            return -1;
          }
      }
    }
  }
  else { // set to identity
    std::cout << "No transform specified, errors will be printed without transforming keypoints." << std::endl;
    bsplineTransform = BSplineTransformType::New();

    BSplineTransformType::RegionType  region;
    BSplineTransformType::SizeType  size;
    size.Fill(10);
    region.SetSize(size);
    bsplineTransform->SetGridRegion( region );
    BSplineTransformType::OriginType  origin;
    origin.Fill ( 100 );
    bsplineTransform->SetGridOrigin ( origin );
    BSplineTransformType::SpacingType  spacing;
    spacing.Fill ( 1.5 );
    bsplineTransform->SetGridSpacing ( spacing );

    BSplineTransformType::ParametersType  parameters( bsplineTransform->GetNumberOfParameters() );
    bsplineTransform->SetParameters( parameters );

    bsplineTransform->SetIdentity();
  }

  // directory for the checker -- make sure it's in the patients folder and has the name of the series
  std::string  fixedVolume = itksys::SystemTools::GetFilenameName( m_FixedImageDir );
  std::string  movingVolume = itksys::SystemTools::GetFilenameName( m_MovingImageDir );
  std::string  fixedPatient = itksys::SystemTools::GetFilenameName( itksys::SystemTools::GetParentDirectory( m_FixedImageDir.c_str() ) );
  std::string  movingPatient = itksys::SystemTools::GetFilenameName( itksys::SystemTools::GetParentDirectory( m_MovingImageDir.c_str() ) );
  std::string  checkerDir;
  if( m_Prefix == "matches" )
    checkerDir = "locreg" + fixedPatient + "_" + fixedVolume + "-" + movingPatient + "_" + movingVolume;
  else
    checkerDir = itksys::SystemTools::GetFilenameWithoutLastExtension( m_Prefix ) + "reg" + fixedPatient + "_" + fixedVolume + "-" + movingPatient + "_" + movingVolume;
  std::string  patientDir = itksys::SystemTools::GetParentDirectory( itksys::SystemTools::GetParentDirectory( m_FixedImageDir.c_str() ).c_str() );

  if( patientDir != "" ) patientDir += "/";
  locRegDir = patientDir + checkerDir;
  std::cout << "Locreg dir: " <<  locRegDir << std::endl;

  itksys::SystemTools::MakeDirectory( locRegDir.c_str() );

  vcl_string  movingDesc = itksys::SystemTools::GetFilenameName( m_MovingImageDir ) + "desc";
  vcl_string  fixedDesc = itksys::SystemTools::GetFilenameName( m_FixedImageDir ) + "desc";
  std::string  matchesDir = itksys::SystemTools::GetFilenameWithoutLastExtension( m_Prefix ) + fixedPatient + "_" + fixedDesc + "-" + movingPatient + "_" + movingDesc;
  matchesDirPath = patientDir + matchesDir;

  if( !itksys::SystemTools::FileExists( matchesDirPath.c_str() ) ) {

    std::cout << "Matches dir does not exist: " << matchesDirPath << std::endl;
    return 1;
  }

  return 0;
}


cdcl_trans< LocationRegistration::dim, LocationRegistration::dof >::sptr
LocationRegistration
::GetCDCTransform()
{
  TransformType::ParametersType  currentParameters( m_FinalTransform->GetParameters() );

  unsigned NDimensions = 3;
  unsigned int par = 0;
  vnl_matrix_fixed< double, 3, 3 >  A;
  for(unsigned int row=0; row<NDimensions; row++) 
    {
    for(unsigned int col=0; col<NDimensions; col++) 
      {
      A[row][col] = currentParameters[par];
      ++par;
      }
    }

  // Transfer the rotation center 
  vnl_vector_fixed< double, 3 >  center;
  for(unsigned int i=0; i<NDimensions; i++) 
    {
    center[i] = currentParameters[par];
    ++par;
    }
  
  // Transfer the translation
  vnl_vector_fixed< double, 3 >  transl;
  for(unsigned int k=0; k<NDimensions; k++) 
    {
    transl[k] = currentParameters[par];
    ++par;
    }

  // CenteredAffineTransform  is:  P' = R x ( P - C ) + ( C + T ),
  // but in cdcl it's  P' = R x ( P - C ) + T1  so correct for that  ( T1 = C + T )
  vnl_vector_fixed< double, 3 >  translWCenter = transl + center;

  return new trans_type( A, translWCenter, center );
}


void LocationRegistration::UpdateITKTransform()
{
  TransformType::ParametersType  currentParameters( m_FinalTransform->GetNumberOfParameters() );

  unsigned NDimensions = 3;
  unsigned int par = 0;
  vnl_matrix_fixed< double, 3, 3 > const &  A = estimation->get_transform()->get_A();
  for(unsigned int row=0; row<NDimensions; row++) 
    {
    for(unsigned int col=0; col<NDimensions; col++) 
      {
      currentParameters[par] = A[row][col];
      ++par;
      }
    }

  // Transfer the rotation center 
  vnl_vector_fixed< double, 3 > const &  center = estimation->get_transform()->center_moving_;
  for(unsigned int i=0; i<NDimensions; i++) 
    {
    currentParameters[par] = center[i];
    ++par;
    }

  // CenteredAffineTransform is:  P' = R x ( P - C ) + ( C + T ), so need to add the center
  // but in cdcl it's  P' = R x ( P - C ) + T1  so correct for that  ( T1 = C + T )
  // Transfer the translation
  vnl_vector_fixed< double, 3 > const &  transl = estimation->get_transform()->get_translation();
  for(unsigned int m=0; m<NDimensions; m++) 
    {
    currentParameters[par] = transl[m] - center[m];
    ++par;
    }


  m_FinalTransform->SetParameters( currentParameters );
}


void LocationRegistration::SetupFinalTransform()
{
  // Setup Transform using descriptor match
  //
  // reset the transform using the current match (make sure it's properly centered)
  // CenteredAffineTransform is:  P' = R x ( P - C ) + ( C + T )
  // in our case:  P' = queryKeypoint
  //               C  = m_SelectedPointMoving
  //               P  = m_SelectedPointMoving
  //               T + C  = queryKeypoint
  TransformType::Pointer  m_KeypointTransform = TransformType::New();

  {
  m_KeypointTransform->SetIdentity();
  TransformType::InputPointType  center;
  center[0] = moving_keypoint->location_[0];
  center[1] = moving_keypoint->location_[1];
  center[2] = moving_keypoint->location_[2];
  m_KeypointTransform->SetCenter( center );

  TransformType::OutputVectorType  translation;
  translation[0] = -center[0] + fixed_keypoint->location_[0];
  translation[1] = -center[1] + fixed_keypoint->location_[1];
  translation[2] = -center[2] + fixed_keypoint->location_[2];
  m_KeypointTransform->SetTranslation( translation );

  vnl_matrix_fixed< double, 3, 3 >  movingR;
  movingR.set_column( 0, moving_keypoint->normal_ );
  movingR.set_column( 1, moving_keypoint->binormal_ );
  vnl_vector_fixed< double, 3 >  moving_normal_binormal = vnl_cross_3d( moving_keypoint->normal_, moving_keypoint->binormal_ );
  movingR.set_column( 2, moving_normal_binormal );

  vnl_matrix_fixed< double, 3, 3 >  fixedR;
  fixedR.set_column( 0, fixed_keypoint->normal_ );
  fixedR.set_column( 1, fixed_keypoint->binormal_ );
  vnl_vector_fixed< double, 3 >  fixed_normal_binormal = vnl_cross_3d( fixed_keypoint->normal_, fixed_keypoint->binormal_ );
  fixedR.set_column( 2, fixed_normal_binormal );

  // fixedR = R * movingR
  // R = fixedR * movingR^-1
  //vnl_matrix_fixed< double, 3, 3 >  R = movingR * vnl_inverse( fixedR );
  vnl_matrix_fixed< double, 3, 3 >  R = fixedR * vnl_inverse( movingR );
  itk::Matrix< double, 3, 3 >  Ritk = R.as_ref(); 
  m_KeypointTransform->SetMatrix( Ritk );
  }


  #if 0
  itk::Point< double, 3 >  movP;
  movP[0] = moving_keypoint->location_[0];
  movP[1] = moving_keypoint->location_[1];
  movP[2] = moving_keypoint->location_[2];

  itk::Point< double, 3 >  fixP;
  fixP[0] = fixed_keypoint->location_[0];
  fixP[1] = fixed_keypoint->location_[1];
  fixP[2] = fixed_keypoint->location_[2];

  WriteROIImage<ImageType>( m_MovingImage, movP, "slice_moving_key_moving_image.png" );
  WriteROIImage<ImageType>( m_FixedImage, fixP, "slice_fixed_key_fixed_image.png" );
  WriteROIImage<ImageType>( m_MovingImage, fixP, "slice_fixed_key_moving_image.png" );
  WriteROIImage<ImageType>( m_FixedImage, movP, "slice_moving_key_fixed_image.png" );
  #endif


  // Setup final transform (same as keypoint transform, centered on the query point). This and above can be combined.
  //
  m_FinalTransform = TransformType::New();

  TransformType::Pointer  m_KeypointInverseTransform = TransformType::New();
  TransformType::CenterType  KeypointTransformInverseCenter = m_KeypointTransform->TransformPoint( m_KeypointTransform->GetCenter() );
  m_KeypointInverseTransform->SetCenter( KeypointTransformInverseCenter );
  m_KeypointTransform->GetInverse( m_KeypointInverseTransform );

  queryPointMapped = m_KeypointInverseTransform->TransformPoint( queryPoint );


  m_FinalTransform->SetIdentity();
  m_FinalTransform->SetCenter( queryPointMapped );

  TransformType::OutputVectorType  translation = m_KeypointTransform->GetTranslation();
  m_FinalTransform->SetTranslation( translation );

  m_FinalTransform->SetMatrix( m_KeypointTransform->GetMatrix() );


  std::cout << "Query point: " << queryPoint << " Mapped: " << queryPointMapped << std::endl;

  //WriteROIImage<ImageType>( m_FixedImage, queryPointMapped, "slice_querymapped_fixed_image.png" );
  //WriteROIImage<ImageType>( m_MovingImage, queryPointMapped, "slice_querymapped_moving_image.png" );
}


void LocationRegistration::ReadFinalTransform( std::string TransformFile )
{
  itk::TransformFileReader::Pointer reader;
  reader = itk::TransformFileReader::New();
  itk::TransformFactory<TransformType>::RegisterTransform();

  reader->SetFileName( TransformFile );

  try
    {
    reader->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << "Error while reading the transform file" << std::endl;
    std::cerr << excp << std::endl;
    std::cerr << "[FAILED]" << std::endl;
    return;
    }

  typedef itk::TransformFileReader::TransformListType * TransformListType;
  TransformListType transforms = reader->GetTransformList();
  std::cout << "Number of transforms = " << transforms->size() << std::endl;

  for( itk::TransformFileReader::TransformListType::const_iterator it = transforms->begin(); it != transforms->end(); ++it ) {
    if(!strcmp((*it)->GetNameOfClass(),"CenteredAffineTransform"))
      {
      m_FinalTransform = static_cast<TransformType*>((*it).GetPointer());
      m_FinalTransform->Print(std::cout);
      }
  }
  if( !m_FinalTransform ) {
    std::cout << "Initializing transform was not assigned properly." << std::endl;
    return;
  }

  {
  TransformType::Pointer  m_FinalInverseTransform = TransformType::New();
  TransformType::CenterType  FinalTransformInverseCenter = m_FinalTransform->TransformPoint( m_FinalTransform->GetCenter() );
  m_FinalInverseTransform->SetCenter( FinalTransformInverseCenter );
  m_FinalTransform->GetInverse( m_FinalInverseTransform );

  queryPointMapped = m_FinalInverseTransform->TransformPoint( queryPoint );

  std::cout << "Query point: " << queryPoint << " mapped: " << queryPointMapped << std::endl;
  }
}


namespace {
// Function object for sorting vector of features and associated region ID's based on distance from a given point
class DistanceAtIndexIsGreater {
public:
  typedef std::pair< LocationRegistration::feature_sptr_type, LocationRegistration::SegmentedImageType::PixelType >  FeatureIDPair;

  DistanceAtIndexIsGreater( itk::Point< double, 3 > const &  ReferencePoint )
    : m_ReferencePoint( ReferencePoint ) {}
  bool operator()( FeatureIDPair const &  i, FeatureIDPair const &  j )
    {
      vnl_vector< double >  distancei = i.first->location_ - m_ReferencePoint.GetVnlVector();
      vnl_vector< double >  distancej = j.first->location_ - m_ReferencePoint.GetVnlVector();
      return distancei.magnitude() < distancej.magnitude();
    }
private:
  itk::Point< double, 3 > const &  m_ReferencePoint;
};


bool AnglesSmaller( const std::pair< double, unsigned int > &  left,
                    const std::pair< double, unsigned int > &  right )
{
  return left.first < right.first;
}


}


void LocationRegistration::FeaturesInROIs( vcl_vector< feature_sptr_type > &  moving_inside,
                                           vcl_vector< feature_sptr_type > &  fixed_inside,
                                           vnl_vector_fixed< double, dim >  const &  moving_x0,
                                           vnl_vector_fixed< double, dim >  const &  moving_x1,
                                           vnl_vector_fixed< double, dim >  const &  fixed_x0,
                                           vnl_vector_fixed< double, dim >  const &  fixed_x1 )
{
  moving_inside.clear();
  fixed_inside.clear();

  typedef vcl_vector< feature_sptr_type >  feature_vector_type;

  // only consider points that are in the current region
  //
  // extract fixed features from the current region only
  // create new copies so that normalization of each dataset individually is possible
  for( feature_vector_type::size_type  i = 0; i < m_FeaturesFixedAll[0].size(); ++i ) {
    vnl_vector_fixed< double, dim >  curr = m_FeaturesFixedAll[0][i]->location_;
    bool inside = true;
    for( unsigned int d = 0; d < dim; ++d ) if( fixed_x0[d] > curr[d] || curr[d] > fixed_x1[d] ) inside = false;

    if( inside ) {
      fixed_inside.push_back( new feature_type( *(m_FeaturesFixedAll[0][i]) ) );
		}
	}


  for( feature_vector_type::size_type  i = 0; i < m_FeaturesMovingAll[0].size(); ++i ) {
    vnl_vector_fixed< double, dim >  curr = m_FeaturesMovingAll[0][i]->location_;
    bool inside = true;
    for( unsigned int d = 0; d < dim; ++d ) if( moving_x0[d] > curr[d] || curr[d] > moving_x1[d] ) inside = false;

    if( inside ) {
      moving_inside.push_back( new feature_type( *(m_FeaturesMovingAll[0][i]) ) );
    }
  }
}


void LocationRegistration::FeaturesInRegions( vcl_vector< feature_sptr_type > &  moving_inside,
                                              vcl_vector< feature_sptr_type > &  fixed_inside,
                                              vnl_vector_fixed< double, dim >  const &  moving_x0,
                                              vnl_vector_fixed< double, dim >  const &  moving_x1,
                                              vnl_vector_fixed< double, dim >  const &  fixed_x0,
                                              vnl_vector_fixed< double, dim >  const &  fixed_x1 )
{
  moving_inside.clear();
  fixed_inside.clear();

  //WriteROIImage<SegmentedImageType>( m_MovingImageSegmented, queryPointMapped, "MovingSegmentedQueryPointMapped.png" );
  //WriteROIImage<SegmentedImageType>( m_FixedImageSegmented, queryPoint, "FixedSegmentedQueryPoint.png" );

  // features and their region ID
  typedef DistanceAtIndexIsGreater::FeatureIDPair  FeatureIDPair;

  typedef vcl_vector< FeatureIDPair >  feature_ID_vector_type;

  typedef vcl_vector< feature_sptr_type >  feature_vector_type;
  
  // get features in current region of interest
  feature_ID_vector_type  moving_in_roi;
  feature_ID_vector_type  fixed_in_roi;


  typedef itk::NearestNeighborInterpolateImageFunction< SegmentedImageType, double >  SegmentedImageInterpolatorType;
  SegmentedImageInterpolatorType::Pointer  interpolator = SegmentedImageInterpolatorType::New();
  interpolator->SetInputImage( m_FixedImageSegmented );

  SegmentedImageInterpolatorType::InputType  pointLocation;

  typedef SegmentedImageType::PixelType  SegmentedPixelType;

  // only consider points that are in the current region
  //
  // extract fixed features from the current region only
  // create new copies so that normalization of each dataset individually is possible
  for( feature_vector_type::size_type  i = 0; i < m_FeaturesFixedAll[0].size(); ++i ) {
    vnl_vector_fixed< double, dim >  curr = m_FeaturesFixedAll[0][i]->location_;
    bool inside = true;
    for( unsigned int d = 0; d < dim; ++d ) if( fixed_x0[d] > curr[d] || curr[d] > fixed_x1[d] ) inside = false;

    if( inside ) {
      // get current region ID
      pointLocation[0] = curr[0];
      pointLocation[1] = curr[1];
      pointLocation[2] = curr[2];
      SegmentedImageInterpolatorType::OutputType  regionID = interpolator->Evaluate( pointLocation );
      SegmentedImageInterpolatorType::IndexType   indexOfPoint;
      m_FixedImageSegmented->TransformPhysicalPointToIndex( pointLocation, indexOfPoint );
      //std::cout << pointLocation << "  " << indexOfPoint << "  " << regionID << "  " << m_FixedImageSegmented->GetPixel( indexOfPoint ) << std::endl;
//FIX: this does not need to create a new feature if moving_in_roi is not used outside this function -- but see comment at the beginning of this block
      fixed_in_roi.push_back( FeatureIDPair( new feature_type( *(m_FeaturesFixedAll[0][i]) ), SegmentedPixelType( regionID ) ) );
    }
  }
  
  typedef SegmentedImageInterpolatorType::OutputType  LabelType;
  LabelType  queryRegionID = interpolator->Evaluate( queryPoint );


  interpolator->SetInputImage( m_MovingImageSegmented );

  // extract moving features from the current region only
  // create new copies so that normalization of each dataset individually is possible
  itk::Point< double, 3 >  fixedRegionCentroid;
  fixedRegionCentroid.Fill( 0.0 );
  double numForCentroid = 0.0;
  for( feature_vector_type::size_type  i = 0; i < m_FeaturesMovingAll[0].size(); ++i ) {
    vnl_vector_fixed< double, dim >  curr = m_FeaturesMovingAll[0][i]->location_;
    bool inside = true;
    for( unsigned int d = 0; d < dim; ++d ) if( moving_x0[d] > curr[d] || curr[d] > moving_x1[d] ) inside = false;

    if( inside ) {
      // get current region ID
      pointLocation[0] = curr[0];
      pointLocation[1] = curr[1];
      pointLocation[2] = curr[2];
      SegmentedImageInterpolatorType::OutputType  regionID = interpolator->Evaluate( pointLocation );
//FIX: this does not need to create a new feature if moving_in_roi is not used outside this function -- but see comment at the beginning of this block
      feature_sptr_type  toAdd = new feature_type( *(m_FeaturesMovingAll[0][i]) );
      moving_in_roi.push_back( FeatureIDPair( toAdd, SegmentedPixelType( regionID ) ) );
      fixedRegionCentroid[0] += toAdd->location_[0];
      fixedRegionCentroid[1] += toAdd->location_[1];
      fixedRegionCentroid[2] += toAdd->location_[2];
      ++numForCentroid;
    }
  }
  fixedRegionCentroid[0] /= numForCentroid;
  fixedRegionCentroid[1] /= numForCentroid;
  fixedRegionCentroid[2] /= numForCentroid;

  //LabelType  queryMappedRegionID = interpolator->Evaluate( queryPointMapped );
  // the region in the other volume is found by mapping centroid of the region in the current volume
  LabelType  queryMappedRegionID = interpolator->Evaluate( fixedRegionCentroid );


  // sort features based on distance from the region center
  DistanceAtIndexIsGreater  distanceCompareFixed( queryPoint );
  std::sort( fixed_in_roi.begin(), fixed_in_roi.end(), distanceCompareFixed );

  DistanceAtIndexIsGreater  distanceCompareMoving( queryPointMapped );
  std::sort( moving_in_roi.begin(), moving_in_roi.end(), distanceCompareMoving );

  //std::cout << "Segmented moving: " << m_MovingImageSegmented << std::endl;
  //std::cout << "Segmented fixed: " << m_FixedImageSegmented << std::endl;

  // from features in the current region of interest, only select those in the oversegmented regions
  // first, add the region where the query point is
//  interpolator->SetInputImage( m_FixedImageSegmented );
  std::cout << "Query point region ID: " << queryRegionID << std::endl;

  // set of (unique) labels sorted based on distance from the center of the region
  // since set is sorted unique associative container and we need to preserve the ordering,
  // use two containers (the vector will have the preserved ordering)
  std::set< LabelType >  labels_set;
  std::vector< LabelType >  labels;
  for( feature_ID_vector_type::iterator  it = fixed_in_roi.begin(); it != fixed_in_roi.end(); ++it ) {
    std::pair< std::set< LabelType >::iterator, bool > result = labels_set.insert( it->second );
    // if this label is unique, add it
    if( result.second ) labels.push_back( it->second );
  }

  std::vector< LabelType >::iterator  label_it = labels.begin();
  std::cout << "First Label ID: " << *label_it << std::endl;
  while( fixed_inside.size() < 1500 && label_it != labels.end() ) {  // this threshold must be the same as the one outside this function
    queryRegionID = *label_it;
    for( feature_ID_vector_type::iterator  it = fixed_in_roi.begin(); it != fixed_in_roi.end(); ++it ) {
      //std::cout << "curr region ID: " << it->second << std::endl;
      if( queryRegionID == it->second ) fixed_inside.push_back( new feature_type( *(it->first) ) );
//std::cout << it->first->location_ << std::endl;
    }
    ++label_it;
  }
  std::cout << "Total in region: " << fixed_in_roi.size() << " Added: " << fixed_inside.size() << " fixed points." << std::endl;

//  interpolator->SetInputImage( m_MovingImageSegmented );
  std::cout << "Query point mapped region ID: " << queryMappedRegionID << std::endl;

  // set of (unique) labels sorted based on distance from the center of the region
  labels_set.clear();
  labels.clear();
  for( feature_ID_vector_type::iterator  it = moving_in_roi.begin(); it != moving_in_roi.end(); ++it ) {
    std::pair< std::set< LabelType >::iterator, bool > result = labels_set.insert( it->second );
    // if this label is unique, add it
    if( result.second ) labels.push_back( it->second );
  }

  label_it = labels.begin();
  std::cout << "First Label ID: " << *label_it << std::endl;
  while( moving_inside.size() < 1500 && label_it != labels.end() ) { // this threshold must be the same as the one outside this function
    queryMappedRegionID = *label_it;
    for( feature_ID_vector_type::iterator  it = moving_in_roi.begin(); it != moving_in_roi.end(); ++it ) {
      //std::cout << "curr region ID: " << it->second << std::endl;
      if( queryMappedRegionID == it->second ) moving_inside.push_back( new feature_type( *(it->first) ) );
    }
    ++label_it;
  }
  std::cout << "Total in region: " << moving_in_roi.size() << " Added: " << moving_inside.size() << " moving points." << std::endl;

}


bool LocationRegistration::SetupEstimation()
{
  // Setup moving ROI (region of interest)
  //
  vnl_vector_fixed< double, dim >  moving_x0;
  vnl_vector_fixed< double, dim >  moving_x1;
  {
  // convert from physical (size) to image coordinates
  //
  itk::ImageRegion<3>::SizeType  DesiredSize;

  // compute origin of the ROI given size and selected point
  itk::Point< double, 3 >  ROIOrigin;
  ImageType::SpacingType  ImageSpacing = m_MovingImage->GetSpacing();
  for( unsigned int d = 0; d < 3; ++d ) {
    ROIOrigin[d] = queryPointMapped[d] - m_ROISize[d] / 2;
    DesiredSize[d] = ImageType::SizeValueType( m_ROISize[d] / ImageSpacing[d] );
  }
  // convert ROI origin (in physical coordinates) to the corresponding index in the image
  itk::ImageRegion<3>::IndexType  DesiredIndex;
  m_MovingImage->TransformPhysicalPointToIndex( ROIOrigin, DesiredIndex );

  m_ROImoving.SetSize( DesiredSize );
  m_ROImoving.SetIndex( DesiredIndex );

  moving_x0 = ROIOrigin.GetVnlVector();
  moving_x1[0] = moving_x0[0] + m_ROISize[0];
  moving_x1[1] = moving_x0[1] + m_ROISize[1];
  moving_x1[2] = moving_x0[2] + m_ROISize[2];

  //vcl_cout << "Index: " << m_ROImoving.GetIndex() << " Size: " << m_ROImoving.GetSize() << vcl_endl;

  // adjust ROI so that it lies within the moving image
  if( m_ROImoving.Crop( m_MovingImage->GetLargestPossibleRegion() ) ) {
    //vcl_cout << "Making sure inverse mapped ROI lies within the moving image" << vcl_endl;
    itk::ImageRegion<3>::SizeType  DesiredSize = m_ROImoving.GetSize();
    itk::ImageRegion<3>::IndexType DesiredIndex = m_ROImoving.GetIndex();
    //vcl_cout << "ROI: " << DesiredIndex[0] << ".." << DesiredIndex[0]+DesiredSize[0] << ", "
    //                    << DesiredIndex[1] << ".." << DesiredIndex[1]+DesiredSize[1] << ", "
    //                    << DesiredIndex[2] << ".." << DesiredIndex[2]+DesiredSize[2] << vcl_endl;

  }
  }


  // Setup Fixed ROI
  //
  vnl_vector_fixed< double, dim >  fixed_x0;
  vnl_vector_fixed< double, dim >  fixed_x1;
  {
  // convert from physical (size) to image coordinates
  //
  itk::ImageRegion<3>::SizeType  DesiredSize;

  // compute origin of the ROI given size and selected point
  itk::Point< double, 3 >  ROIOrigin;
  ImageType::SpacingType  ImageSpacing = m_FixedImage->GetSpacing();
  for( unsigned int d = 0; d < 3; ++d ) {
    ROIOrigin[d] = queryPoint[d] - m_ROISize[d] / 2;
    DesiredSize[d] = ImageType::SizeValueType( m_ROISize[d] / ImageSpacing[d] );
  }
//ROIOrigin[0] =- 20;
  // convert ROI origin (in physical coordinates) to the corresponding index in the image
  itk::ImageRegion<3>::IndexType  DesiredIndex;
  m_FixedImage->TransformPhysicalPointToIndex( ROIOrigin, DesiredIndex );

  m_ROI.SetSize( DesiredSize );
  m_ROI.SetIndex( DesiredIndex );


  // adjust ROI so that it lies within the fixed image
  if( m_ROI.Crop( m_FixedImage->GetLargestPossibleRegion() ) ) {
    vcl_cout << "Making sure inverse mapped ROI lies within the moving image" << vcl_endl;
    itk::ImageRegion<3>::SizeType  DesiredSize = m_ROI.GetSize();
    itk::ImageRegion<3>::IndexType DesiredIndex = m_ROI.GetIndex();
    vcl_cout << "ROI: " << DesiredIndex[0] << ".." << DesiredIndex[0]+DesiredSize[0] << ", "
                        << DesiredIndex[1] << ".." << DesiredIndex[1]+DesiredSize[1] << ", "
                        << DesiredIndex[2] << ".." << DesiredIndex[2]+DesiredSize[2] << vcl_endl;

  }

  fixed_x0 = ROIOrigin.GetVnlVector();
  fixed_x1[0] = fixed_x0[0] + m_ROISize[0];
  fixed_x1[1] = fixed_x0[1] + m_ROISize[1];
  fixed_x1[2] = fixed_x0[2] + m_ROISize[2];
  }



  // Set up Covariance Driven Correspondence object
  //
  vcl_vector< double > fixed_spacing;
  fixed_spacing.push_back( 1.0 );

#if SYMMETRIC_MATCHING
	// if segmented volumes supplied
  if( m_FixedImageSegmented && m_MovingImageSegmented ) {  // use oversegmentations to get features
    std::cout << "Fixed and moving segmented images loaded, so using them to get features in segmented regions" << std::endl;
//    vcl_vector< feature_type::sptr >  moving_inside;
//    vcl_vector< feature_type::sptr >  fixed_inside;

    this->FeaturesInRegions( moving_inside, fixed_inside, moving_x0, moving_x1, fixed_x0, fixed_x1 );
  }
	else {
    this->FeaturesInROIs( moving_inside, fixed_inside, moving_x0, moving_x1, fixed_x0, fixed_x1 );
	}

  if( moving_inside.size() < 1000 || fixed_inside.size() < 1000 ) return false;

  // if already ran and initialized with fixed data, set only new moving data from the different ROI
  if( estimation ) {
    estimation->set_features( moving_inside, fixed_inside );
    estimation->set_transform( transform );
  }
  else {
    estimation = new estimation_symmetric_type( m_FeaturesMovingAll[0], m_FeaturesFixedAll[0], transform, transform->inverse(), moving_inside, fixed_inside );
    estimation->weight_by_strength( true );
    estimation->normalize_matches();
    estimation->set_voronoi_maps( m_MovingVoronoiMap, m_FixedVoronoiMap );
  }
#elif CDC
  this->GetFeaturesInMovingROI();

  // FORWARD, full covariance
  // FIX: Unncecessary copying
  vcl_vector< cdcl_feature<dim>::sptr >  m_FeaturesFixedAllWOutShape;
  for( vcl_vector< feature_sptr_type >::size_type  i = 0; i < m_FeaturesFixedAll[0].size(); ++i ) {
    feature_sptr_type  p = m_FeaturesFixedAll[0][i];
    m_FeaturesFixedAllWOutShape.push_back( p->clone() );
  }
  vcl_vector< vcl_vector< feature_sptr_type > >  m_FeaturesFixedAll0;
  m_FeaturesFixedAll0.push_back( m_FeaturesFixedAllWOutShape );

  vcl_vector< feature_sptr_type >  m_FeaturesMovingROIWoutShape;
  for( vcl_vector< feature_sptr_type >::size_type  i = 0; i < m_FeaturesMovingROI[0].size(); ++i ) {
    feature_sptr_type  p = m_FeaturesMovingROI[0][i];
    m_FeaturesMovingROIWoutShape.push_back( p->clone() );
  }
  vcl_vector< vcl_vector< cdcl_feature< 3 >::sptr > >  m_FeaturesMovingROI0;
  m_FeaturesMovingROI0.push_back( m_FeaturesMovingROIWoutShape );


  // get features from the fixed ROI
  vcl_vector< cdcl_feature<dim>::sptr >  m_FeaturesFixedROIWOutShape;
  for( vcl_vector< feature_sptr_type >::size_type  i = 0; i < m_FeaturesFixedAll[0].size(); ++i ) {
    // moving point p
    feature_sptr_type  p = m_FeaturesFixedAll[0][i];

    // make sure current point is inside region of interest
    itk::Point< double, 3 >  p_itk( p->location_.data_block() );
    itk::Index< 3 >  point_index;
    m_FixedImage->TransformPhysicalPointToIndex( p_itk, point_index );

    if( m_ROI.IsInside( point_index ) ) {
      m_FeaturesFixedROIWOutShape.push_back( p->clone() );
    }

  }
  vcl_vector< vcl_vector< cdcl_feature< 3 >::sptr > >  m_FeaturesFixedROI0;
  m_FeaturesFixedROI0.push_back( m_FeaturesFixedROIWOutShape );


  //estimation = new cdcl_estimation<dim, dof>( m_FeaturesMovingROI0, m_FeaturesFixedAll0, transform, fixed_spacing );
  estimation = new cdcl_estimation<dim, dof>( m_FeaturesMovingROI0, m_FeaturesFixedROI0, transform, fixed_spacing );
  estimation->normalize_matches();
#else
  this->GetFeaturesInMovingROI();

  if( estimation ) {
    estimation->set_moving_data( m_FeaturesMovingROI );
    estimation->set_transform( transform );
  }
  else {
    estimation = new estimation_type( m_FeaturesMovingROI, m_FeaturesFixedAll, transform, fixed_spacing );
    estimation->weight_by_strength( true );
    estimation->normalize_matches();
  }
#endif

  estimation->initialize();

  return true;
}


void LocationRegistration::GetFeaturesInMovingROI()
{
  // Get moving features in the region of interest
  //
  m_FeaturesMovingROI.clear();
  vcl_vector< feature_type::sptr >  one_level;
  itk::Point< double, 3 >  p_itk;
  itk::Index< 3 >  point_index;
  for( vcl_vector< feature_sptr_type >::size_type  i = 0; i < m_FeaturesMovingAll[0].size(); ++i ) {
    // moving point p
    feature_sptr_type  p = m_FeaturesMovingAll[0][i];

    // make sure current point is inside region of interest
    p_itk[0] = p->location_[0];
    p_itk[1] = p->location_[1];
    p_itk[2] = p->location_[2];
    m_MovingImage->TransformPhysicalPointToIndex( p_itk, point_index );

    if( m_ROImoving.IsInside( point_index ) ) {
      feature_type::sptr  feature = dynamic_cast< feature_type* >( p->clone().as_pointer() );
      one_level.push_back( feature );
    }

  }
  m_FeaturesMovingROI.push_back( one_level );
}   


void LocationRegistration::WarpVolume()
{
  // WarpVolume
  //
  // Setup transform
  //
  // In ITK resampling is done from the output image to the input, so we need to inverse the transform
  TransformType::Pointer  finalTransformInverse = TransformType::New();
  TransformType::CenterType  finalTransformInverseCenter = m_FinalTransform->TransformPoint( m_FinalTransform->GetCenter() );
  finalTransformInverse->SetCenter( finalTransformInverseCenter );
  m_FinalTransform->GetInverse( finalTransformInverse );


  //typedef itk::ResampleImageFilter< ImageType, OutputImageType >    ResampleFilterType;
  typedef itk::ResampleImageFilter< ImageType, ImageType >    ResampleFilterType;
  ResampleFilterType::Pointer  resample = ResampleFilterType::New();

  // If we have initialized (i.e. have region of interest m_ROI), then only warp image in m_ROI
  typedef itk::RegionOfInterestImageFilter< ImageType, ImageType >  RegionInterestType;
  RegionInterestType::Pointer  RegionInterest;


  // Update m_ROImoving
  //
  itk::Vector<unsigned int>  m_ROISizeVector( m_ROISize );
  TransformType::OutputVectorType  mappedROISizeVector = finalTransformInverse->TransformVector( m_ROISizeVector );
  unsigned int  mappedROISize[3];
  for( unsigned int d = 0; d < 3; ++d ) mappedROISize[d] = (unsigned int) mappedROISizeVector[d];

  // convert from physical (size) to image coordinates
  //
  itk::ImageRegion<3>::SizeType  DesiredSize;
  // compute origin of the ROI given size and selected point
  itk::Point< double, 3 >  ROIOrigin;
  ImageType::SpacingType  ImageSpacing = m_MovingImage->GetSpacing();
  for( unsigned int d = 0; d < 3; ++d ) {
    ROIOrigin[d] = queryPointMapped[d] - mappedROISize[d] / 2;
    DesiredSize[d] = ImageType::SizeValueType( mappedROISize[d] / ImageSpacing[d] );
  }
  // convert ROI origin (in physical coordinates) to the corresponding index in the image
  itk::ImageRegion<3>::IndexType  DesiredIndex;
  m_MovingImage->TransformPhysicalPointToIndex( ROIOrigin, DesiredIndex );

  m_ROImoving.SetSize( DesiredSize );
  m_ROImoving.SetIndex( DesiredIndex );

  //vcl_cout << "Index: " << m_ROImoving.GetIndex() << " Size: " << m_ROImoving.GetSize() << vcl_endl;

  // adjust ROI so that it lies within the moving image
  if( m_ROImoving.Crop( m_MovingImage->GetLargestPossibleRegion() ) ) {
    //vcl_cout << "Making sure inverse mapped ROI lies within the moving image" << vcl_endl;
    itk::ImageRegion<3>::SizeType  DesiredSize = m_ROImoving.GetSize();
    itk::ImageRegion<3>::IndexType DesiredIndex = m_ROImoving.GetIndex();
    //vcl_cout << "ROI: " << DesiredIndex[0] << ".." << DesiredIndex[0]+DesiredSize[0] << ", "
    //                    << DesiredIndex[1] << ".." << DesiredIndex[1]+DesiredSize[1] << ", "
    //                    << DesiredIndex[2] << ".." << DesiredIndex[2]+DesiredSize[2] << vcl_endl;
  }
  
  #if ( WRITE_PANELS || LAST_ITER )
  RegionInterest = RegionInterestType::New();
  RegionInterest->SetRegionOfInterest( m_ROImoving );
  RegionInterest->SetInput( m_MovingImage );
  //resample->SetInput( m_MovingImage );
  resample->SetInput( RegionInterest->GetOutput() );


  resample->SetTransform( finalTransformInverse );
  resample->UseReferenceImageOn();
  resample->SetOutputParametersFromImage( m_FixedImage );
  ImageType::PixelType  defaultPixelValue = 0;
  resample->SetDefaultPixelValue( defaultPixelValue );
  resample->Update();

  m_MovingMappedImage = resample->GetOutput();
  #endif
}


void LocationRegistration::ComputeDecisionMeasurements( std::vector< double > &  measurements )
{
#if CDC
  double max_weighted_error = 1.0;
  double max_maxTrace       = 0.02;
  double max_sheet_angles   = 10.0;
  double max_tube_angles    = 10.0;
  double min_weighted_error = 1.0;
  double min_maxTrace       = 0.02;
  double min_sheet_angles   = 10.0;
  double min_tube_angles    = 10.0;
#else
  double sheet_angles = estimation->sheet_angles() * 180.0 / vnl_math::pi;
  double tube_angles = estimation->tube_angles() * 180.0 / vnl_math::pi;
  double weighted_error = estimation->weighted_error();

  vnl_matrix_fixed< double, 3, 3 >  covarianceJ_from_samples( 0.0 );
  double maxTrace = 0.0;
  double maxEvalue = 0.0;
  ComputeTransferErrorCovariance( m_MovingImage, m_ROImoving, estimation->get_transform(), covarianceJ_from_samples, maxTrace, maxEvalue );

  vnl_matrix_fixed< double, 3, 3 >  covarianceJ_from_samples_backward( 0.0 );
  double maxTraceBackward = 0.0;
  double maxEvalueBackward = 0.0;
  estimation->estimate_LS_backward( false );
  ComputeTransferErrorCovariance( m_FixedImage, m_ROI, estimation->get_transform_backward(), covarianceJ_from_samples_backward, maxTraceBackward, maxEvalueBackward );

  double sheet_angles_backward = estimation->sheet_angles_backward() * 180.0 / vnl_math::pi;
  double tube_angles_backward = estimation->tube_angles_backward() * 180.0 / vnl_math::pi;
  double weighted_error_backward = estimation->weighted_error_backward();

  double max_weighted_error = std::max( weighted_error, weighted_error_backward );
  double max_maxTrace       = std::max( maxTrace, maxTraceBackward );
  double max_sheet_angles   = std::max( sheet_angles, sheet_angles_backward );
  double max_tube_angles    = std::max( tube_angles, tube_angles_backward );
  double min_weighted_error = std::min( weighted_error, weighted_error_backward );
  double min_maxTrace       = std::min( maxTrace, maxTraceBackward );
  double min_sheet_angles   = std::min( sheet_angles, sheet_angles_backward );
  double min_tube_angles    = std::min( tube_angles, tube_angles_backward );
#endif

  measurements.clear();
  measurements.resize( 8 );
  measurements[0] = max_weighted_error;
  measurements[1] = max_maxTrace;
  measurements[2] = max_sheet_angles;
  measurements[3] = max_tube_angles;
  measurements[4] = min_weighted_error;
  measurements[5] = min_maxTrace;
  measurements[6] = min_sheet_angles;
  measurements[7] = min_tube_angles;

  std::cout << "Measurements: ";
  for( unsigned int i = 0; i < 8; ++i ) std::cout << measurements[i] << " ";
  std::cout << std::endl;
}


void LocationRegistration::WriteResults( cdcl_trans<dim,dof>::sptr const &  initial_transform,
                                         double const &  initialRMS,
                                         bool converged,
                                         unsigned int k, // iteration
                                         unsigned int candidateNum,
                                         std::vector< double > const &  measurements,
                                         std::string  prefix )
{
  // we need to use inverse, since we are mapping from the fixed coordinate system to the moving coordinate system
  // (we have a deformation vector (forming the field) for each fixed point on a grid
  cdcl_trans<dim,dof>::sptr  transform_inverse = transform->inverse();

  //std::cout << "CDC trans inverse: " << std::endl;
  //transform_inverse->print( std::cout );

  vnl_vector_fixed< double, 3 >  t_centered = transform_inverse->get_translation() \
                                            + transform_inverse->get_A() * fixed_keypoint->location_ \
                                            - transform_inverse->get_A() * transform_inverse->center_moving_ \
                                            - moving_keypoint->location_;

  // do the same for the initial transform
  //
  cdcl_trans<dim,dof>::sptr  initial_transform_inverse = initial_transform->inverse();

  //std::cout << "CDC trans inverse: " << std::endl;
  //transform_inverse->print( std::cout );

  vnl_vector_fixed< double, 3 >  initial_t_centered = initial_transform_inverse->get_translation() \
                                                    + initial_transform_inverse->get_A() * fixed_keypoint->location_ \
                                                    - initial_transform_inverse->get_A() * initial_transform_inverse->center_moving_ \
                                                    - moving_keypoint->location_;


  vcl_ostringstream  resultsName;
  resultsName << locRegDir << "//" << prefix << "results" << vcl_setw( 6 ) << vcl_setfill( '0' ) << k << "_" << vcl_setw( 2 ) << candidateNum << ".txt";
  std::ofstream  outStream( resultsName.str().c_str() );


  if( deformationField ) {
    // The method is a least squares fit: Xa = b,
    // where _a_ are the affine parameters X and _b_ composed using moving and fixed points.
    vnl_matrix_fixed< double, Dimension, Dimension >  ALS;
    vnl_vector_fixed< double, Dimension >  tLS;
    this->FitAffineTransformToDeformationField( deformationField, m_ROI, fixed_keypoint->location_, moving_keypoint->location_, ALS, tLS );

    double medianDeformationComponentLSError;
    double scaleDeformationComponentLSError;
    this->MeanAndScaleOfAffineApproximationError( deformationField, m_ROI, fixed_keypoint->location_, moving_keypoint->location_, ALS, tLS, medianDeformationComponentLSError, scaleDeformationComponentLSError );

    outStream << "Affine Approximation Squared Error: " << medianDeformationComponentLSError << " scale: " << scaleDeformationComponentLSError << std::endl;

    outStream << std::endl << "Affine transform fitted to deformation: " << std::endl;
    outStream << ALS << std::endl << tLS << std::endl << std::endl;


    double medianSquaredInitialCDCError;
    double scaleInitialCDCError;
    this->MeanAndScaleOfAffineApproximationError( deformationField, m_ROI, fixed_keypoint->location_, moving_keypoint->location_, initial_transform_inverse->get_A(), t_centered, medianSquaredInitialCDCError, scaleInitialCDCError );

    outStream << "CDC Initial Error: " << std::sqrt( medianSquaredInitialCDCError ) << " scale: " << std::sqrt( scaleInitialCDCError ) << std::endl;


    double medianSquaredFinalCDCError;
    double scaleFinalCDCError;
    this->MeanAndScaleOfAffineApproximationError( deformationField, m_ROI, fixed_keypoint->location_, moving_keypoint->location_, transform_inverse->get_A(), t_centered, medianSquaredFinalCDCError, scaleFinalCDCError );

    outStream << "CDC Final Error: " << std::sqrt( medianSquaredFinalCDCError ) << " scale: " << std::sqrt( scaleFinalCDCError ) << std::endl;

    outStream << "CDC Error reduced: " << ( (medianSquaredInitialCDCError > medianSquaredFinalCDCError) ? 1 : 0 ) << std::endl;
  }


  outStream << std::endl << "Final estimated transform: " << std::endl;
  outStream << estimation->get_transform()->get_A() << std::endl << t_centered << std::endl << std::endl;

  outStream << "Final estimated covariance: " << std::endl;
  outStream << estimation->get_transform()->get_covariance() << std::endl << std::endl;

  estimation->compute_transfer_error_covariance();
  outStream << "Final estimated covarianceJ from moving points: " << std::endl;
  outStream << estimation->get_transform()->get_covarianceJ() << std::endl;

  const vnl_matrix< double >  covarianceJ_from_pts = estimation->get_transform()->get_covarianceJ();
  vnl_matrix< double >  V( 3, 3 );
  vnl_vector< double >  D( 3 );
  vnl_symmetric_eigensystem_compute( covarianceJ_from_pts, V, D );
  outStream << "covarianceJ from points evalues: " << D( 0 ) << "  " << D( 1 ) << "  " << D( 2 ) << std::endl << std::endl;

  //vnl_symmetric_eigensystem_compute( covarianceJ_from_samples, V, D );
  //outStream << "covarianceJ from samples evalues: " << D( 0 ) << "  " << D( 1 ) << "  " << D( 2 ) << std::endl;

  // RMS is computed using features
  double finalRMS = estimation->RMS_error();

  outStream << "Fixed keypoint location: " << fixed_keypoint->location_ << std::endl;
  outStream << "Moving keypoint location: " << moving_keypoint->location_ << std::endl;

  outStream << "Converged: " << (converged ? 1 : 0) << std::endl;
  outStream << "Initial RMS: " << initialRMS << std::endl;
  outStream << "Final RMS: " << finalRMS << std::endl;
  outStream << "RMS reduced: " << ( (initialRMS > finalRMS) ? 1 : 0 ) << std::endl;
  //outStream << "Maximum largest evalue: " << maxEvalue << std::endl;


  // map query point
  {
  TransformType::Pointer  m_InverseFinalTransform = TransformType::New();
  TransformType::CenterType  TransformInverseCenter = m_FinalTransform->TransformPoint( m_FinalTransform->GetCenter() );
  m_InverseFinalTransform->SetCenter( TransformInverseCenter );
  m_FinalTransform->GetInverse( m_InverseFinalTransform );

  queryPointMapped = m_InverseFinalTransform->TransformPoint( queryPoint );
  outStream << "Query point: " << queryPoint << std::endl;
  outStream << "Query point mapped: " << queryPointMapped << std::endl;

  if( deformationField ) {
    // get the deformation vector from the field
    DeformationFieldType::IndexType  queryIndex;
    deformationField->TransformPhysicalPointToIndex( queryPoint, queryIndex );
    VectorPixelType  deformationVector = deformationField->GetPixel( queryIndex );

    itk::Point< double, 3 >  queryPointMappedDiff = queryPoint + deformationVector;

    vnl_vector_fixed< double, 3 >  distance = queryPointMappedDiff.GetVnlVector() - queryPointMapped.GetVnlVector();
    outStream << "Query point deformed: " << queryPointMappedDiff << std::endl;
    outStream << "Distance between mapped and deformed: " << distance.magnitude() << std::endl;
  }

  // save final estimated transform
  //
  vcl_ostringstream  transform_file_name;
  transform_file_name << locRegDir << "//transform" << vcl_setw( 6 ) << vcl_setfill( '0' ) << k << "_" << vcl_setw( 2 ) << candidateNum << ".vtk";
  itk::TransformFileWriter::Pointer  transformWriter;
  transformWriter = itk::TransformFileWriter::New();

  transformWriter->SetInput( m_FinalTransform );
  transformWriter->SetFileName( transform_file_name.str().c_str() );
  transformWriter->Update();

  // save final estimated inverse transform
  //
  vcl_ostringstream  transform_inverse_file_name;
  transform_inverse_file_name << locRegDir << "//transform_inverse" << vcl_setw( 6 ) << vcl_setfill( '0' ) << k << "_" << vcl_setw( 2 ) << candidateNum << ".vtk";
  transformWriter->SetInput( m_InverseFinalTransform );
  transformWriter->SetFileName( transform_inverse_file_name.str().c_str() );
  transformWriter->Update();
  }


  outStream << "Final alignment error : Maximum trace : Sheet angles (deg): Tube angles (deg)" << std::endl;
  outStream << "Classification:\t";
  for( unsigned int ind = 0; ind < 8; ++ind ) outStream << measurements[ind] << "\t";
  outStream << std::endl;


  outStream.close();
}


int LocationRegistration::Run()
{
  this->ReadImages();

  this->ReadSegmentedImages();
  
  this->ReadFeatures();

  int returnStatus = this->ReadVoronoiMaps();
  // returnStatus == 0  if success
  if( returnStatus ) return returnStatus;

  this->ReadMovingKeypointDescriptors();

  this->ReadFixedKeypointDescriptors();

  int cantReadGT = this->SetupGroundTruthTransform() ;
  if( cantReadGT ) return cantReadGT;



  unsigned int k = 0;

  typedef itk::RegularExpressionSeriesFileNames    NameGeneratorType;

  NameGeneratorType::Pointer nameGenerator = NameGeneratorType::New();

  // find out how many initializations we have
  //
  const unsigned int subMatch = 0;
  std::string  regularExpression = "fixed_descriptor.*_00.*";
  nameGenerator->SetRegularExpression( regularExpression );
  nameGenerator->SetSubMatch( subMatch );
  nameGenerator->SetDirectory( matchesDirPath );
  unsigned int fixedDescFilesSize = nameGenerator->GetFileNames().size();

  regularExpression = "moving_descriptor.*_00.*";
  nameGenerator->SetRegularExpression( regularExpression );
  unsigned int movingDescFilesSize = nameGenerator->GetFileNames().size();

  if( fixedDescFilesSize == 0 ) {
    std::cout << "No keypoint matches found." << std::endl;
    return 1;
  }
  if( movingDescFilesSize == 0 ) {
    std::cout << "No keypoint matches found." << std::endl;
    return 1;
  }
  if( fixedDescFilesSize != movingDescFilesSize ) {
    std::cout << "Wrong keypoint matches: Different number of fixed and moving descriptors." << std::endl;
    return 1;
  }


  // verification
  struct svm_model* model;
  // model file is assumed to be in the curent directory
  char* model_file = "training_data_all.model";
  bool predict_probability = 1; 
  load_model( model, model_file, predict_probability );

  // measurment vector size
  const unsigned int size = 8;

  struct svm_node *x = new svm_node[size+1];  // need one more for the end of array
  for( unsigned int ind = 0; ind < size; ++ind ) {
    x[ind].index = ind;
    x[ind].value = 0.0;
  }

  // end of array
  x[8].index = -1;
  x[8].value = 0.0;

  // copied from scaling_factors saved by svm-scale (minimum and maximum possible values of measurements)
  double feature_min[size] = {0.823246, 0.17355, 13.5656, 14.8583, 0.776814, 0.124, 10.3695, 13.049};
  double feature_max[size] = {16.1025, 139.276, 62.7776, 63.5239, 1.64664, 60.897, 59.8739, 60.8141};


#define TEST_SOME 0

#if !TEST_SOME

  unsigned int numInitializations = 0;
  std::vector< std::string >  transformFiles;

  if( DECISION_DEBUGGING ) {
    typedef itk::RegularExpressionSeriesFileNames    NameGeneratorType;

    NameGeneratorType::Pointer transformNameGenerator = NameGeneratorType::New();

    std::string  transformRegularExpression = "transform*.*";
    transformNameGenerator->SetRegularExpression( transformRegularExpression );
    transformNameGenerator->SetSubMatch( subMatch );
    transformNameGenerator->SetDirectory( locRegDir );

    transformFiles = transformNameGenerator->GetFileNames();

    numInitializations = transformFiles.size();
  }
  else {
    numInitializations = std::min( fixedDescFilesSize, 100u ); // 100 matches
  }

//k=41;
  while( k < numInitializations) {

#else
  const unsigned int numInitializations = 8;
  const unsigned int testThese[numInitializations] = {63, 13, 36, 49, 52, 76, 80, 87};
  unsigned int q = 0;

  while( q < numInitializations) {
    unsigned int k = testThese[q];
#endif


    itk::TimeProbe timer;
    timer.Start();

    std::vector< std::string >  movingDescFiles;
    std::vector< std::string >  fixedDescFiles;

    unsigned int candidateNum = 0;
    unsigned int clickNum = 0;
    if( DECISION_DEBUGGING ) {
      // given: transform000010_00.vtk will return: 000010_00
      std::string  initNumber = itksys::SystemTools::GetFilenameName( transformFiles[k] ).substr( 9, 9 );

      vcl_istringstream  initNumberStream;
      initNumberStream.str( initNumber );
      char underscore = '_';
      initNumberStream >> clickNum >> underscore >> candidateNum;

      std::stringstream  nextTransformStream;
      nextTransformStream << "transform" << vcl_setw( 6 ) << vcl_setfill( '0' ) << clickNum << "_" << vcl_setw( 2 ) << candidateNum+1 << ".vtk";

      std::string  nextTransformPath = locRegDir + "/" + nextTransformStream.str();

      std::cout << "NEXT: " << nextTransformPath << std::endl;

      //if( itksys::SystemTools::FileExists( nextTransformPath.c_str() ) ) {
      //  ++k;
      //  continue;
      //}
    
      std::string  movingDescFile = matchesDirPath + "/moving_descriptor" + initNumber + ".vtk";
      std::string  fixedDescFile  = matchesDirPath + "/fixed_descriptor" + initNumber + ".vtk";

      movingDescFiles.push_back( movingDescFile );
      fixedDescFiles.push_back( fixedDescFile );
    }
    else {
      clickNum = k;
      // get filenames of candidate initializations for this location
      //
      const unsigned int subMatch = 0;
      std::stringstream  fixedRegularExpressionStream;
      fixedRegularExpressionStream << "fixed_descriptor" << std::setfill( '0' ) << std::setw( 6 ) << k << ".*";
      nameGenerator->SetRegularExpression( fixedRegularExpressionStream.str() );
      nameGenerator->SetSubMatch( subMatch );
      nameGenerator->SetDirectory( matchesDirPath );

      fixedDescFiles = nameGenerator->GetFileNames();

      std::stringstream  movingRegularExpressionStream;
      movingRegularExpressionStream << "moving_descriptor" << std::setfill( '0' ) << std::setw( 6 ) << k << ".*";
      nameGenerator->SetRegularExpression( movingRegularExpressionStream.str() );

      movingDescFiles = nameGenerator->GetFileNames();

      if( fixedDescFiles.size() == 0 ) {
        std::cout << "No keypoint matches found for this location." << std::endl;
        return 1;
      }
      if( movingDescFiles.size() == 0 ) {
        std::cout << "No keypoint matches found for this location." << std::endl;
        return 1;
      }
      if( fixedDescFiles.size() != movingDescFiles.size() ) {
        std::cout << "Wrong keypoint matches: Different number of fixed and moving descriptors for this location." << std::endl;
        return 1;
      }
    }

    vcl_vector< cdcl_keypoint<3>::sptr >  moving_keypoints;
    vcl_vector< cdcl_keypoint<3>::sptr >  fixed_keypoints;
  
    vcl_vector< vnl_vector< float > >  moving_descriptors;
    vcl_vector< vnl_vector< float > >  fixed_descriptors;

    vcl_vector< itk::Point< double, 3 > >  queryPoints;


    typedef std::pair< double, unsigned int >  CandidateIndType;
    std::vector< CandidateIndType >  candidates;
    unsigned int candidateInd = 0;
    while( candidateInd < movingDescFiles.size() ) {
      if( !DECISION_DEBUGGING ) candidateNum = candidateInd;
    
      std::cout << "match files: " << movingDescFiles[candidateInd] << "  " << fixedDescFiles[candidateInd] << std::endl;

      // read the moving descriptor
      vcl_vector< cdcl_keypoint<3>::sptr >  moving_keypoint_vect;
      vcl_vector< vnl_vector< float > >  moving_descriptor_vect;
      cdcl_read_keypoint_descriptors_VTK( moving_keypoint_vect, moving_descriptor_vect, movingDescFiles[candidateInd] );
      moving_keypoint = moving_keypoint_vect[0];
      vnl_vector< float >  moving_descriptor = moving_descriptor_vect[0];
      moving_keypoints.push_back( moving_keypoint );
      moving_descriptors.push_back( moving_descriptor );

      // read the fixed descriptor
      vcl_vector< cdcl_keypoint<3>::sptr >  fixed_keypoint_vect;
      vcl_vector< vnl_vector< float > >  fixed_descriptor_vect;
      cdcl_read_keypoint_descriptors_VTK( fixed_keypoint_vect, fixed_descriptor_vect, fixedDescFiles[candidateInd] );
      fixed_keypoint = fixed_keypoint_vect[0];
      vnl_vector< float >  fixed_descriptor = fixed_descriptor_vect[0];
      fixed_keypoints.push_back( fixed_keypoint );
      fixed_descriptors.push_back( fixed_descriptor );

      // read clicked query location
      vcl_string  queryLocationFileName = fixedDescFiles[candidateInd];
      itksys::SystemTools::ReplaceString( queryLocationFileName, "fixed_descriptor", "query_location" );
      itksys::SystemTools::ReplaceString( queryLocationFileName, ".vtk", ".txt" );
      vcl_ifstream  queryLocationFile( queryLocationFileName.c_str() );
      queryLocationFile >> queryPoint[0] >> queryPoint[1] >> queryPoint[2];
      queryPoints.push_back( queryPoint );
      queryLocationFile.close();


      // initialize transform and GetCDCTransform from ITK transform
      if( DECISION_DEBUGGING ) {
        this->ReadFinalTransform( transformFiles[k] );
      }
      else {
        this->SetupFinalTransform();
      }
      transform = this->GetCDCTransform();
      cdcl_trans<dim,dof>::sptr  initial_transform = transform->clone();


      //const unsigned int sizes[2] = {50, 70};
      const unsigned int sizes[1] = {50};
      //const unsigned int sizes[3] = {30, 50, 70};
      //for( unsigned int s = 0; s < 2; ++s ) {

      bool estimationSetup = false;
      unsigned int s = 0;
      while( !estimationSetup && s < 1 ) {

        m_ROISize[0] = sizes[s];
        m_ROISize[1] = sizes[s];
        m_ROISize[2] = sizes[s];

        estimationSetup = this->SetupEstimation();

        // if it was not possible to setup the estimation, try with larger region
//        if( !this->SetupEstimation() ) {
//          ++candidateInd;
//          continue;
//        }

        ++s;
      }

      if( !estimationSetup ) {
        ++candidateInd;
        continue;
      }

      estimation->find_closest_euclidean();

      //estimation->estimate_scale_and_assign_weight();

      double sheet_angles = estimation->sheet_angles() * 180.0 / vnl_math::pi;
      double tube_angles = estimation->tube_angles() * 180.0 / vnl_math::pi;

      double sheet_angles_backward = estimation->sheet_angles_backward() * 180.0 / vnl_math::pi;
      double tube_angles_backward = estimation->tube_angles_backward() * 180.0 / vnl_math::pi;

      double angles_average = ( sheet_angles + tube_angles + sheet_angles_backward + tube_angles_backward ) / 4;

      std::cout << "angles: " << candidateNum << ": " << sheet_angles << "  " << tube_angles << "  "
                                                    << sheet_angles_backward << "  " << tube_angles_backward << "  "
                                                    << angles_average << std::endl;

      candidates.push_back( CandidateIndType( angles_average, candidateNum ) );

      ++candidateInd;
//      }

    }

    std::sort( candidates.begin(), candidates.end(), AnglesSmaller );

    bool verified = false;
    candidateInd = 0;
//candidateNum=4;
//candidateNum=19;
    while( !verified && candidateInd < candidates.size() ) {
    
      if( DECISION_DEBUGGING ) candidateNum = 0; // temporarily set the candidateNum to 0 to do the indexing below right 
		  else candidateNum = candidates[candidateInd].second;

      moving_keypoint = moving_keypoints[candidateNum];
      vnl_vector< float >  moving_descriptor = moving_descriptors[candidateNum];

      fixed_keypoint = fixed_keypoints[candidateNum];
      vnl_vector< float >  fixed_descriptor = fixed_descriptors[candidateNum];

      queryPoint = queryPoints[candidateNum];

      candidateNum = candidates[candidateInd].second;

      std::cout << "Processing initialization: " << k << " candidate num: " << candidateNum << std::endl;



      if( DECISION_DEBUGGING ) {
        this->ReadFinalTransform( transformFiles[k] );
      }
      else {
        this->SetupFinalTransform();
      }

      // initialize transform and GetCDCTransform from ITK transform
      transform = this->GetCDCTransform();
      cdcl_trans<dim,dof>::sptr  initial_transform = transform->clone();


      // RMS error after initialization
      double initialRMS = -1.0;


      double probability_being_good = 0.0;
      //const unsigned int sizes[2] = {50, 70};
      const unsigned int sizes[1] = {50};
      //const unsigned int sizes[3] = {30, 50, 70};
//      for( unsigned int s = 0; s < 2; ++s ) {


      bool estimationSetup = false;
      unsigned int s = 0;
      while( !estimationSetup && s < 1 ) {

        m_ROISize[0] = sizes[s];
        m_ROISize[1] = sizes[s];
        m_ROISize[2] = sizes[s];

        estimationSetup = this->SetupEstimation();

        // if it was not possible to setup the estimation, try with larger region
//        if( !this->SetupEstimation() ) {
//          ++candidateInd;
//          continue;
//        }

        ++s;
      }

      if( !estimationSetup ) {
	      std::cout << "Warning: estimation not setup for candidate: " << candidateNum << " but it should have been." << std::endl;
        continue;
      }


        #if ( WRITE_PANELS || LAST_ITER )
        cdcl_extract_data_callback3d< dim, dof, feature_type >  extractDataCallback;
        //extractDataCallback.initialize_moving_mapped_features( m_FeaturesMovingROI[0], transform );
        extractDataCallback.initialize_moving_mapped_features( moving_inside, transform );

        vtkSmartPointer< vtkPolyData >  polyDataFixed = vtkSmartPointer< vtkPolyData >::New();
        cdcl_features_to_poly_data( m_FeaturesFixedAll[0], polyDataFixed );
//cdcl_features_to_poly_data( fixed_inside, polyDataFixed ); -- this does not display features, why? (would be more efficient though)
        #endif

        
        unsigned int max_iterations;
        if( DECISION_DEBUGGING ) max_iterations = 1;
        else max_iterations = 70;

        unsigned int iteration = 0;

        bool converged = false;
        //while( iteration < max_iterations ) {
	    while( !converged && iteration < max_iterations ) {

          #if WRITE_PANELS
          #if CDC
          converged = estimation->one_iteration( (void*) &extractDataCallback, &cdcl_extract_data_callback3d< dim, dof >::display_points, iteration );  
          #else
          converged = estimation->one_iteration( (void*) &extractDataCallback, &cdcl_extract_data_callback3d< dim, dof, feature_type >::display_points, iteration );  
          #endif
          if( iteration == 0 && initialRMS == -1.0 ) initialRMS = estimation->RMS_error();

          //// write only only first iteration, when converged or the last iteration
          //if( !converged && iteration > 0 && ( iteration % 10 != 0 ) ) {
          //  ++iteration;
          //  continue;
          //}

          this->WarpVolume();

          vtkSmartPointer< vtkPolyData >  polyDataMovingMapped = extractDataCallback.GetPolyDataMovingROIMapped();

          OutputImageType::Pointer  m_SlicesImage;
          this->GenerateSlices( polyDataMovingMapped, polyDataFixed, m_MovingMappedImage, m_FixedImage, m_ROI, m_SlicesImage );

          #else

          #if CDC
          converged = estimation->one_iteration( 0, 0, iteration );  
          #else
          converged = estimation->one_iteration( 0, 0, iteration );  
          #endif
          if( iteration == 0 && initialRMS == -1.0 ) initialRMS = estimation->RMS_error();

          this->WarpVolume();
          #endif

          // UpdateITKTransform
          // Do this at the end because polydata get written in the callback before the transform is updated.
          // If we put this after estimation, then image gets resampled with the updated transform, but points are still mapped
          // using the previous transform estimate.
          this->UpdateITKTransform();


          #if ( WRITE_PANELS )
          // Write the checkerboard
          //
          typedef itk::ImageFileWriter<OutputImageType>                               WriterType;

          //  Writer
          WriterType::Pointer writer = WriterType::New();
          writer->SetInput( m_SlicesImage );

          vcl_ostringstream  checkerName;
          if( deformationType != NONE ) {
            checkerName << locRegDir << "//checkergood" << vcl_setw( 6 ) << vcl_setfill( '0' ) << k << "_" << vcl_setw( 2 ) << candidateNum << "_" << vcl_setw( 2 ) << iteration << ".png";
          }
          else {
            checkerName << locRegDir << "//checker" << vcl_setw( 6 ) << vcl_setfill( '0' ) << k << "_" << vcl_setw( 2 ) << candidateNum << "_" << vcl_setw( 2 ) << iteration << ".png";
          }
          writer->SetFileName( checkerName.str().c_str() );
          writer->Update();
          #endif

          ++iteration;

          if( !validTransform( estimation->get_transform() ) ) break;
      
        }



//converged=true;
        // no need to grow more if current have not converged
        // this is right when we already have grown a little
        // however, if we haven't grown and first region failed, larger region could still converge
        // it might be better though to try first all initializations and then try them again with larger initialization region
        //if( !DECISION_DEBUGGING && !converged ) break;
        //the condition above is for growing, the below is w/out growing
        if( !DECISION_DEBUGGING && !converged ) {
			    ++candidateInd;
				  continue;
			  }

 #if !CDC
        std::cout << "converged: " << converged << " oscillated: " << estimation->oscillated() << std::endl;
        //if( estimation->oscillated() ) break;
        //the condition above is for growing, the below is w/out growing
        if( estimation->oscillated() ) {
					++candidateInd;
					continue;
				}
#endif
        std::vector< double >  measurements;
        this->ComputeDecisionMeasurements( measurements );

        // testing vector
        for( unsigned int ind = 0; ind < 8; ++ind ) x[ind].value = measurements[ind];

        // run SVM verification
        const unsigned int nr_class = 2;
        double probability_estimates[nr_class];
        double result = predict( model, x, feature_min, feature_max, size, predict_probability, probability_estimates );

        std::cout << "REGION: " << s << " ver: " << result << " prob: " << probability_estimates[0] << " :";
        for( unsigned int ind = 0; ind < 8; ++ind ) std::cout << x[ind].value << " ";
        std::cout << std::endl;

        //if( result == 1.0 ) verified = true;
        // probability of being aligned must be greated than specified value
        //if( probability_estimates[0] > 0.85 ) verified = true;
//probability_estimates[0] = 0.99;
        if( probability_estimates[0] > 0.5 ) verified = true;
        else if( verified ) break;  // verified has been set to true before, but now verification failed, so do not continue with growing
//verified=true;
        std::cout << "verified: " << verified << std::endl;
        // write the results only if probability of being good is greater than in the previous run
        // this will overwrite results that used smaller region
        if( verified ) {
          if( probability_estimates[0] > probability_being_good ) { 
            if( DECISION_DEBUGGING )
              this->WriteResults( initial_transform, initialRMS, converged, clickNum, candidateNum, measurements, "decision_" );
            else
              this->WriteResults( initial_transform, initialRMS, converged, clickNum, candidateNum, measurements );
            probability_being_good = probability_estimates[0];
          }
          else break; // if probability got worse do not grow anymore
        }


//     }


        // if we want only last checkerboard, we need to run one more iteration to generate it, this is the same as in the inner loop above
        #if LAST_ITER
          m_ROISize[0] = sizes[s] + 30;
          m_ROISize[1] = sizes[s] + 30;
          m_ROISize[2] = sizes[s] + 30;

          #if CDC
          converged = estimation->one_iteration( (void*) &extractDataCallback, &cdcl_extract_data_callback3d< dim, dof >::display_points, iteration );  
          #else
          converged = estimation->one_iteration( (void*) &extractDataCallback, &cdcl_extract_data_callback3d< dim, dof, feature_type >::display_points, iteration );  
          #endif

          this->WarpVolume();

          vtkSmartPointer< vtkPolyData >  polyDataMovingMapped = extractDataCallback.GetPolyDataMovingROIMapped();

          OutputImageType::Pointer  m_SlicesImage;
          this->GenerateSlices( polyDataMovingMapped, polyDataFixed, m_MovingMappedImage, m_FixedImage, m_ROI, m_SlicesImage );

          // Write the checkerboard
          //
          typedef itk::ImageFileWriter<OutputImageType>                               WriterType;

          //  Writer
          WriterType::Pointer writer = WriterType::New();
          writer->SetInput( m_SlicesImage );

          vcl_ostringstream  checkerName;
          checkerName << locRegDir << "//checkerfinal" << vcl_setw( 6 ) << vcl_setfill( '0' ) << k << "_" << vcl_setw( 2 ) << candidateNum << "_" << vcl_setw( 2 ) << iteration << ".png";
          writer->SetFileName( checkerName.str().c_str() );
          writer->Update();

        #endif



      ++candidateInd;
    }
    // process another location
    ++k;
    timer.Stop();
    std::cout << "Timer: This location took " << timer.GetMeanTime() << " seconds.\n";
    #if TEST_SOME
    ++q;
    #endif
  }

  delete [] x;
  svm_destroy_model(model);


  return 0;
}


// The method is a least squares fit: Xa = b,
// where _a_ are the affine parameters X and _b_ composed using moving and fixed points.
void LocationRegistration::FitAffineTransformToDeformationField( DeformationFieldType::Pointer const &  deformationField,
                                           itk::ImageRegion< 3 > const &  m_ROI,
                                           vnl_vector_fixed< double, 3 > const &  fixed_keypoint_location,
                                           vnl_vector_fixed< double, 3 > const &  moving_keypoint_location,
                                           vnl_matrix_fixed< double, 3, 3 > &  A,
                                           vnl_vector_fixed< double, 3 > &  t )
{
  const unsigned int Dimension = LocationRegistration::dim;
  const unsigned int Dof = LocationRegistration::dof;

  //std::cout << "fixed: " << fixed_keypoint_location << "  " << moving_keypoint_location << std::endl;

  vnl_matrix_fixed< double, Dof, Dof >  XtX( 0.0 );
  vnl_vector_fixed< double, Dof >  Xtq( 0.0 );

  double num = 0.0;
  typedef itk::ImageRegionConstIterator< DeformationFieldType > FieldIterator;
  FieldIterator  fi( deformationField, m_ROI );


  for( fi.GoToBegin(); !fi.IsAtEnd(); ++fi ) {

    itk::Point< double, Dimension >  sourcePoint;
    deformationField->TransformIndexToPhysicalPoint( fi.GetIndex(), sourcePoint );

    vnl_vector_fixed< double, Dimension >  sourcePointVnl = sourcePoint.Get_vnl_vector();

    DeformationFieldType::PixelType  deformationVector = fi.Value();
    itk::Point< double, Dimension >  deformedPoint = sourcePoint + deformationVector;
    vnl_vector_fixed< double, Dimension >  deformedPointVnl = deformedPoint.Get_vnl_vector();
    

    // Jacobian w.r.t. parameters.
    vnl_matrix_fixed< double, Dimension, Dof >  jacobian_wrt_par_;

    // center the point on the region
    vnl_vector_fixed< double, Dimension >  const  location_centered = sourcePointVnl - fixed_keypoint_location;


    jacobian_wrt_par_.fill( 0.0 );

    // derivative w.r.t. each of the affine matrix parameters
    unsigned int p = 0;
    for( unsigned int d1 = 0; d1 < Dimension; ++d1 ) 
      for( unsigned int d2 = 0; d2 < Dimension; ++d2 ) {
        jacobian_wrt_par_( d1, p ) = location_centered[d2];
        ++p;
      }

    // derivative w.r.t. each of the translation vector parameters
    for( unsigned int d = 0; d < Dimension; ++d, ++p ) {
      jacobian_wrt_par_( d, p ) = 1.0;
    }

    // weight equals 1/N
    XtX += jacobian_wrt_par_.transpose() * jacobian_wrt_par_;

    //std::cout << XtX << std::endl << std::endl;

    // weight equals 1/N
    Xtq += jacobian_wrt_par_.transpose() * ( deformedPointVnl - moving_keypoint_location );

    ++num;
    //std::cout << Xtq << std::endl << std::endl;

  }

  Xtq /= num;
  XtX /= num;

  vnl_svd< double > svd( XtX );
  vnl_matrix_fixed< double, Dof, Dof >  invXtX = svd.inverse();


  //std::cout << "XtX: " << XtX << std::endl << std::endl;
  //std::cout << "XtX^-1: " << invXtX << std::endl << std::endl;
  //std::cout << "Xtq: " << Xtq << std::endl << std::endl;

  vnl_vector_fixed< double, Dof >  Aparams = invXtX * Xtq;

  //std::cout << "LS estimated: " << Aparams << std::endl;

  A.set_row( 0, Aparams.extract( 3, 0 ) );
  A.set_row( 1, Aparams.extract( 3, 3 ) );
  A.set_row( 2, Aparams.extract( 3, 6 ) );

  t = Aparams.extract( 3, 9 );

  std::cout << "LS estimated: " << std::endl << A << std::endl << t << std::endl;
}


void LocationRegistration::MeanAndScaleOfAffineApproximationError( DeformationFieldType::Pointer const &  deformationField,
                                             itk::ImageRegion< 3 > const &  m_ROI,
                                             vnl_vector_fixed< double, 3 > const &  fixed_keypoint_location,
                                             vnl_vector_fixed< double, 3 > const &  moving_keypoint_location,
                                             vnl_matrix_fixed< double, 3, 3 > const &  A,
                                             vnl_vector_fixed< double, 3 > const &  t,
                                             double &  medianDeformationComponentLSError,
                                             double &  scaleDeformationComponentLSError )
{
  const unsigned int Dimension = LocationRegistration::dim;

  // Get center point by interpolating the field at the keypoint location.
  typedef itk::VectorLinearInterpolateImageFunction< DeformationFieldType, double >  FieldInterpolatorType;
  FieldInterpolatorType::Pointer  interpField = FieldInterpolatorType::New();
  interpField->SetInputImage( deformationField );

  FieldInterpolatorType::InputType  fixedPointCenter;
  fixedPointCenter[0] = fixed_keypoint_location[0];
  fixedPointCenter[1] = fixed_keypoint_location[1];
  fixedPointCenter[2] = fixed_keypoint_location[2];
  FieldInterpolatorType::OutputType  deformationAtCenterVector = interpField->Evaluate( fixedPointCenter );

  vnl_vector_fixed< double, Dimension >  deformationAtCenterVectorVnl;
  deformationAtCenterVectorVnl[0] = deformationAtCenterVector[0];
  deformationAtCenterVectorVnl[1] = deformationAtCenterVector[1];
  deformationAtCenterVectorVnl[2] = deformationAtCenterVector[2];


  //// Get center point by using the deformation vector right at the center of the ROI.
  //// fixed_keypoint_location mapped
  //typename DeformationFieldType::IndexType  centerIndex = m_ROI.GetIndex();
  //centerIndex[0] = m_ROI.GetSize()[0] / 2;
  //centerIndex[1] = m_ROI.GetSize()[1] / 2;
  //centerIndex[2] = m_ROI.GetSize()[2] / 2;
  //
  //itk::Point< double, 3 >  centerPoint;
  //deformationField->TransformIndexToPhysicalPoint( centerIndex, centerPoint );
  //vnl_vector_fixed< double, Dimension >  fixed_keypoint_location = centerPoint.Get_vnl_vector();

  //typename DeformationFieldType::PixelType  deformationAtCenterVector = deformationField->GetPixel( centerIndex );

  //vnl_vector_fixed< double, Dimension >  deformationAtCenterVectorVnl;
  //deformationAtCenterVectorVnl[0] = deformationAtCenterVector[0];
  //deformationAtCenterVectorVnl[1] = deformationAtCenterVector[1];
  //deformationAtCenterVectorVnl[2] = deformationAtCenterVector[2];

  //vnl_vector_fixed< double, Dimension >  fixed_keypoint_location_mapped = fixed_keypoint_location + deformationAtCenterVectorVnl;


  //std::cout << "fixed: " << fixed_keypoint_location << "  " << moving_keypoint_location << std::endl;

  typedef itk::ImageRegionConstIterator< DeformationFieldType > FieldIterator;
  FieldIterator  fi( deformationField, m_ROI );

  // Compute different errors caused by using kernel transform as opposed to deformation field
  // and how much is the affine acomponent when using kernel transform and least squres estimated transform.
  std::vector< double >  deformationComponentsLS;
  for( fi.GoToBegin(); !fi.IsAtEnd(); ++fi ) {

    itk::Point< double, 3 >  sourcePoint;
    deformationField->TransformIndexToPhysicalPoint( fi.GetIndex(), sourcePoint );

    vnl_vector_fixed< double, Dimension >  sourcePointVnl = sourcePoint.Get_vnl_vector();

    DeformationFieldType::PixelType  deformationVector = fi.Value();
    vnl_vector_fixed< double, Dimension >  deformationVectorVnl;
    deformationVectorVnl[0] = deformationVector[0];
    deformationVectorVnl[1] = deformationVector[1];
    deformationVectorVnl[2] = deformationVector[2];

    //std::cout << sourcePointVnl << " " << deformationVectorVnl << std::endl;


    vnl_vector_fixed< double, Dimension >  deformationAffineLS = A * (sourcePointVnl - fixed_keypoint_location) + t - sourcePointVnl + moving_keypoint_location;
    // It does not matter if it is centered or not, since we only need relative component (affine transformation) -- the two centers would get subtracted
    //vnl_vector_fixed< double, Dimension >  deformationAffineLS = A * (sourcePointVnl ) + t - sourcePointVnl ;

    //std::cout << sourcePointVnl << "  " << tLS << std::endl;


    double deformationComponentLS = ( deformationVectorVnl - deformationAffineLS ).two_norm();
    deformationComponentsLS.push_back( deformationComponentLS  );

    //std::cout << deformationVectorVnl.two_norm() << " " << deformationAffineLS.two_norm() << " " << deformationComponentLS << std::endl;

  }
  //std::cout << "Original Deformation Shift | Shift from LS affine | Deformation - LS Affine" << std::endl;

  rrel_util_median_and_scale( deformationComponentsLS.begin(), deformationComponentsLS.end(), medianDeformationComponentLSError, scaleDeformationComponentLSError, 1 );

  std::cout << "Sample size: " << deformationComponentsLS.size() << std::endl;

}


void LocationRegistration::ComputeTransferErrorCovariance( itk::Image< signed short, dim >::Pointer const &  Image,
                                     itk::ImageRegion< dim > const &  m_ROImoving,
                                     vbl_smart_ptr< cdcl_trans<dim, dof> > const & transform,
                                     vnl_matrix_fixed< double, dim, dim > &  covarianceJ,
                                     double &  maxTrace,
                                     double &  maxEval )
{
  typedef itk::Image< signed short, dim >  ImageType;

  vnl_matrix_fixed< double, dof, dof >  covariance = transform->get_covariance();

  ImageType::SpacingType  spacing = Image->GetSpacing();
  itk::Vector< double, dim >  size;
  size[0] = m_ROImoving.GetSize()[0] * spacing[0];
  size[1] = m_ROImoving.GetSize()[1] * spacing[1];
  size[2] = m_ROImoving.GetSize()[2] * spacing[2];

  itk::Point< double, dim >  start;
  Image->TransformIndexToPhysicalPoint( m_ROImoving.GetIndex(), start );
  // center the region
  start[0] = start[0];
  start[1] = start[1];
  start[2] = start[2];

  itk::Point< double, dim >  end = start + size - spacing;

  itk::Point< double, dim >  center = start + 0.5 * size;

  //vcl_cout << "Centers: " << transform->center_moving_ << "  " << center << vcl_endl;

  covarianceJ.fill( 0.0 );
  double interval = 5.0;  // 5mm
  unsigned int count = 0;
  maxTrace = 0.0;
  maxEval = 0.0;

  // the region could be sampled only at the boundaries for efficiency reasons
  for( double k = start[2]; k <= end[2]; k += interval ) {
    for( double j = start[1]; j <= end[1]; j += interval ) {
      for( double i = start[0]; i <= end[0]; i += interval ) {
        vnl_vector_fixed< double, 3 >  point;
        point[0] = i;
        point[1] = j;
        point[2] = k;
        vnl_matrix_fixed< double, dim, dof >  Jth = transform->jacobian_wrt_par( point );
        vnl_matrix_fixed< double, dim, dim >  covarianceJCurrent = Jth * covariance * Jth.transpose();
        
        covarianceJ += covarianceJCurrent;
        
        double currTrace = vnl_trace( covarianceJCurrent );
        if( currTrace > maxTrace ) maxTrace = currTrace;
        
        vnl_matrix< double >  V( 3, 3 );
        vnl_vector< double >  D( 3 );
        vnl_symmetric_eigensystem_compute( covarianceJCurrent, V, D );
        
        if( D( 2 ) > maxEval ) maxEval = D( 2 );
        /* vcl_cout << "D:;t" << (point-center.GetVnlVector()).two_norm() << ";t" << D( 0 ) << ";t" << D( 1 ) << ";t" << D( 2 ) << vcl_endl; */   \
        
        ++count;
      }
    }
  }

  covarianceJ *= ( 1.0 / count );


//#define COMPUTE_TRANSFER_ERROR_COV( i, j, k ) \
//    vnl_vector_fixed< double, 3 >  point;   \
//    point[0] = i;   \
//    point[1] = j;   \
//    point[2] = k;   \
//    vnl_matrix_fixed< double, dim, dof >  Jth = transform->jacobian_wrt_par( point );   \
//    vnl_matrix_fixed< double, dim, dim >  covarianceJCurrent = Jth * covariance * Jth.transpose();   \
//    \
//    covarianceJ += covarianceJCurrent;   \
//    \
//    double currTrace = vnl_trace( covarianceJCurrent );   \
//    if( currTrace > maxTrace ) maxTrace = currTrace;   \
//    \
//    vnl_matrix< double >  V( 3, 3 );   \
//    vnl_vector< double >  D( 3 );   \
//    vnl_symmetric_eigensystem_compute( covarianceJCurrent, V, D );   \
//    \
//    if( D( 2 ) > maxEval ) maxEval = D( 2 );   \
//    /* vcl_cout << "D:;   \t" << (point-center.GetVnlVector()).two_norm() << ";   \t" << D( 0 ) << ";   \t" << D( 1 ) << ";   \t" << D( 2 ) << vcl_endl; */   \
//    \
//    ++count;
//
//
//  // the region is sampled only at the boundaries for efficiency reasons
//  for( double k = start[2]; k <= end[2]; k += (size[2]-spacing[2]) ) {
//    for( double j = start[1]; j <= end[1]; j += interval ) {
//      for( double i = start[0]; i <= end[0]; i += interval ) {
//        COMPUTE_TRANSFER_ERROR_COV( i, j, k );
//      }
//    }
//  }
//
//  // the region is sampled only at the boundaries for efficiency reasons
//  for( double k = start[2]; k <= end[2]; k += interval ) {
//    for( double j = start[1]; j <= end[1]; j += (size[1]-spacing[1]) ) {
//      for( double i = start[0]; i <= end[0]; i += interval ) {
//        COMPUTE_TRANSFER_ERROR_COV( i, j, k )
//      }
//    }
//  }
//
//  // the region is sampled only at the boundaries for efficiency reasons
//  for( double k = start[2]; k <= end[2]; k += interval ) {
//    for( double j = start[1]; j <= end[1]; j += interval ) {
//      for( double i = start[0]; i <= end[0]; i += (size[0]-spacing[0]) ) {
//        COMPUTE_TRANSFER_ERROR_COV( i, j, k )
//      }
//    }
//  }
//  covarianceJ *= ( 1.0 / count );

}


bool LocationRegistration::validTransform( vbl_smart_ptr< cdcl_trans< dim, dof > > const &  trans )
{
  vnl_svd< double >  svd( trans->get_A() );

  vnl_diag_matrix< vnl_svd< double >::singval_t >  W = svd.W();

  // the scaling given as singular values of affine part A must be within bounds
  for( unsigned int d = 0; d < dim; ++d ) {
    if( W( d, d ) > 10.0 || W( d, d ) < 0.1 ) return false;
  }

  return true;
}



}


#endif
