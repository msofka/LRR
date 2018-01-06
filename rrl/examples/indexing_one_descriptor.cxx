#include <vul/vul_file.h>
#include <vul/vul_arg.h>

#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include "itkBSplineDeformableTransform.h"
#include "itkTransformFileReader.h"
#include "itkAffineTransform.h"
#include "itkTransformFactory.h"
#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkMaximumProjectionImageFilter.h"
#include "itkResampleImageFilter.h"
#include "itkTileImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"

#include "itkImageSeriesReader.h"
#include "itkImageSeriesWriter.h"
#include "itkVTKImageIO.h"
#include "itkPNGImageIO.h"

#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"

#include "itkKdTree.h"
#include "itkKdTreeGenerator.h"
#include "itkListSample.h"
#include "itkHistogram.h"
#include "itkBoundingBox.h"

#include "itkTimeProbe.h"

#include <cdcl/cdcl_feature.h>
#include <cdcl/cdcl_utils_VTK.h>

#include <rsdl/rsdl_kd_tree.h>

#include <vnl/algo/vnl_symmetric_eigensystem.h>
#include <vnl/algo/vnl_generalized_eigensystem.h>

#include <vcl_utility.h>
#include <vcl_algorithm.h>

#include "vtkXMLPolyDataReader.h"
#include "cdcl/io/vtkPolyDataToFeaturesWithShapeFilter.h"


//:
// \file
// \brief  Invariant indexing with shape context in 3D, use ITK kd-tree for indexing.
//         There is some unfortunate naming here. When using synthetic BSpline transform, fixed image is the original
//         and moving image is the deformed image (it was assumed that the deformed image will be registered with the original).
//         Therefore indexing is done from the fixed to moving (don't have inverse for BSpline).
//         When using deformation field as a result of deformable registration, the transform is used to map moving image
//         to the fixed. So the indexing will be from moving to fixed. This would be a better naming had we switched the role
//         of fixed/moving for the BSpline case. Now when talking about warping field, fixed/moving are switched (in the code).
// \author Michal Sofka
// \date   January 2008


// write debugging image chips for each descriptor
#define OUTPUT_IMAGES 0


namespace {

// Kd-tree for descriptor vectors
//static const unsigned int descriptorSize = 771;
//typedef itk::Vector< float, descriptorSize >                      DescriptorVectorType;
typedef itk::Array< float >                                       DescriptorVectorType;
typedef itk::Statistics::ListSample< DescriptorVectorType >       DescriptorSampleType;
typedef itk::Statistics::KdTreeGenerator< DescriptorSampleType >  DescriptorTreeGeneratorType;
typedef DescriptorTreeGeneratorType::KdTreeType                   DescriptorTreeType;


// Kd-tree for keypoint locations
typedef itk::Vector< float, 3 >                                 KeypointVectorType;
typedef itk::Statistics::ListSample< KeypointVectorType >       KeypointSampleType;
typedef itk::Statistics::KdTreeGenerator< KeypointSampleType >  KeypointTreeGeneratorType;
typedef KeypointTreeGeneratorType::KdTreeType                   KeypointTreeType;


struct match {
  match() : moving_keypoint_( 0 ), fixed_keypoint_( 0 ), distance_( 0.0 ) {}
  match( cdcl_keypoint< 3 >::sptr const &  moving_keypoint,
         vnl_vector< float > const &  moving_descriptor,
         cdcl_keypoint< 3 >::sptr const &  fixed_keypoint,
         vnl_vector< float > const &  fixed_descriptor,
         double const &  distance )
    : moving_keypoint_( moving_keypoint ),
      moving_descriptor_( moving_descriptor ),
      fixed_keypoint_( fixed_keypoint ),
      fixed_descriptor_( fixed_descriptor ),
      distance_( distance ) {}
  cdcl_keypoint< 3 >::sptr  moving_keypoint_;
  vnl_vector< float >  moving_descriptor_;
  cdcl_keypoint< 3 >::sptr  fixed_keypoint_;
  vnl_vector< float >  fixed_descriptor_;
  double distance_;
};

bool DistanceAtIndexIsSmaller( match const &  left, match const &  right )
{
  return left.distance_ < right.distance_;
}



// Function object for sorting vector of keypoint indices based on distance
class LocationAtIndexIsGreater {
public:
  LocationAtIndexIsGreater( KeypointVectorType &  FixedKeypoint,
                            KeypointTreeType::Pointer const &  MovingKeypointTree )
    : m_FixedKeypoint( FixedKeypoint ),
      m_MovingKeypointTree( MovingKeypointTree ) {}
  bool operator()( int i, int j )
    {
      KeypointVectorType  movingKeypoint_i = m_MovingKeypointTree->GetMeasurementVector( i );
      vnl_vector< float >  distancei = m_FixedKeypoint.GetVnlVector() - movingKeypoint_i.GetVnlVector();
      KeypointVectorType  movingKeypoint_j = m_MovingKeypointTree->GetMeasurementVector( j );
      vnl_vector< float >  distancej = m_FixedKeypoint.GetVnlVector() - movingKeypoint_j.GetVnlVector();
      return distancei.magnitude() < distancej.magnitude();
    }
private:
  KeypointVectorType const &  m_FixedKeypoint;
  KeypointTreeType::Pointer const &  m_MovingKeypointTree;
};


//bool DistanceForMeasurementIsGreater( vcl_pair< DescriptorVectorType, double > const &  left,
//                                      vcl_pair< DescriptorVectorType, double > const &  right )
//{
//  return left.second < right.second;
//}

class DistanceID {
  public:
    DescriptorTreeType::InstanceIdentifier  id_;
    DescriptorVectorType  vector_;
    double distance_;
};

bool DistanceForMeasurementIsGreater( DistanceID const &  left,
                                      DistanceID const &  right )
{
  return left.distance_ < right.distance_;
}


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

  size[0] = 70 / spacing[0];
  size[1] = 70 / spacing[1];
  size[2] = 70 / spacing[2];

  index[0] = index[0] - size[0] / 2;
  index[1] = index[1] - size[1] / 2;
  index[2] = index[2] - size[2] / 2;

  m_ROI.SetSize( size );
  m_ROI.SetIndex( index );

  bool canCrop = m_ROI.Crop( Image->GetLargestPossibleRegion() );
  if( !canCrop || !Image->GetLargestPossibleRegion().IsInside( index ) ) {
    std::cout << "Error: point outside the image. Or the ROI outside of image." << std::endl;
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
  bool canCrop = extractionRegion.Crop( Image->GetLargestPossibleRegion() );
  if( !canCrop || !Image->GetLargestPossibleRegion().IsInside( index ) ) {
    std::cout << "Error: point outside the image. Or the ROI outside of image." << std::endl;
    return;
  }


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
  bool canCrop = extractionRegion.Crop( Image->GetLargestPossibleRegion() );
  if( !canCrop || !Image->GetLargestPossibleRegion().IsInside( index ) ) {
    std::cout << "Error: point outside the image. Or the ROI outside of image." << std::endl;
    return;
  }

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
  movingMappedIndex[1] = movingMappedIndex[1] * resamplingFactor;
  movingMappedResampler->SetOutputStartIndex( movingMappedIndex);

  typename SliceImageType::SizeType     movingMappedSize    = extractedYZ->GetLargestPossibleRegion().GetSize();
  movingMappedSize[1] = movingMappedSize[1] * resamplingFactor;
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
  bool canCrop = extractionRegion.Crop( Image->GetLargestPossibleRegion() );
  if( !canCrop || !Image->GetLargestPossibleRegion().IsInside( index ) ) {
    std::cout << "Error: point outside the image. Or the ROI outside of image." << std::endl;
    return;
  }

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
  movingMappedIndex[1] = movingMappedIndex[1] * resamplingFactor;
  movingMappedResampler->SetOutputStartIndex( movingMappedIndex);

  typename SliceImageType::SizeType     movingMappedSize    = extractedXZ->GetLargestPossibleRegion().GetSize();
  movingMappedSize[1] = movingMappedSize[1] * resamplingFactor;
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
  writer->SetImageIO( io );
  writer->SetInput( rescaler->GetOutput() );

  writer->SetFileName( FileName.c_str() );
  writer->Update();

}


template <class ImageType>
void WriteROIImage( typename ImageType::Pointer const &  Image,
                    cdcl_keypoint< 3 >::sptr const &  Keypoint,
                    vcl_string  FileName )
{
  itk::Point< double, 3 >  Point;
  Point[0] = Keypoint->location_[0];
  Point[1] = Keypoint->location_[1];
  Point[2] = Keypoint->location_[2];

  WriteROIImage<ImageType>( Image, Point, FileName );
}

}


int
main( int argc, char* argv[] )
{
  // deformation field case:
  // suppose we have a deformation field which warps image1 to image2
  // then warper in itk warps image image1 (output) to image2 input, s.t. input = output + deformation
  // therefore movingDescriptorFile is the one of image2 and fixedDescriptorFile is the one of image1
  //
  // e.g. indexing_one_descriptor.exe 241/4804/000002croppeddesc0.vtk 241/4802/000002croppeddesc0.vtk -trans 241/4802_000002cropped-4804_000002cropped.vtk -mov 241/4804/000002cropped -fix 241/4802/000002cropped
  // e.g. indexing_one_descriptor.exe 2505/4742/000002croppeddesc1.vtk 2505/4741/000002croppeddesc1.vtk -trans 2505/4741_000002cropped-4742_000002cropped.vtk -mov 2505/4742/000002cropped -fix 2505/4741/000002cropped
  vul_arg< const char* >        movingDescriptorFile     ( 0, "VTK moving descriptor file (of a synthetically deformed image)" );
  vul_arg< const char* >        fixedDescriptorFile      ( 0, "VTK fixed descriptor file" );
  vul_arg< const char* >        transformFile            ( "-trans", "ITK transform file (transform used to warp points)", 0 );
  vul_arg< const char* >        movingFeatureFilePrefix  ( "-fea", "VTK moving feature file prefix (of a deformed image) to test keypoint repeatability" );
  vul_arg< unsigned int >       levelInd                 ( "-ind", "level index of the feature file (-fea)" );
  vul_arg< const char* >        movingImageFile          ( "-mov", "moving image file (of a deformed image)" );
  vul_arg< const char* >        fixedImageFile           ( "-fix", "fixed image file" );
  vul_arg< const char* >        locationsFile            ( "-locs", "File with a list of initial locations (X Y Z on each line)", 0 );

  vul_arg_parse( argc, argv );

  vcl_cout << "Command: ";
  for( int i = 0; i < argc; ++i )
    vcl_cout << argv[i] << " ";
  vcl_cout << vcl_endl;


  // Read image volumes
  //
  bool volumesRead = false;
  typedef float PixelType;
  const     unsigned int   Dimension = 3;
  typedef itk::Image< PixelType, Dimension >  FixedImageType;
  typedef itk::Image< PixelType, Dimension >  MovingImageType;
  FixedImageType::Pointer  fixedImage;
  MovingImageType::Pointer  movingImage;
  if( movingImageFile.set() && fixedImageFile.set() ) {


    // Try using image reader first (for mhd file type, for example).
    //
    typedef itk::ImageFileReader< FixedImageType > FixedImageReaderType;
    FixedImageReaderType::Pointer fixedImageReader = FixedImageReaderType::New();
    fixedImageReader->SetFileName( fixedImageFile() );

    fixedImageReader->Update();

    fixedImage = fixedImageReader->GetOutput();


    typedef itk::ImageFileReader< MovingImageType > MovingImageReaderType;
    MovingImageReaderType::Pointer movingImageReader = MovingImageReaderType::New();
    movingImageReader->SetFileName( movingImageFile() );

    movingImageReader->Update();

    movingImage = movingImageReader->GetOutput();

    volumesRead = fixedImage.IsNotNull() && movingImage.IsNotNull();


    if( !volumesRead ) {
      // Try reading Dicom

      const unsigned int Dimension = 3;
      //typedef signed short PixelType;
      typedef float PixelType;

      typedef   itk::ImageSeriesReader< FixedImageType >   FixedReaderType;
      typedef   itk::ImageSeriesReader< MovingImageType >  MovingReaderType;


      // Read the input image.
      // USE GDCMImageIO, DICOMImageIO2 is OLD
      typedef itk::GDCMImageIO                        ImageIOType;
      typedef itk::GDCMSeriesFileNames                NamesGeneratorType;
      
      ImageIOType::Pointer io = ImageIOType::New();

      // Get the DICOM filenames from the directory
      NamesGeneratorType::Pointer  fixedNames = NamesGeneratorType::New();
      fixedNames->SetDirectory( fixedImageFile() );
      const FixedReaderType::FileNamesContainer &  fixedFileNames = 
                                fixedNames->GetInputFileNames();

      FixedReaderType::Pointer  fixedImageReader = FixedReaderType::New();
      fixedImageReader->SetFileNames( fixedFileNames );
      fixedImageReader->SetImageIO( io );
      // WRITER DOES NOT HAVE ReverseOrderOn FUNCTION! -> It would get flipped
      //fixedReader->ReverseOrderOn();
      std::cout << fixedNames;


      // Get the DICOM filenames from the directory
      NamesGeneratorType::Pointer  movingNames = NamesGeneratorType::New();
      movingNames->SetDirectory( movingImageFile() );
      const FixedReaderType::FileNamesContainer &  movingFileNames = 
                                movingNames->GetInputFileNames();

      MovingReaderType::Pointer  movingImageReader = MovingReaderType::New();
      movingImageReader->SetFileNames( movingFileNames );
      movingImageReader->SetImageIO( io );
      std::cout << movingNames;

      fixedImage = fixedImageReader->GetOutput();
      movingImage = movingImageReader->GetOutput();

      fixedImageReader->Update();
      movingImageReader->Update();
      std::cout << "Fixed spacing: " << fixedImage->GetSpacing() << std::endl;
      std::cout << "Moving spacing: " << movingImage->GetSpacing() << std::endl;

      volumesRead = fixedImage.IsNotNull() && movingImage.IsNotNull();
    }
  }


#define USE_ALSO_ITK_KDTREE 0


  // setup indexing
#if USE_ALSO_ITK_KDTREE
  DescriptorSampleType::Pointer  sampleMovingDescriptors = DescriptorSampleType::New();
#endif
  KeypointSampleType::Pointer    sampleMovingKeypoints = KeypointSampleType::New();

  vcl_vector< cdcl_keypoint< 3 >::sptr >                movingKeypoints;
  // read moving descriptors
// can be dropped below if rsdl kd-tree is taken out, also comment in the line below: empty_vector.swap( movingDescriptors );
vcl_vector< vnl_vector< float > >                    movingDescriptors;
  {

    std::cout << "Reading moving keypoint descriptors..." << std::endl;
    cdcl_read_keypoint_descriptors_VTK( movingKeypoints, movingDescriptors, movingDescriptorFile() );

#if USE_ALSO_ITK_KDTREE
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
#endif    
    // recommended hack on releasing memory is to swap in empty vector
    vcl_vector< vnl_vector< float > >  empty_vector;
//empty_vector.swap( movingDescriptors );


    if( movingFeatureFilePrefix.set() && levelInd.set() ) {
      vcl_cout << "Using features for keypoint testing." << vcl_endl;

      // read features
      vcl_vector< vcl_vector< vbl_smart_ptr< cdcl_feature_with_shape< 3 > > > >  multires_points;
      vtkSmartPointer< vtkXMLPolyDataReader >  featureReader = vtkSmartPointer< vtkXMLPolyDataReader >::New();
      vcl_string  filename = vcl_string( movingFeatureFilePrefix() ) + "_00.vtk";
      featureReader->SetFileName( filename.c_str() );
      featureReader->Update();  // need this update (bug in the converting filter)

      // create the converting filter
      vtkSmartPointer< vtkPolyDataToFeaturesWithShapeFilter >  polyDataToFeaturesFilter = vtkSmartPointer< vtkPolyDataToFeaturesWithShapeFilter >::New();
      polyDataToFeaturesFilter->SetInput( featureReader->GetOutput() );
      polyDataToFeaturesFilter->Update();
      vtkSmartPointer< vtkFeatureWithShapeAttributeSet >  featureAttributeSetVTK = polyDataToFeaturesFilter->GetOutput();

      typedef vtkFeatureWithShapeAttributeSet::FeatureSetType  FeatureSetType;
      FeatureSetType &  features = featureAttributeSetVTK->GetPoints();
      multires_points.push_back( features );


      KeypointVectorType  mv_k;
      for( vcl_vector< vbl_smart_ptr< cdcl_feature< 3 > > >::size_type  n = 0; n < multires_points[levelInd()].size(); ++n ) {
        mv_k[0] = multires_points[levelInd()][n]->location_[0];
        mv_k[1] = multires_points[levelInd()][n]->location_[1];
        mv_k[2] = multires_points[levelInd()][n]->location_[2];
        sampleMovingKeypoints->PushBack( mv_k );
      }
    }
    else {
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

  }

#if USE_ALSO_ITK_KDTREE
  // moving keypoint descriptor vectors
  DescriptorTreeGeneratorType::Pointer  descriptorTreeGenerator = DescriptorTreeGeneratorType::New();

  descriptorTreeGenerator->SetSample( sampleMovingDescriptors );
  descriptorTreeGenerator->SetBucketSize( 16 );
  descriptorTreeGenerator->Update();

  DescriptorTreeType::Pointer  kdTreeMovingDescriptors = descriptorTreeGenerator->GetOutput();

  vcl_cout << "Kd tree size: " << kdTreeMovingDescriptors->Size() << vcl_endl;
#endif

  // The same moving descriptor tree but with rsdl_kdtree
  //
  // Obtain the number of cartesian, nc, and number of angular, na,
  // from the first invariant feature
  unsigned int nc = movingDescriptors[0].size();
  unsigned int na = 0;

  // Construct a kd-tree for the set of search points from the fixed_set
  vcl_vector< rsdl_point > rsdlTreeMovingDescriptors( movingDescriptors.size() );
  for( unsigned int pt = 0; pt < movingDescriptors.size(); ++pt ) {
    // Set number of expected cartesian and angular values
    rsdlTreeMovingDescriptors[pt].resize( nc, na );
    // unfortunatelly, rsdl does not work with floats, so we need to convert
    vnl_vector< double >  movingDescriptors_double( movingDescriptors[pt].size() );
    for( unsigned int d = 0; d < movingDescriptors[pt].size(); ++d ) {
      movingDescriptors_double[d] = movingDescriptors[pt][d];
    }
    rsdlTreeMovingDescriptors[pt].set_cartesian( movingDescriptors_double );
  }
  // use default parameters: each node has 4 data points
  //
  rsdl_kd_tree kd_tree( rsdlTreeMovingDescriptors, 0, 4 );

  double param_max_leaves_ratio_ = 0.002;
  unsigned int param_leaves_lower_bound_ = 15;
  unsigned int param_leaves_upper_bound_ = 200;

  unsigned int max_leaves = (unsigned int)( param_max_leaves_ratio_ * (double)rsdlTreeMovingDescriptors.size());

  // bound the max_leaves between 15 and 200 (based off the yang-smith heuristic)
  // recall each leaf node in Kd tree has 4 data points by default.
  //
  max_leaves = max_leaves < param_leaves_lower_bound_ ? param_leaves_lower_bound_ : max_leaves;
  max_leaves = max_leaves > param_leaves_upper_bound_ ? param_leaves_upper_bound_ : max_leaves;

  //max_leaves = 200;
  std::cout << "max_leaves: " << max_leaves << std::endl;

  // moving keypoint locations
  KeypointTreeGeneratorType::Pointer  keypointTreeGenerator = KeypointTreeGeneratorType::New();

  keypointTreeGenerator->SetSample( sampleMovingKeypoints );
  keypointTreeGenerator->SetBucketSize( 16 );
  keypointTreeGenerator->Update();

  KeypointTreeType::Pointer  kdTreeMovingKeypoints = keypointTreeGenerator->GetOutput();


  // read fixed descriptors
  //
  vcl_vector< cdcl_keypoint< 3 >::sptr >                fixedKeypoints;
  vcl_vector< vnl_vector< float > >                    fixedDescriptors;
  std::cout << "Reading fixed keypoint descriptors..." << std::endl;
  cdcl_read_keypoint_descriptors_VTK( fixedKeypoints, fixedDescriptors, fixedDescriptorFile() );

  std::cout << "Have read: " << fixedKeypoints.size() << " keypoints " << std::endl;

  KeypointSampleType::Pointer    sampleFixedKeypoints = KeypointSampleType::New();

  // read fixed keypoints
  {
    // fixed keypoint locations
    //
    KeypointVectorType  fi_k;
    for( unsigned int n = 0 ; n < fixedKeypoints.size(); ++n ) {
      fi_k[0] = fixedKeypoints[n]->location_[0];
      fi_k[1] = fixedKeypoints[n]->location_[1];
      fi_k[2] = fixedKeypoints[n]->location_[2];
      //mv.Set_vnl_vector( m_FixedKeypoints[n]->location_.as_vector() );
      sampleFixedKeypoints->PushBack( fi_k );
      //std::cout << fixedKeypoints[n]->location_ << std::endl;
    }

  }

  // fixed keypoint locations
  KeypointTreeGeneratorType::Pointer  fixedKeypointTreeGenerator = KeypointTreeGeneratorType::New();

  fixedKeypointTreeGenerator->SetSample( sampleFixedKeypoints );
  fixedKeypointTreeGenerator->SetBucketSize( 16 );
  fixedKeypointTreeGenerator->Update();

  KeypointTreeType::Pointer  kdTreeFixedKeypoints = fixedKeypointTreeGenerator->GetOutput();



  // setup transforms

  typedef itk::AffineTransform<double,3> AffineTransformType;

  const unsigned int ImageDimension = 3;
	const unsigned int SpaceDimension = ImageDimension;
	const unsigned int SplineOrder = 3;
	typedef double CoordinateRepType;

	typedef itk::BSplineDeformableTransform<
													CoordinateRepType,
													SpaceDimension,
													SplineOrder >     BSplineTransformType;

  BSplineTransformType::Pointer  bsplineTransform;

  // Deformation field
  typedef itk::Vector< float, ImageDimension >           VectorPixelType;
  typedef itk::Image<  VectorPixelType, ImageDimension > DeformationFieldType;
  DeformationFieldType::Pointer  deformationField;

  enum deformationTypeEnum { NONE, BSPLINE, FIELD };
  deformationTypeEnum  deformationType = NONE;

  if( transformFile.set() ) {
    vcl_string  transformExtension = itksys::SystemTools::GetFilenameLastExtension( transformFile() );
    //std::cout << "Transform extension: " << transformExtension << std::endl;
    // if warp field supplied
    if( transformExtension == ".vtk" || transformExtension == ".mhd" ) {
      typedef   itk::ImageFileReader< DeformationFieldType >  FieldReaderType;

      std::cout << "Reading deformation field: " << transformFile();

      typedef itk::VTKImageIO                        ImageIOType;
      ImageIOType::Pointer io = ImageIOType::New();

      if( itksys::SystemTools::FileExists( transformFile() ) && io->CanReadFile( transformFile() ) ) {
        deformationType = FIELD;

        FieldReaderType::Pointer fieldReader = FieldReaderType::New();
        fieldReader->SetFileName( transformFile() );
        fieldReader->SetImageIO( io );
        fieldReader->Update();
        deformationField = fieldReader->GetOutput();

        std::cout << " of size " << deformationField->GetLargestPossibleRegion().GetSize() << std::endl;
      }
      else {
        std::cout << "Warning: The field " << transformFile() << " cannot be read." << std::endl;
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
      reader->SetFileName( transformFile() );
    
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
    deformationType = BSPLINE;
  }


#if OUTPUT_IMAGES
  typedef itk::Image<float, 2>  ImageType;
  //  Get information about origin, spacing and region of the images
  ImageType::PointType    imageOrigin;
  imageOrigin.Fill( 0.0 );
  ImageType::SpacingType  imageSpacing;
  imageSpacing.Fill( 1.0 );
  ImageType::RegionType   imageRegion;
  ImageType::RegionType::SizeType  imageSize;
  //imageSize[0] = 3*8*8;//(unsigned int) ( vcl_ceil( std::sqrt( float( fixedDescriptors[0].size() ) ) ) );   // when the directions are not flipped
  //imageSize[0] = 3*4*4;//(unsigned int) ( vcl_ceil( std::sqrt( float( fixedDescriptors[0].size() ) ) ) );   // gradient mirroring
  imageSize[0] = 2*2*2*6;//(unsigned int) ( vcl_ceil( std::sqrt( float( fixedDescriptors[0].size() ) ) ) );   // covariance context
  imageSize[1] = fixedDescriptors[0].size() / imageSize[0];
  imageRegion.SetSize( imageSize );
  ImageType::RegionType::IndexType imageIndex;
  imageIndex.Fill( 0 );
  imageRegion.SetIndex( imageIndex );

  //  Set up and allocate the response image
  ImageType::Pointer  descIm = ImageType::New();
  descIm->SetOrigin( imageOrigin );
  descIm->SetSpacing( imageSpacing );
  descIm->SetRegions( imageRegion );
  descIm->Allocate();
  descIm->FillBuffer( 0.0f );

  typedef itk::Image<vxl_byte, 2>                                            WriteImageType;
  typedef itk::ImageFileWriter<WriteImageType>                               WriterType;

  typedef itk::RescaleIntensityImageFilter<ImageType, WriteImageType>  RescalerType;
  //  Rescaler to set the final image to 0..255
  RescalerType::Pointer rescaler = RescalerType::New();
  rescaler->SetOutputMinimum( 0 );
  rescaler->SetOutputMaximum( 255 );
  rescaler->SetInput( descIm );

  //  Writer
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( rescaler->GetOutput() );
  

  typedef itk::ImageRegionIterator< ImageType >           IteratorType;
  IteratorType  descItr( descIm, descIm->GetLargestPossibleRegion() );
#endif


  // histograms of failures
  const unsigned int num_angles = 19;
  const unsigned int num_distances = 4;
  unsigned int normals_distances[num_angles][num_distances];
  unsigned int binormals_distances[num_angles][num_distances];
  for( unsigned int a = 0; a < num_angles; ++a )
    for( unsigned int d = 0; d < num_distances; ++d ) {
      normals_distances[a][d] = 0;
      binormals_distances[a][d] = 0;
    }

  unsigned int keypoints_normals_distances[num_angles][num_distances];
  unsigned int keypoints_binormals_distances[num_angles][num_distances];
  for( unsigned int a = 0; a < num_angles; ++a )
    for( unsigned int d = 0; d < num_distances; ++d ) {
      keypoints_normals_distances[a][d] = 0;
      keypoints_binormals_distances[a][d] = 0;
    }




  typedef  itk::BoundingBox< unsigned long, 3, double >  BoundingBoxType;
  typedef BoundingBoxType::PointType  PointType;

  BoundingBoxType::Pointer  fixedBoundingBox = BoundingBoxType::New();
  BoundingBoxType::PointsContainerPointer  fixedPoints = BoundingBoxType::PointsContainer::New();
  
  for( vcl_vector< cdcl_keypoint< 3 >::sptr >::size_type  i = 0; i < fixedKeypoints.size(); ++i ) {
    PointType  point;
    point[0] = fixedKeypoints[i]->location_[0];
    point[1] = fixedKeypoints[i]->location_[1];
    point[2] = fixedKeypoints[i]->location_[2];
    fixedPoints->push_back( point );
  }
  fixedBoundingBox->SetPoints( fixedPoints );

  PointType  minFixedPoint = fixedBoundingBox->GetMinimum();
  PointType  maxFixedPoint = fixedBoundingBox->GetMaximum();



  BoundingBoxType::Pointer  movingBoundingBox = BoundingBoxType::New();
  BoundingBoxType::PointsContainerPointer  movingPoints = BoundingBoxType::PointsContainer::New();
  
  for( vcl_vector< cdcl_keypoint< 3 >::sptr >::size_type  i = 0; i < movingKeypoints.size(); ++i ) {
    PointType  point;
    point[0] = movingKeypoints[i]->location_[0];
    point[1] = movingKeypoints[i]->location_[1];
    point[2] = movingKeypoints[i]->location_[2];
    movingPoints->push_back( point );
  }
  movingBoundingBox->SetPoints( movingPoints );

  PointType  minMovingPoint = movingBoundingBox->GetMinimum();
  PointType  maxMovingPoint = movingBoundingBox->GetMaximum();

  // transform current keypoint and find out whether it is a good one or a bad one
  //
  // moving image was obtained by deforming the fixed image, so deform the fixed keypoint with the same transform
  // and compare the two points
  //BSplineTransformType::OutputPointType  movingMapped = bsplineTransform->TransformPoint( moving_keypoint->location_.data_block() );
  vnl_vector_fixed< double, 3 >  minFixedPointMappedVnl;
  vnl_vector_fixed< double, 3 >  maxFixedPointMappedVnl;
  if( deformationType != NONE ) {

    // mapping error
    vnl_vector_fixed< double, 3 >  distance;
    if( deformationType == BSPLINE ) {
      // moving image was obtained by deforming the fixed image, so deform the fixed keypoint with the same transform
      // and compare the two points
      BSplineTransformType::OutputPointType  minFixedPointMapped = bsplineTransform->TransformPoint( minFixedPoint );
      minFixedPointMappedVnl = minFixedPointMapped.GetVnlVector();
      BSplineTransformType::OutputPointType  maxFixedPointMapped = bsplineTransform->TransformPoint( maxFixedPoint );
      maxFixedPointMappedVnl = maxFixedPointMapped.GetVnlVector();

    }
    else if( deformationType == FIELD ) {

      // get the deformation vector from the field
      DeformationFieldType::IndexType  minFixedIndex;
      deformationField->TransformPhysicalPointToIndex( minFixedPoint, minFixedIndex );
      VectorPixelType  deformationVector = deformationField->GetPixel( minFixedIndex );

      itk::Point< double, 3 >  minFixedPointMapped = minFixedPoint + deformationVector;
      minFixedPointMappedVnl = minFixedPointMapped.GetVnlVector();

      // get the deformation vector from the field
      DeformationFieldType::IndexType  maxFixedIndex;
      deformationField->TransformPhysicalPointToIndex( maxFixedPoint, maxFixedIndex );
      deformationVector = deformationField->GetPixel( maxFixedIndex );

      itk::Point< double, 3 >  maxFixedPointMapped = maxFixedPoint + deformationVector;
      maxFixedPointMappedVnl = maxFixedPointMapped.GetVnlVector();
    }
  }

  std::cout << "Bounding box, moving: " << minMovingPoint << ", " << maxMovingPoint 
                          << " fixed: " << minFixedPoint << ", " << maxFixedPoint
                          << " fixed mapped: " << minFixedPointMappedVnl << ", " << maxFixedPointMappedVnl << std::endl;



  typedef itk::Vector<PointType::CoordRepType, PointType::PointDimension>  VectorType;
  VectorType  maxTranslation = maxFixedPoint - minMovingPoint;
  VectorType  minTranslation = minFixedPoint - maxMovingPoint;


  std::string  prefix;  // prefix with the generated matches
  std::ifstream  noduleFile;
  if( locationsFile.set() ) {
    noduleFile.open( locationsFile() );
    if( !noduleFile ) {
      std::cout << "Specified nodule file cannot be opened: " << locationsFile() << std::endl;
      return 1;
    }
    prefix = itksys::SystemTools::GetFilenameWithoutExtension( locationsFile() ); //"nodules";
  }
  else prefix = "matches";


  // setup directory for descriptor matches
  std::string  fixedDesc = itksys::SystemTools::GetFilenameWithoutLastExtension( fixedDescriptorFile() );
  std::string  movingDesc = itksys::SystemTools::GetFilenameWithoutLastExtension( movingDescriptorFile() );
  std::string  fixedPatient = itksys::SystemTools::GetFilenameName( itksys::SystemTools::GetParentDirectory( fixedDescriptorFile() ) );
  std::string  movingPatient = itksys::SystemTools::GetFilenameName( itksys::SystemTools::GetParentDirectory( movingDescriptorFile() ) );
  std::string  matchesDir = prefix + fixedPatient + "_" + fixedDesc + "-" + movingPatient + "_" + movingDesc;
  std::string  patientDir = itksys::SystemTools::GetParentDirectory( itksys::SystemTools::GetParentDirectory( fixedDescriptorFile() ).c_str() );
  std::cout << patientDir << std::endl << matchesDir << std::endl;
  if( patientDir != "" ) patientDir += "/";
  std::string matchesDirPath = patientDir + matchesDir;
  itksys::SystemTools::MakeDirectory( matchesDirPath.c_str() );

  unsigned int numberOfKeypointsTested = 100;

  // total matches tested = numberOfKeypointsTested * numberOfMatchesTest, but numberOfMatchesTest is not always 20 (they might not be available if the neighborhood is small)
  unsigned int totalMatchesTested = 0;

  srand( (unsigned)time( NULL ) );


  std::cout << "Keypoint Distance\tDescriptor Difference\t" << std::endl;
  unsigned int goodAtFirst = 0;
  unsigned int goodIn10 = 0;
  unsigned int goodAtAll = 0;
  unsigned int goodCandidateKeypoints = 0;
  unsigned int badOrientationsAtFirst = 0;
  unsigned int badOrientationsAtAll = 0;

  float goodnessThreshold = 9.0f;
  float keypointGoodnessThreshold = 9.0f;

  // histogram, where each bin with index N is the number of times the good initialization was found N distance away
  std::vector< unsigned int >  distanceToFirstGood( 30, 0 );


//numberOfKeypointsTested = 12;
//vcl_vector< KeypointVectorType >  queryPoints( numberOfKeypointsTested );
//
//KeypointVectorType::ValueType  valueArr0[3] = {-52.619,	-104.727,	-294.702};
//queryPoints[0] = valueArr0;
//
//KeypointVectorType::ValueType  valueArr1[3] = {-76.3974,	-122.002,	-291.855};
//queryPoints[1] = valueArr1;
//
//KeypointVectorType::ValueType  valueArr2[3] = {-107.773,	-157.609,	-281.637};
//queryPoints[2] = valueArr2;
//
//KeypointVectorType::ValueType  valueArr3[3] = {-86.1652,	-107.833,	-287.36};
//queryPoints[3] = valueArr3;
//
//KeypointVectorType::ValueType  valueArr4[3] = {-85.776,	-161.692,	-198.453};
//queryPoints[4] = valueArr4;
//
//KeypointVectorType::ValueType  valueArr5[3] = {-86.8065,	-149.925,	-283.16};
//queryPoints[5] = valueArr5;
//
//KeypointVectorType::ValueType  valueArr6[3] = {-166.281,	-104.558,	-281.872};
//queryPoints[6] = valueArr6;
//
//KeypointVectorType::ValueType  valueArr7[3] = {-169.669,	-163.349,	-188.519};
//queryPoints[7] = valueArr7;
//
//KeypointVectorType::ValueType  valueArr8[3] = {-85.4691,	-109.791,	-250.002};
//queryPoints[8] = valueArr8;
//
//KeypointVectorType::ValueType  valueArr9[3] = {-53.2987,	-99.2357,	-279.551};
//queryPoints[9] = valueArr9;
//
//KeypointVectorType::ValueType  valueArr10[3] = {-64.1848,	-115.782,	-276.457};
//queryPoints[10] = valueArr10;
//
//KeypointVectorType::ValueType  valueArr11[3] = {-42.747,	-168.006,	-253.842};
//queryPoints[11] = valueArr11;
//

  std::string  matchesFile = fixedPatient + "_" + fixedDesc + "-" + movingPatient + "_" + movingDesc + "good.txt";
  std::string matchesFilePath = patientDir + matchesFile;

  std::ofstream  goodMatchesFile( matchesFilePath.c_str() );


  itk::TimeProbe  timer;
  itk::TimeProbe  timerITK;

  itk::TimeProbe  timerOneClick;

  // ****************
  // Keypoint testing
  // ****************

  for( unsigned int t = 0; t < numberOfKeypointsTested; ++t ) {
    timerOneClick.Start();
//if( t != 0 && t != 2 && t != 5 && t != 7 && t != 11 && t != 12 ) continue;
//KeypointVectorType  queryPoint = queryPoints[t];

    KeypointVectorType  queryPoint;
    if( locationsFile.set() ) { // load a location from a file
      if( !( noduleFile >> queryPoint[0] >> queryPoint[1] >> queryPoint[2] ) ) { // if we can't read anymore, finish
        t = numberOfKeypointsTested;
        continue;
      }
      else {
        std::cout << "Read nodule location: " << queryPoint << std::endl;
      }
    }
    else { // generate a random location
      for( unsigned int d = 0;  d < 3; ++d ) {
        // generate a number in a specific range
        // subtract radius of the shape context
        double context_radius = 30;
        double RANGE_MIN = minFixedPoint[d] + context_radius;
        double RANGE_MAX = maxFixedPoint[d] - context_radius;

        double randNum = ((double) rand() / (double) RAND_MAX) * (RANGE_MAX - RANGE_MIN) + RANGE_MIN;

        queryPoint[d] = randNum;
      }
    }
    std::cout << "Clicked: " << queryPoint << " ";

    double radius = 30.0;
    KeypointTreeType::InstanceIdentifierVectorType  fixedKeypointNeighbors;
    kdTreeFixedKeypoints->Search( queryPoint, radius, fixedKeypointNeighbors );

    if( fixedKeypointNeighbors.size() == 0 ) {
      std::cout << "No keypoints were found in the vicinity of the click, radius: " << radius << std::endl;
      continue;
    }

    // sort the keypoints based on the distance from the given location
    LocationAtIndexIsGreater  keypointCompare( queryPoint, kdTreeFixedKeypoints );
    // sort keypoint neighbors based on the distance from the given point
    std::sort( fixedKeypointNeighbors.begin(), fixedKeypointNeighbors.end(), keypointCompare );

    // Go through all fixed keypoints, find the closest descriptor match, compute candidate transformation (translation)
    // and enter it into the histogram. Record index of this histogram entry.
    //
    // indices to the histogram for each keypoint
    unsigned int numberOfKeypointNeighbors = fixedKeypointNeighbors.size();
    unsigned int goodKeypoints = 0;
    double  distanceToFirstGoodKeypoint = 0.0;

    vcl_vector< match >  matches;
    for( unsigned int l = 0; l < numberOfKeypointNeighbors; ++l ) {
      KeypointVectorType  fixedKeypoint = kdTreeFixedKeypoints->GetMeasurementVector( fixedKeypointNeighbors[l] );
      cdcl_keypoint< 3 >::sptr  fixed_keypoint  = fixedKeypoints[fixedKeypointNeighbors[l]];
      //std::cout << "fk: " << fixed_keypoint << std::endl;
     
      vnl_vector< float >  fixed_descriptor = fixedDescriptors[fixedKeypointNeighbors[l]];
      DescriptorVectorType  fixedDescriptor( fixed_descriptor.size() );
      //fixedDescriptor.SetData( fixed_descriptor.data_block() );
        for( unsigned int d = 0; d < fixed_descriptor.size(); ++d ) {
          //mv.SetSize( movingDescriptors[n].size() );
          //mv.SetData( movingDescriptors[n].data_block() );
          fixedDescriptor.SetElement( d, fixed_descriptor[d] );
        }

      // find descriptor in the moving image closest to the descriptor in the fixed image
      unsigned int numberOfDescriptorNeighbors = 1;
      DescriptorTreeType::InstanceIdentifierVectorType  descriptorNeighbors;

#if USE_ALSO_ITK_KDTREE
      timerITK.Start();
      kdTreeMovingDescriptors->Search( fixedDescriptor, numberOfDescriptorNeighbors, descriptorNeighbors ) ;
      timerITK.Stop();
      //std::cout << "ITK Tree query time: " << timerITK.GetMeanTime() << " seconds.\n";
#endif


#define SORT_ON_RATIO 0

      // rsdlTreeMovingDescriptors query 
      // Create a query point from the invariants of the current constellation
      rsdl_point query_pt( nc, na );
      // unfortunatelly, rsdl does not work with floats, so we need to convert
      vnl_vector< double >  moving_set_pt_double( fixed_descriptor.size() );
      for( unsigned int d = 0; d < fixed_descriptor.size(); ++d ) {
        moving_set_pt_double[d] = fixed_descriptor[d];
      }
      query_pt.set_cartesian( moving_set_pt_double );

      // Find the 2 nearest points
      vcl_vector< rsdl_point >  near_neighbor_pts;
      vcl_vector< int >  near_neighbor_indices;

      //kd_tree.n_nearest( query_pt, 2, near_neighbor_pts, near_neighbor_indices);
      timer.Start();
#if SORT_ON_RATIO
      kd_tree.n_nearest( query_pt, 2, near_neighbor_pts, near_neighbor_indices, true, max_leaves );
#else
      kd_tree.n_nearest( query_pt, 1, near_neighbor_pts, near_neighbor_indices, true, max_leaves );
#endif
      //kd_tree.n_nearest( query_pt, 1, near_neighbor_pts, near_neighbor_indices );
      timer.Stop();
      //std::cout << near_neighbor_indices[0] << "  " << descriptorNeighbors[0] << std::endl;
      descriptorNeighbors.push_back( near_neighbor_indices[0] );
      //std::cout << "rsdl Tree query time: " << timer.GetMeanTime() << " seconds.\n";



      // pick closest descriptor
      cdcl_keypoint< 3 >::sptr  moving_keypoint = movingKeypoints[descriptorNeighbors[0]];
#if USE_ALSO_ITK_KDTREE
      vnl_vector< float >  moving_descriptor = kdTreeMovingDescriptors->GetMeasurementVector( descriptorNeighbors[0] );
#else
      vnl_vector< float >  moving_descriptor = movingDescriptors[ descriptorNeighbors[0] ];
#endif

      vnl_vector_fixed< double, 3 >  fixedNormalMappedVnl;
      vnl_vector_fixed< double, 3 >  fixedBinormalMappedVnl;

      // transform current keypoint and find out whether it is a good one or a bad one
      //
      // moving image was obtained by deforming the fixed image, so deform the fixed keypoint with the same transform
      // and compare the two points
      //BSplineTransformType::OutputPointType  movingMapped = bsplineTransform->TransformPoint( moving_keypoint->location_.data_block() );
      vnl_vector_fixed< double, 3 >  fixedMappedVnl;
      if( deformationType != NONE ) {

        vnl_vector< double >  fixed_normal = fixed_keypoint->normal_;
        vnl_vector< double >  fixed_binormal = fixed_keypoint->binormal_;

        // mapping error
        vnl_vector_fixed< double, 3 >  distance;
        if( deformationType == BSPLINE ) {
          // moving image was obtained by deforming the fixed image, so deform the fixed keypoint with the same transform
          // and compare the two points
          BSplineTransformType::OutputPointType  fixedMapped = bsplineTransform->TransformPoint( fixed_keypoint->location_.data_block() );
          fixedMappedVnl = fixedMapped.GetVnlVector();

          // map normal and binormal
          BSplineTransformType::OutputPointType  fixedMappedNor = bsplineTransform->TransformPoint( ( fixed_keypoint->location_ + fixed_normal ).data_block() );
          BSplineTransformType::OutputPointType  fixedMappedBin = bsplineTransform->TransformPoint( ( fixed_keypoint->location_ + fixed_binormal ).data_block() );

          fixedNormalMappedVnl = ( fixedMappedNor.GetVnlVector() - fixedMappedVnl ).normalize();
          fixedBinormalMappedVnl = ( fixedMappedBin.GetVnlVector() - fixedMappedVnl ).normalize();
        }
        else if( deformationType == FIELD ) {

          // deformation field was obtained from registration 
          itk::Point< double, 3 >  fixedPoint;
          fixedPoint[0] = fixed_keypoint->location_[0];
          fixedPoint[1] = fixed_keypoint->location_[1];
          fixedPoint[2] = fixed_keypoint->location_[2];

          // get the deformation vector from the field
          DeformationFieldType::IndexType  fixedIndex;
          deformationField->TransformPhysicalPointToIndex( fixedPoint, fixedIndex );
          VectorPixelType  deformationVector = deformationField->GetPixel( fixedIndex );

          itk::Point< double, 3 >  fixedPointMapped = fixedPoint + deformationVector;
          fixedMappedVnl = fixedPointMapped.GetVnlVector();

          // map normal
          //
          // deformation field was obtained from registration 
          itk::Point< double, 3 >  fixedPointNor;
          fixedPointNor[0] = fixed_keypoint->location_[0] + fixed_normal[0];
          fixedPointNor[1] = fixed_keypoint->location_[1] + fixed_normal[1];
          fixedPointNor[2] = fixed_keypoint->location_[2] + fixed_normal[2];

          // get the deformation vector from the field
          DeformationFieldType::IndexType  fixedIndexNor;
          deformationField->TransformPhysicalPointToIndex( fixedPointNor, fixedIndexNor );
          VectorPixelType  deformationVectorNor = deformationField->GetPixel( fixedIndexNor );

          itk::Point< double, 3 >  fixedPointNorMapped = fixedPointNor + deformationVectorNor;
          fixedNormalMappedVnl = ( fixedPointNorMapped.GetVnlVector() - fixedMappedVnl ).normalize();


          // map binormal
          //
          // deformation field was obtained from registration 
          itk::Point< double, 3 >  fixedPointBin;
          fixedPointBin[0] = fixed_keypoint->location_[0] + fixed_binormal[0];
          fixedPointBin[1] = fixed_keypoint->location_[1] + fixed_binormal[1];
          fixedPointBin[2] = fixed_keypoint->location_[2] + fixed_binormal[2];

          // get the deformation vector from the field
          DeformationFieldType::IndexType  fixedIndexBin;
          deformationField->TransformPhysicalPointToIndex( fixedPointBin, fixedIndexBin );
          VectorPixelType  deformationVectorBin = deformationField->GetPixel( fixedIndexBin );

          itk::Point< double, 3 >  fixedPointBinMapped = fixedPointBin + deformationVectorBin;
          fixedBinormalMappedVnl = ( fixedPointBinMapped.GetVnlVector() - fixedMappedVnl ).normalize();

        }
      }


      // find moving keypoints which are closest to the mapped fixed keypoint
      KeypointVectorType  fixedKeypointMapped;
      for( unsigned int d = 0; d < 3; ++d ) {
        fixedKeypointMapped.SetElement( d, fixedMappedVnl[d] );
      }

      vnl_vector_fixed< double, 3 >  distance = fixedMappedVnl - moving_keypoint->location_;
      //vnl_vector< float >  difference = moving_descriptor - fixed_descriptor;
      //std::cout << " Match dist: " << distance.magnitude() << std::endl;
      //std::cout << difference.magnitude() << "\t";
      if( distance.magnitude() < goodnessThreshold ) {
        double normal_difference = vcl_acos( dot_product( moving_keypoint->normal_, fixedNormalMappedVnl ) ) * 180.0 / vnl_math::pi;
        double binormal_difference = vcl_acos( dot_product( moving_keypoint->binormal_, fixedBinormalMappedVnl ) ) * 180.0 / vnl_math::pi;
        if( normal_difference < 20.0
          && binormal_difference < 20.0 ) {
          if( goodKeypoints == 0 ) {
            vnl_vector_fixed< double, 3 >  vectorToFirstGoodKeypoint;
            vectorToFirstGoodKeypoint[0] = fixed_keypoint->location_[0] - queryPoint.GetVnlVector()[0];
            vectorToFirstGoodKeypoint[1] = fixed_keypoint->location_[1] - queryPoint.GetVnlVector()[1];
            vectorToFirstGoodKeypoint[2] = fixed_keypoint->location_[2] - queryPoint.GetVnlVector()[2];
            distanceToFirstGoodKeypoint = vectorToFirstGoodKeypoint.magnitude();
          }
          ++goodKeypoints;
        }
      }





      vnl_vector_fixed< double, 3 >  candidateTranslation = fixed_keypoint->location_ - moving_keypoint->location_;

      // if the current vote results in a candidate translation which would put matched point
      // outside the bounding box of the moving points, discard the vote
      itk::Point< double, 3 >  movingPoint;
      movingPoint[0] = fixed_keypoint->location_[0] - candidateTranslation[0];
      movingPoint[1] = fixed_keypoint->location_[1] - candidateTranslation[1];
      movingPoint[2] = fixed_keypoint->location_[2] - candidateTranslation[2];
      if( !movingBoundingBox->IsInside( movingPoint ) ) {
        std::cout << "DISCARDING MATCH" << std::endl;
        continue;
      }

      double descriptor_distance = (fixed_descriptor - moving_descriptor).magnitude();

#if !SORT_ON_RATIO
      match  one_match( moving_keypoint, moving_descriptor, fixed_keypoint, fixed_descriptor, descriptor_distance );
#else // store ratio instead of distance (the fields is used only in sorting, so it's okay
      vnl_vector< float >  moving_descriptor_second_closest = movingDescriptors[ near_neighbor_indices[1] ];
      double distance_to_second = (fixed_descriptor - moving_descriptor_second_closest).magnitude();
      match  one_match( moving_keypoint, moving_descriptor, fixed_keypoint, fixed_descriptor, descriptor_distance / distance_to_second );
#endif
      matches.push_back( one_match );

      //if( l > 20 ) break;
    }

    std::cout << "good keypoints: " << goodKeypoints << " out of: " << fixedKeypointNeighbors.size() << " first good: " << distanceToFirstGoodKeypoint << std::endl;



    // ******************
    // Descriptor testing
    // ******************

    vcl_sort( matches.begin(), matches.end(), DistanceAtIndexIsSmaller );

    //vcl_cout << "Ratios: ";
    //for( unsigned int b = 0; b < matches.size(); ++b ) {
    //  vcl_cout << matches[b].distance_ << " ";
    //}
    //vcl_cout << vcl_endl;


// consider only first good?
// if so, all matches within the neighborhood are tested, not just first (20) numberOfMatchesTest
#define CONSIDER_FIRST_ONLY 0

// write matches to files?
#define WRITE_MATCHES 1

    bool goodMatchFound = false;
    bool goodMatchIn10Found = false;

    vcl_vector< match >::size_type  index_to_first_good = 0;
    // after this for loop is finished, fixed_keypoint and moving_keypoint will contain the closest match in terms of descriptor distance
    cdcl_keypoint< 3 >::sptr  fixed_keypoint;
    cdcl_keypoint< 3 >::sptr  moving_keypoint;
    KeypointVectorType  fixedKeypointMapped;
    vnl_vector_fixed< double, 3 >  fixedNormalMappedVnl;
    vnl_vector_fixed< double, 3 >  fixedBinormalMappedVnl;

#if !CONSIDER_FIRST_ONLY
    vcl_vector< match >::size_type  numberOfMatchesTest = 20;
    numberOfMatchesTest = vcl_min( numberOfMatchesTest, matches.size() );
#else
    vcl_vector< match >::size_type  numberOfMatchesTest = matches.size();
#endif
    totalMatchesTested += numberOfMatchesTest;
    for( vcl_vector< match >::size_type  f = 0; f < numberOfMatchesTest; ++f ) {
      vcl_vector< match >::size_type  c = numberOfMatchesTest - f - 1;  // if this is in for statement and unsigned int, then there is overflow to positive numbers
      // select the keypoints as the match with the smallest descriptor distance
      fixed_keypoint = matches[c].fixed_keypoint_;
      moving_keypoint = matches[c].moving_keypoint_;

      //for( vcl_vector< match >::size_type  s = 0; s < 10 && s < matches.size(); ++s ) {
      //  //vcl_cout << matches[s].fixed_keypoint_->location_ << "  " << matches[s].moving_keypoint_->location_ << "  " << matches[s].distance_ << "\t";
      //  vcl_cout << matches[s].distance_ << "\t";
      //}
      //vcl_cout << vcl_endl;

      std::cout << c << ": fixed: " << fixed_keypoint->location_
                << "  moving: " << moving_keypoint->location_;



      itk::Point< double, 3 >  movingPoint;
      movingPoint[0] = moving_keypoint->location_[0];
      movingPoint[1] = moving_keypoint->location_[1];
      movingPoint[2] = moving_keypoint->location_[2];
      if( !movingBoundingBox->IsInside( movingPoint ) ) {
        std::cout << "RESULT: " << movingPoint << " IS OUT OF BOUNDS: " << minMovingPoint << "  " << maxMovingPoint << std::endl;
      }

      // now analyze the keypoint match
      //
      // moving image was obtained by deforming the fixed image, so deform the fixed keypoint with the same transform
      // and compare the two points
      //BSplineTransformType::OutputPointType  movingMapped = bsplineTransform->TransformPoint( moving_keypoint->location_.data_block() );
      vnl_vector_fixed< double, 3 >  fixedMappedVnl;
      vnl_vector< double >  fixed_normal   = fixed_keypoint->normal_;
      vnl_vector< double >  fixed_binormal = fixed_keypoint->binormal_;
      if( deformationType != NONE ) {

        // mapping error
        vnl_vector_fixed< double, 3 >  distance;
        if( deformationType == BSPLINE ) {
          // moving image was obtained by deforming the fixed image, so deform the fixed keypoint with the same transform
          // and compare the two points
          BSplineTransformType::OutputPointType  fixedMapped = bsplineTransform->TransformPoint( fixed_keypoint->location_.data_block() );
          fixedMappedVnl = fixedMapped.GetVnlVector();

          // map normal and binormal
          BSplineTransformType::OutputPointType  fixedMappedNor = bsplineTransform->TransformPoint( ( fixed_keypoint->location_ + fixed_normal ).data_block() );
          BSplineTransformType::OutputPointType  fixedMappedBin = bsplineTransform->TransformPoint( ( fixed_keypoint->location_ + fixed_binormal ).data_block() );

          fixedNormalMappedVnl = ( fixedMappedNor.GetVnlVector() - fixedMappedVnl ).normalize();
          fixedBinormalMappedVnl = ( fixedMappedBin.GetVnlVector() - fixedMappedVnl ).normalize();
        }
        else if( deformationType == FIELD ) {

          // map point
          //
          // deformation field was obtained from registration 
          itk::Point< double, 3 >  fixedPoint;
          fixedPoint[0] = fixed_keypoint->location_[0];
          fixedPoint[1] = fixed_keypoint->location_[1];
          fixedPoint[2] = fixed_keypoint->location_[2];

          // get the deformation vector from the field
          DeformationFieldType::IndexType  fixedIndex;
          deformationField->TransformPhysicalPointToIndex( fixedPoint, fixedIndex );
          VectorPixelType  deformationVector = deformationField->GetPixel( fixedIndex );

          itk::Point< double, 3 >  fixedPointMapped = fixedPoint + deformationVector;
          fixedMappedVnl = fixedPointMapped.GetVnlVector();

          // map normal
          //
          // deformation field was obtained from registration 
          itk::Point< double, 3 >  fixedPointNor;
          fixedPointNor[0] = fixed_keypoint->location_[0] + fixed_normal[0];
          fixedPointNor[1] = fixed_keypoint->location_[1] + fixed_normal[1];
          fixedPointNor[2] = fixed_keypoint->location_[2] + fixed_normal[2];

          // get the deformation vector from the field
          DeformationFieldType::IndexType  fixedIndexNor;
          deformationField->TransformPhysicalPointToIndex( fixedPointNor, fixedIndexNor );
          VectorPixelType  deformationVectorNor = deformationField->GetPixel( fixedIndexNor );

          itk::Point< double, 3 >  fixedPointNorMapped = fixedPointNor + deformationVectorNor;
          fixedNormalMappedVnl = ( fixedPointNorMapped.GetVnlVector() - fixedMappedVnl ).normalize();


          // map binormal
          //
          // deformation field was obtained from registration 
          itk::Point< double, 3 >  fixedPointBin;
          fixedPointBin[0] = fixed_keypoint->location_[0] + fixed_binormal[0];
          fixedPointBin[1] = fixed_keypoint->location_[1] + fixed_binormal[1];
          fixedPointBin[2] = fixed_keypoint->location_[2] + fixed_binormal[2];

          // get the deformation vector from the field
          DeformationFieldType::IndexType  fixedIndexBin;
          deformationField->TransformPhysicalPointToIndex( fixedPointBin, fixedIndexBin );
          VectorPixelType  deformationVectorBin = deformationField->GetPixel( fixedIndexBin );

          itk::Point< double, 3 >  fixedPointBinMapped = fixedPointBin + deformationVectorBin;
          fixedBinormalMappedVnl = ( fixedPointBinMapped.GetVnlVector() - fixedMappedVnl ).normalize();

        }
      }


      // find moving keypoints which are closest to the mapped fixed keypoint
      for( unsigned int d = 0; d < 3; ++d ) {
        fixedKeypointMapped.SetElement( d, fixedMappedVnl[d] );
      }

      vnl_vector_fixed< double, 3 >  distance = fixedMappedVnl - moving_keypoint->location_;
      //vnl_vector< float >  difference = moving_descriptor - fixed_descriptor;
      std::cout << " Dist: " << distance.magnitude();

      itk::Point< double, 3 >  fixedPointMapped;  
      fixedPointMapped[0] = fixedMappedVnl[0];
      fixedPointMapped[1] = fixedMappedVnl[1];
      fixedPointMapped[2] = fixedMappedVnl[2];
      bool insideBB = movingBoundingBox->IsInside( fixedPointMapped );
      if( !insideBB ) std::cout << "  loc. outside overlap!" << std::endl;
      else std::cout << std::endl;

      // compute histogram values for normals, binormals and distances; use descriptor match
      double normal_difference = vcl_acos( dot_product( moving_keypoint->normal_, fixedNormalMappedVnl ) ) * 180.0 / vnl_math::pi;
      unsigned int normal_bin_number = vnl_math_rnd( normal_difference / 10.0 );

      double binormal_difference = vcl_acos( dot_product( moving_keypoint->binormal_, fixedBinormalMappedVnl ) ) * 180.0 / vnl_math::pi;
      unsigned int binormal_bin_number = vnl_math_rnd( binormal_difference / 10.0 );

      //std::cout << difference.magnitude() << "\t";
      if( c == 0 ) {
        if( distance.magnitude() < goodnessThreshold 
          && normal_difference < 20.0
          && binormal_difference < 20.0 ) ++goodAtFirst;

        unsigned int distance_bin_number;
        if( distance.magnitude() > 9.0 ) {
          distance_bin_number = 3;
        }
        else {
          distance_bin_number = (unsigned int) ( vcl_floor( distance.magnitude() / 3.0 ) );
        }

        ++normals_distances[normal_bin_number][distance_bin_number];
        ++binormals_distances[binormal_bin_number][distance_bin_number];
      }


      if( normal_difference  > 45 || binormal_difference > 45 ) {
        if( c == 0 ) ++badOrientationsAtFirst;
        else ++badOrientationsAtAll;
      }
      else
      if( distance.magnitude() < goodnessThreshold 
          && normal_difference < 20.0
          && binormal_difference < 20.0 ) {

        index_to_first_good = c;

        goodMatchesFile << vcl_setw( 6 ) << vcl_setfill( '0' ) << t << "_" << vcl_setw( 2 ) << c << vcl_endl;

        if( !goodMatchFound ) {
               ++goodAtAll;
               goodMatchFound = true;
         }

         if( !goodMatchIn10Found && c < 10 ) {
             ++goodIn10;
             goodMatchIn10Found = true;
         }


#if !CONSIDER_FIRST_ONLY
         vnl_vector_fixed< double, 3 >  vectorToFirstGood;
         vectorToFirstGood[0] = fixed_keypoint->location_[0] - queryPoint[0];
         vectorToFirstGood[1] = fixed_keypoint->location_[1] - queryPoint[1];
         vectorToFirstGood[2] = fixed_keypoint->location_[2] - queryPoint[2];
         unsigned int distanceToGood = (unsigned int) ( vnl_math_rnd( vectorToFirstGood.magnitude() ) );
         if( distanceToGood >= distanceToFirstGood.size() ) distanceToFirstGood.resize( distanceToGood+1, 0 );
         ++distanceToFirstGood[distanceToGood];
#endif

      }

#if WRITE_MATCHES
      // save the current match
      //
      // save moving keypoint descriptor
      vcl_ostringstream  movingDescriptorFile;
      movingDescriptorFile << matchesDirPath << "//moving_descriptor" << vcl_setw( 6 ) << vcl_setfill( '0' ) << t << "_" << vcl_setw( 2 ) << c << ".vtk";

      vcl_vector< vbl_smart_ptr< cdcl_keypoint< 3 > > >  movingKeypointOne;
      movingKeypointOne.push_back( moving_keypoint );
      vcl_vector< vnl_vector< float > >  movingDescriptorOne;
      movingDescriptorOne.push_back( matches[c].moving_descriptor_ );

      cdcl_write_keypoint_descriptors_VTK( movingKeypointOne, movingDescriptorOne, movingDescriptorFile.str() );

      // save fixed keypoint descriptor
      vcl_ostringstream  fixedDescriptorFile;
      fixedDescriptorFile << matchesDirPath << "//fixed_descriptor" << vcl_setw( 6 ) << vcl_setfill( '0' ) << t << "_" << vcl_setw( 2 ) << c << ".vtk";

      vcl_vector< vbl_smart_ptr< cdcl_keypoint< 3 > > >  fixedKeypointOne;
      fixedKeypointOne.push_back( fixed_keypoint );
      vcl_vector< vnl_vector< float > >  fixedDescriptorOne;
      fixedDescriptorOne.push_back( matches[c].fixed_descriptor_ );

      cdcl_write_keypoint_descriptors_VTK( fixedKeypointOne, fixedDescriptorOne, fixedDescriptorFile.str() );

      // save clicked query location
      vcl_ostringstream  queryLocationFileName;
      queryLocationFileName << matchesDirPath << "//query_location" << vcl_setw( 6 ) << vcl_setfill( '0' ) << t << "_" << vcl_setw( 2 ) << c << ".txt";
      vcl_ofstream  queryLocationFile( queryLocationFileName.str().c_str() );
      queryLocationFile << queryPoint[0] << " " << queryPoint[1] << " " << queryPoint[2] << vcl_endl;
      queryLocationFile.close();



      std::cout << "Volumes have been read so outputing keypoint neighborhoods." << std::endl;
      if( volumesRead ) {
        vcl_ostringstream  movingImageFile;
        movingImageFile << matchesDirPath << "//moving_image" << vcl_setw( 6 ) << vcl_setfill( '0' ) << t << "_" << vcl_setw( 2 ) << c << ".png";      

        std::cout << "Moving keypoint wrote: " << movingImageFile.str() << std::endl;
        WriteROIImage<MovingImageType>( movingImage, moving_keypoint, movingImageFile.str() );

        vcl_ostringstream  fixedImageFile;
        fixedImageFile << matchesDirPath << "//fixed_image" << vcl_setw( 6 ) << vcl_setfill( '0' ) << t << "_" << vcl_setw( 2 ) << c << ".png";      

        std::cout << "Fixed keypoint wrote: " << fixedImageFile.str() << std::endl;
        WriteROIImage<MovingImageType>( fixedImage, fixed_keypoint, fixedImageFile.str() );

        vcl_ostringstream  queryImageFile;
        queryImageFile << matchesDirPath << "//query_image_fixed" << vcl_setw( 6 ) << vcl_setfill( '0' ) << t << "_" << vcl_setw( 2 ) << c << ".png";      

        std::cout << "Query wrote: " << queryImageFile.str() << std::endl;
        itk::Point< double, 3 >  queryPt;
        queryPt[0] = queryPoint[0];
        queryPt[1] = queryPoint[1];
        queryPt[2] = queryPoint[2];
        WriteROIImage<MovingImageType>( fixedImage, queryPt, queryImageFile.str() );
       
        //vcl_ostringstream  queryImageMovingFile;
        //queryImageMovingFile << matchesDirPath << "//query_image_moving" << vcl_setw( 6 ) << vcl_setfill( '0' ) << t << ".png";      
        //WriteROIImage<MovingImageType>( movingImage, queryPt, queryImageMovingFile.str() );
      }
#endif
    }


#if CONSIDER_FIRST_ONLY
    if( goodMatchFound ) {
      cdcl_keypoint< 3 >::sptr  fixed_keypoint_first_good = matches[index_to_first_good].fixed_keypoint_;

      vnl_vector_fixed< double, 3 >  vectorToFirstGood;
      vectorToFirstGood[0] = fixed_keypoint_first_good->location_[0] - queryPoint[0];
      vectorToFirstGood[1] = fixed_keypoint_first_good->location_[1] - queryPoint[1];
      vectorToFirstGood[2] = fixed_keypoint_first_good->location_[2] - queryPoint[2];
      unsigned int distanceToGood = (unsigned int) ( vnl_math_rnd( vectorToFirstGood.magnitude() ) );
      if( distanceToGood >= distanceToFirstGood.size() ) distanceToFirstGood.resize( distanceToGood+1, 0 );
      ++distanceToFirstGood[distanceToGood];
    }
#endif

    timerOneClick.Stop();
    std::cout << "Answer for one click #" << t << " in: " << timerOneClick.GetMeanTime() << " sec." << std::endl;


    // fixed_keypoint and moving_keypoint are now from the first match

    //std::cout << "Volumes have been read so outputing keypoint neighborhoods." << std::endl;
    //if( volumesRead ) {
    //  vcl_ostringstream  movingImageFile;
    //  movingImageFile << matchesDirPath << "//moving_image" << vcl_setw( 6 ) << vcl_setfill( '0' ) << t << ".png";      

    //  std::cout << "Moving keypoint wrote: " << movingImageFile.str() << std::endl;
    //  WriteROIImage<MovingImageType>( movingImage, moving_keypoint, movingImageFile.str() );

    //  vcl_ostringstream  fixedImageFile;
    //  fixedImageFile << matchesDirPath << "//fixed_image" << vcl_setw( 6 ) << vcl_setfill( '0' ) << t << ".png";      

    //  std::cout << "Fixed keypoint wrote: " << fixedImageFile.str() << std::endl;
    //  WriteROIImage<MovingImageType>( fixedImage, fixed_keypoint, fixedImageFile.str() );

    //  vcl_ostringstream  queryImageFile;
    //  queryImageFile << matchesDirPath << "//query_image_fixed" << vcl_setw( 6 ) << vcl_setfill( '0' ) << t << ".png";      

    //  std::cout << "Query wrote: " << queryImageFile.str() << std::endl;
    //  itk::Point< double, 3 >  queryPt;
    //  queryPt[0] = queryPoint[0];
    //  queryPt[1] = queryPoint[1];
    //  queryPt[2] = queryPoint[2];
    //  WriteROIImage<MovingImageType>( fixedImage, queryPt, queryImageFile.str() );
    // 
    //  //vcl_ostringstream  queryImageMovingFile;
    //  //queryImageMovingFile << matchesDirPath << "//query_image_moving" << vcl_setw( 6 ) << vcl_setfill( '0' ) << t << ".png";      
    //  //WriteROIImage<MovingImageType>( movingImage, queryPt, queryImageMovingFile.str() );
    //}


    // output normal and binormal computed from the covariance matrices (as is done in the descriptor computation)
    //
    //// for decomposition to get the eigenvectors from the covariances
    //vnl_matrix< double >  V( 3, 3 );
    //vnl_vector< double >  D( 3 );
    //const vnl_matrix< double >  fixed_covariance = fixed_keypoint->covariance_;
    //vnl_symmetric_eigensystem_compute( fixed_covariance, V, D );

    //vnl_vector< double >  fixed_normal   = V.get_column( 0 );
    //vnl_vector< double >  fixed_binormal = V.get_column( 1 );

    //// limiting to half space for keypoint normal (ellipsoid is mirrored and from the keypoint point of view
    //// the two halves have the same orientation)
    //if( fixed_normal[2] < 0.0 ) fixed_normal[2] *= -1.0;

    //// limiting to half space for keypoint binormal (ellipsoid is mirrored and from the keypoint point of view
    //// the two halves have the same orientation)
    //if( fixed_binormal[1] < 0.0 ) fixed_binormal[1] *= -1.0;

    //vnl_vector< double >  fixed_normal = fixed_keypoint->normal_;
    //vnl_vector< double >  fixed_binormal = fixed_keypoint->binormal_;

    //std::cout << "Fixed: nor: " << fixed_normal << " binor: " << fixed_binormal << std::endl;

    unsigned int numberOfMappedKeypointNeighbors = 1;
    KeypointTreeType::InstanceIdentifierVectorType  keypointNeighbors;
    kdTreeMovingKeypoints->Search( fixedKeypointMapped, numberOfMappedKeypointNeighbors, keypointNeighbors );
    for( unsigned int l = 0; l < numberOfMappedKeypointNeighbors; ++l ) {
      KeypointVectorType  closestKeypoint = kdTreeMovingKeypoints->GetMeasurementVector( keypointNeighbors[l] );
      vnl_vector_fixed< float, 3 >  candidate_distance = closestKeypoint.GetVnlVector() - fixedKeypointMapped.GetVnlVector();

#if USE_ALSO_ITK_KDTREE
      DescriptorVectorType  closestDescriptor = kdTreeMovingDescriptors->GetMeasurementVector( keypointNeighbors[l] );
#else
      const vnl_vector< float > &  closestDescriptor = movingDescriptors[ keypointNeighbors[l] ];
#endif

      std::cout << "Candidate " << l << ": " << candidate_distance.magnitude() << " ";

      cdcl_keypoint<3>::sptr  moving_keypoint_closest = movingKeypoints[keypointNeighbors[l]];

      vnl_vector< double >  moving_normal   = moving_keypoint_closest->normal_;
      vnl_vector< double >  moving_binormal = moving_keypoint_closest->binormal_;

      if( volumesRead ) {
        vcl_ostringstream  movingImageFile;
        movingImageFile << matchesDirPath << "//moving_image_closest" << vcl_setw( 6 ) << vcl_setfill( '0' ) << t << ".png";      

        WriteROIImage<MovingImageType>( movingImage, moving_keypoint_closest, movingImageFile.str() );
      }

      //const vnl_matrix< double >  moving_covariance = moving_keypoint->covariance_;
      //vnl_symmetric_eigensystem_compute( moving_covariance, V, D );

      //vnl_vector< double >  moving_normal   = V.get_column( 0 );
      //vnl_vector< double >  moving_binormal = V.get_column( 1 );

      //// limiting to half space for keypoint normal (ellipsoid is mirrored and from the keypoint point of view
      //// the two halves have the same orientation) -> flip the normal
      //if( moving_normal[2] < 0.0 ) moving_normal *= -1.0;

      //// limiting to half space for keypoint binormal (ellipsoid is mirrored and from the keypoint point of view
      //// the two halves have the same orientation) -> flip the binormal
      //if( moving_binormal[1] < 0.0 ) moving_binormal *= -1.0;

      //std::cout << "Fixed: nor: " << vcl_setprecision( 3 ) << fixed_normal << " binor: " << fixed_binormal << std::endl;


      // compute histogram values for normals, binormals and distances; use closest keypoint
      double normal_difference = vcl_acos( dot_product( moving_normal, fixedNormalMappedVnl ) ) * 180.0 / vnl_math::pi;
      unsigned int normal_bin_number = vcl_floor( normal_difference / 10.0 );

      double binormal_difference = vcl_acos( dot_product( moving_binormal, fixedBinormalMappedVnl ) ) * 180.0 / vnl_math::pi;
      unsigned int binormal_bin_number = vcl_floor( binormal_difference / 10.0 );

      if( candidate_distance.magnitude() < keypointGoodnessThreshold 
          && normal_difference < 20.0
          && binormal_difference < 20.0 ) ++goodCandidateKeypoints;

      unsigned int distance_bin_number;
      if( candidate_distance.magnitude() > 9.0 ) {
        distance_bin_number = 3;
      }
      else {
        distance_bin_number = (unsigned int) ( vcl_floor( candidate_distance.magnitude() / 3.0 ) );
      }

      ++keypoints_normals_distances[normal_bin_number][distance_bin_number];
      ++keypoints_binormals_distances[binormal_bin_number][distance_bin_number];


#if OUTPUT_IMAGES
      // write the candidate descriptor to the image
      unsigned int i;
      for( descItr.GoToBegin(), i = 0; ! descItr.IsAtEnd() && i < closestDescriptor.size(); ++descItr, i+=3 ) {
           descItr.Set( closestDescriptor[i] );
          //std::cout << closestDescriptor[i] << " ";
      }
      for( i = 1; ! descItr.IsAtEnd() && i < closestDescriptor.size(); ++descItr, i+=3 ) {
           descItr.Set( closestDescriptor[i] );
          //std::cout << closestDescriptor[i] << " ";
      }
      for( i = 2; ! descItr.IsAtEnd() && i < closestDescriptor.size(); ++descItr, i+=3 ) {
           descItr.Set( closestDescriptor[i] );
          //std::cout << closestDescriptor[i] << " ";
      }
      descIm->Modified();

      vcl_ostringstream  fileNameMoving;
      fileNameMoving << "desc_candidate" << vcl_setw( 6 ) << vcl_setfill( '0' ) << t << ".png";      
      //vcl_cout << "Saving " << fileNameMoving.str() << vcl_endl;

      writer->SetFileName( fileNameMoving.str().c_str() );
      writer->Update();
#endif
    }
    //std::cout << std::endl;


#if OUTPUT_IMAGES
    // write the descriptor to the image
    unsigned int i;
    for( descItr.GoToBegin(), i = 0; ! descItr.IsAtEnd() && i < moving_descriptor.size(); ++descItr, i+=3 ) {
         descItr.Set( moving_descriptor[i] );
        //std::cout << moving_descriptor[i] << " ";
    }
    for( i = 1; ! descItr.IsAtEnd() && i < moving_descriptor.size(); ++descItr, i+=3 ) {
         descItr.Set( moving_descriptor[i] );
        //std::cout << moving_descriptor[i] << " ";
    }
    for( i = 2; ! descItr.IsAtEnd() && i < moving_descriptor.size(); ++descItr, i+=3 ) {
         descItr.Set( moving_descriptor[i] );
        //std::cout << moving_descriptor[i] << " ";
    }
    descIm->Modified();

    vcl_ostringstream  fileNameMoving;
    fileNameMoving << "desc_moving" << vcl_setw( 6 ) << vcl_setfill( '0' ) << t << ".png";
    //vcl_cout << "Saving " << fileNameMoving.str() << vcl_endl;

    writer->SetFileName( fileNameMoving.str().c_str() );
    writer->Update();

    unsigned int ind;
    for( descItr.GoToBegin(), ind = 0; ! descItr.IsAtEnd() && ind < fixed_descriptor.size(); ++descItr, ind+=3 ) {
        descItr.Set( fixed_descriptor[ind] );
        //std::cout << fixed_descriptor[ind] << " ";
    }
    for( ind = 1; ! descItr.IsAtEnd() && ind < fixed_descriptor.size(); ++descItr, ind+=3 ) {
        descItr.Set( fixed_descriptor[ind] );
        //std::cout << fixed_descriptor[ind] << " ";
    }
    for( ind = 2; ! descItr.IsAtEnd() && ind < fixed_descriptor.size(); ++descItr, ind+=3 ) {
        descItr.Set( fixed_descriptor[ind] );
        //std::cout << fixed_descriptor[ind] << " ";
    }
    descIm->Modified();

    vcl_ostringstream  fileNameFixed;
    fileNameFixed << "desc_fixed" << vcl_setw( 6 ) << vcl_setfill( '0' ) << t << ".png";
    //vcl_cout << "Saving " << fileNameFixed.str() << vcl_endl;

    writer->SetFileName( fileNameFixed.str().c_str() );
    writer->Update();
#endif

  }

  if( locationsFile.set() ) {
    noduleFile.close();
  }
  goodMatchesFile.close();

  std::cout << "moving keypoints: " << movingKeypoints.size() << std::endl
            << "fixed keypoints: " << fixedKeypoints.size() << std::endl
            << "keypoints tested: " << numberOfKeypointsTested << std::endl
            << "good matches on the first hit: " << goodAtFirst << std::endl
            << "good matches in 10 hits: " << goodIn10 << std::endl
            << "good matches in X hits: " << goodAtAll << std::endl
            << "good keypoints (tolerance " << keypointGoodnessThreshold << ", 20, 20): " << goodCandidateKeypoints << std::endl
            << "bad orientations on the first hit: " << badOrientationsAtFirst << std::endl
            << "bad orientations in X hits: " << badOrientationsAtAll << std::endl;

  // histograms of failures
  for( unsigned int a = 0; a < num_angles; ++a ) {
    vcl_cout << a*10 << "-" << (a+1)*10 << "\t";
  }
  vcl_cout << vcl_endl;

  vcl_cout << "nor_dist: " << vcl_endl;
  for( unsigned int d = 0; d < num_distances; ++d ) {
    for( unsigned int a = 0; a < num_angles; ++a ) {
      vcl_cout << normals_distances[a][d] << "\t";
    }
    vcl_cout << vcl_endl;
  }
 
  vcl_cout << "binor_dist: " << vcl_endl;
  for( unsigned int d = 0; d < num_distances; ++d ) {
    for( unsigned int a = 0; a < num_angles; ++a ) {
      vcl_cout << binormals_distances[a][d] << "\t";
    }
    vcl_cout << vcl_endl;
  }

  // histograms of failures
  vcl_cout << "keypoints nor_dist: " << vcl_endl;
  for( unsigned int d = 0; d < num_distances; ++d ) {
    for( unsigned int a = 0; a < num_angles; ++a ) {
      vcl_cout << keypoints_normals_distances[a][d] << "\t";
    }
    vcl_cout << vcl_endl;
  }
 
  vcl_cout << "keypoints binor_dist: " << vcl_endl;
  for( unsigned int d = 0; d < num_distances; ++d ) {
    for( unsigned int a = 0; a < num_angles; ++a ) {
      vcl_cout << keypoints_binormals_distances[a][d] << "\t";
    }
    vcl_cout << vcl_endl;
  }

  vcl_cout << "Distance to first good: ";
  unsigned int binSum = 0;
  for( vcl_vector< unsigned int >::size_type  n = 0; n < distanceToFirstGood.size(); ++n ) {
    binSum += distanceToFirstGood[n];
    if( n % 10 == 0 ) {
      vcl_cout << binSum << " ";
    //  binSum = 0;
    }
  }
  vcl_cout << vcl_endl;
  vcl_cout << "Max Distance: " << distanceToFirstGood.size() - 1 << vcl_endl;
  vcl_cout << "Total matches tested: " << totalMatchesTested << vcl_endl;

  return 0;
}
