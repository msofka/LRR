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
#include "itkBoundingBox.h"

#include "itkImageSeriesReader.h"
#include "itkImageSeriesWriter.h"

#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"

#include "itkKdTree.h"
#include "itkKdTreeGenerator.h"
#include "itkListSample.h"

#include "vtkXMLPolyDataReader.h"
#include <cdcl/io/vtkPolyDataToFeaturesWithShapeFilter.h>

#include <cdcl/cdcl_feature.h>
#include <cdcl/cdcl_utils_VTK.h>

#include <vnl/algo/vnl_symmetric_eigensystem.h>
#include <vnl/algo/vnl_generalized_eigensystem.h>
#include <vnl/vnl_random.h>

#include <vcl_utility.h>
#include <vcl_algorithm.h>


//:
// \file
// \brief  Testing repeatability of features.
// \author Michal Sofka
// \date   April 2008


// write debugging image chips for each descriptor
#define OUTPUT_IMAGES 0


namespace {

template <class ImageType>
void WriteROIImage( typename ImageType::Pointer const &  Image,
                    cdcl_keypoint< 3 >::sptr const &  Keypoint,
                    vcl_string  FileName )
{
  typedef typename itk::ImageRegion<3>  RegionOfInterestType;
  typedef itk::Point< double, 3 >  PointType;
  
  PointType  point;
  point[0] = Keypoint->location_[0];
  point[1] = Keypoint->location_[1];
  point[2] = Keypoint->location_[2];

  RegionOfInterestType  m_ROI;

  RegionOfInterestType::SizeType  size;
  RegionOfInterestType::IndexType  index;
  Image->TransformPhysicalPointToIndex( point, index );

  typename ImageType::SpacingType  spacing = Image->GetSpacing();

  size[0] = 50 / spacing[0];
  size[1] = 50 / spacing[1];
  size[2] = 50 / spacing[2];

  m_ROI.SetSize( size );
  m_ROI.SetIndex( index );

  m_ROI.Crop( Image->GetLargestPossibleRegion() );



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

  //  Writer
  WriterType::Pointer writer = WriterType::New();
  //writer->SetInput( m_SlicesImage );
  writer->SetInput( rescaler->GetOutput() );

  writer->SetFileName( FileName.c_str() );
  writer->Update();

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
  // e.g. feature_repeatability.exe 241/4804/000002cropped_00.vtk 241/4802/000002cropped_00.vtk 241/4802_000002cropped-4804_000002cropped.vtk
  vul_arg< const char* >        movingFeatureFile     ( 0, "VTK moving feature file (of a synthetically deformed image)" );
  vul_arg< const char* >        fixedFeatureFile      ( 0, "VTK fixed feature file" );
  vul_arg< const char* >        transformFile         ( "-trans", "ITK transform file (transform used to warp points)", "" );

  vul_arg_parse( argc, argv );

  vcl_cout << "Command: ";
  for( int i = 0; i < argc; ++i )
    vcl_cout << argv[i] << " ";
  vcl_cout << vcl_endl;



  // Read features
  //
  vcl_vector< vcl_vector< cdcl_feature_with_shape< 3 >::sptr > >  m_FeaturesFixedAll;
  vcl_vector< vcl_vector< cdcl_feature_with_shape< 3 >::sptr > >  m_FeaturesMovingAll;

  {
  vtkSmartPointer< vtkXMLPolyDataReader >  featureReader = vtkSmartPointer< vtkXMLPolyDataReader >::New();
  vcl_string  filename = vcl_string( movingFeatureFile() );
  featureReader->SetFileName( filename.c_str() );
  featureReader->Update();  // need this update (bug in the converting filter)

  // create the converting filter
  vtkSmartPointer< vtkPolyDataToFeaturesWithShapeFilter >  polyDataToFeaturesFilter = vtkSmartPointer< vtkPolyDataToFeaturesWithShapeFilter >::New();
  polyDataToFeaturesFilter->SetInput( featureReader->GetOutput() );
  polyDataToFeaturesFilter->Update();
  vtkSmartPointer< vtkFeatureWithShapeAttributeSet >  featureAttributeSetVTK = polyDataToFeaturesFilter->GetOutput();

  typedef vtkFeatureWithShapeAttributeSet::FeatureSetType  FeatureSetType;
  FeatureSetType &  features = featureAttributeSetVTK->GetPoints();
  m_FeaturesMovingAll.push_back( features );
  }

  {
  vtkSmartPointer< vtkXMLPolyDataReader >  featureReader = vtkSmartPointer< vtkXMLPolyDataReader >::New();
  vcl_string  filename = vcl_string( fixedFeatureFile() );
  featureReader->SetFileName( filename.c_str() );
  featureReader->Update();  // need this update (bug in the converting filter)

  // create the converting filter
  vtkSmartPointer< vtkPolyDataToFeaturesWithShapeFilter >  polyDataToFeaturesFilter = vtkSmartPointer< vtkPolyDataToFeaturesWithShapeFilter >::New();
  polyDataToFeaturesFilter->SetInput( featureReader->GetOutput() );
  polyDataToFeaturesFilter->Update();
  vtkSmartPointer< vtkFeatureWithShapeAttributeSet >  featureAttributeSetVTK = polyDataToFeaturesFilter->GetOutput();

  typedef vtkFeatureWithShapeAttributeSet::FeatureSetType  FeatureSetType;
  FeatureSetType &  features = featureAttributeSetVTK->GetPoints();
  m_FeaturesFixedAll.push_back( features );
  }




  // Kd-tree for moving feature locations
  typedef itk::Vector< float, 3 >                                FeatureVectorType;
  typedef itk::Statistics::ListSample< FeatureVectorType >       FeatureSampleType;
  typedef itk::Statistics::KdTreeGenerator< FeatureSampleType >  FeatureTreeGeneratorType;
  typedef FeatureTreeGeneratorType::KdTreeType                   FeatureTreeType;



  // setup indexing
  FeatureSampleType::Pointer    sampleMovingFeatures = FeatureSampleType::New();

  // moving feature vectors
  //
  sampleMovingFeatures->SetMeasurementVectorSize( m_FeaturesMovingAll[0].size() );

  FeatureVectorType  mv;
  for( unsigned int n = 0 ; n < m_FeaturesMovingAll[0].size(); ++n ) {
    for( unsigned int d = 0; d < m_FeaturesMovingAll[0][n]->location_.size(); ++d ) {
      mv.SetElement( d, m_FeaturesMovingAll[0][n]->location_[d] );
    }
    sampleMovingFeatures->PushBack( mv );
  }
  
  //// recommended hack on releasing memory is to swap in empty vector
  //vcl_vector< vnl_vector< float > >  empty_vector;
  //empty_vector.swap( m_FeaturesMovingAll );


  // moving feature locations
  FeatureTreeGeneratorType::Pointer  featureTreeGenerator = FeatureTreeGeneratorType::New();

  featureTreeGenerator->SetSample( sampleMovingFeatures );
  featureTreeGenerator->SetBucketSize( 16 );
  featureTreeGenerator->Update();

  FeatureTreeType::Pointer  kdTreeMovingFeatures = featureTreeGenerator->GetOutput();


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
    vcl_string  transformExtension = itksys::SystemTools::GetFilenameExtension( transformFile() );
    //std::cout << "Transform extension: " << transformExtension << std::endl;
    // if warp field supplied
    if( transformExtension == ".vtk" || transformExtension == ".mhd" ) {
      deformationType = FIELD;

      typedef   itk::ImageFileReader< DeformationFieldType >  FieldReaderType;

      std::cout << "Reading deformation field: " << transformFile();

      FieldReaderType::Pointer fieldReader = FieldReaderType::New();
      fieldReader->SetFileName( transformFile() );
      fieldReader->Update();
      deformationField = fieldReader->GetOutput();

      std::cout << " of size " << deformationField->GetLargestPossibleRegion().GetSize() << std::endl;
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



  //typedef  itk::BoundingBox< unsigned long, 3, double >  BoundingBoxType;
  //typedef BoundingBoxType::PointType  PointType;

  //BoundingBoxType::Pointer  fixedBoundingBox = BoundingBoxType::New();
  //BoundingBoxType::PointsContainerPointer  fixedPoints = BoundingBoxType::PointsContainer::New();
  //
  //for( vcl_vector< cdcl_keypoint< 3 >::sptr >::size_type  i = 0; i < fixedKeypoints.size(); ++i ) {
  //  PointType  point;
  //  point[0] = fixedKeypoints[i]->location_[0];
  //  point[1] = fixedKeypoints[i]->location_[1];
  //  point[2] = fixedKeypoints[i]->location_[2];
  //  fixedPoints->push_back( point );
  //}
  //fixedBoundingBox->SetPoints( fixedPoints );

  //PointType  minFixedPoint = fixedBoundingBox->GetMinimum();
  //PointType  maxFixedPoint = fixedBoundingBox->GetMaximum();



  //BoundingBoxType::Pointer  movingBoundingBox = BoundingBoxType::New();
  //BoundingBoxType::PointsContainerPointer  movingPoints = BoundingBoxType::PointsContainer::New();
  //
  //for( vcl_vector< cdcl_keypoint< 3 >::sptr >::size_type  i = 0; i < movingKeypoints.size(); ++i ) {
  //  PointType  point;
  //  point[0] = movingKeypoints[i]->location_[0];
  //  point[1] = movingKeypoints[i]->location_[1];
  //  point[2] = movingKeypoints[i]->location_[2];
  //  movingPoints->push_back( point );
  //}
  //movingBoundingBox->SetPoints( movingPoints );

  //PointType  minMovingPoint = movingBoundingBox->GetMinimum();
  //PointType  maxMovingPoint = movingBoundingBox->GetMaximum();



  srand( (unsigned)time( NULL ) );
  vnl_random rng;
  rng.reseed(333248);
  // Usually, you will want to generate a number in a specific range
  int RANGE_MIN = 0;
  int RANGE_MAX = m_FeaturesFixedAll.size()-1;


  // histograms of failures
  const unsigned int num_angles = 19;
  const unsigned int num_distances = 4;
  unsigned int normals_distances[num_angles][num_distances];
  for( unsigned int a = 0; a < num_angles; ++a )
    for( unsigned int d = 0; d < num_distances; ++d ) {
      normals_distances[a][d] = 0;
    }


  std::cout << "Computing feature repeatability ..." << std::endl;
  unsigned int goodFeatures = 0;

  float goodnessThreshold = 2.0f;

  for( unsigned int t = 0; t < m_FeaturesFixedAll[0].size(); ++t ) {

    cdcl_feature_with_shape< 3 >::sptr  fixed_feature = m_FeaturesFixedAll[0][t];

    vcl_vector< vnl_vector_fixed< double, 3 > >  fixedNormalsMappedVnl;
    
    // moving image was obtained by deforming the fixed image, so deform the fixed keypoint with the same transform
    // and compare the two points
    //BSplineTransformType::OutputPointType  movingMapped = bsplineTransform->TransformPoint( moving_keypoint->location_.data_block() );
    vnl_vector_fixed< double, 3 >  fixedMappedVnl;
    if( deformationType != NONE ) {

      for( unsigned int d = 0; d < fixed_feature->directions_.size(); ++ d ) {

        vnl_vector< double >  fixed_normal = fixed_feature->directions_[d];

        // mapping error
        vnl_vector_fixed< double, 3 >  distance;
        if( deformationType == BSPLINE ) {
          // moving image was obtained by deforming the fixed image, so deform the fixed keypoint with the same transform
          // and compare the two points
          BSplineTransformType::OutputPointType  fixedMapped = bsplineTransform->TransformPoint( fixed_feature->location_.data_block() );
          fixedMappedVnl = fixedMapped.GetVnlVector();

          // map normal
          BSplineTransformType::OutputPointType  fixedMappedNor = bsplineTransform->TransformPoint( ( fixed_feature->location_ + fixed_normal ).data_block() );

          fixedNormalsMappedVnl.push_back( ( fixedMappedNor.GetVnlVector() - fixedMappedVnl ).normalize() );
        }
        else if( deformationType == FIELD ) {

          // map point
          //
          // deformation field was obtained from registration 
          itk::Point< double, 3 >  fixedPoint;
          fixedPoint[0] = fixed_feature->location_[0];
          fixedPoint[1] = fixed_feature->location_[1];
          fixedPoint[2] = fixed_feature->location_[2];

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
          fixedPointNor[0] = fixed_feature->location_[0] + fixed_normal[0];
          fixedPointNor[1] = fixed_feature->location_[1] + fixed_normal[1];
          fixedPointNor[2] = fixed_feature->location_[2] + fixed_normal[2];

          // get the deformation vector from the field
          DeformationFieldType::IndexType  fixedIndexNor;
          deformationField->TransformPhysicalPointToIndex( fixedPointNor, fixedIndexNor );
          VectorPixelType  deformationVectorNor = deformationField->GetPixel( fixedIndexNor );

          itk::Point< double, 3 >  fixedPointNorMapped = fixedPointNor + deformationVectorNor;
          fixedNormalsMappedVnl.push_back( ( fixedPointNorMapped.GetVnlVector() - fixedMappedVnl ).normalize() );

        }
      }
    }
    //else {
    //  fixedMappedVnl = fixed_feature->location_;
    //}


    // find moving keypoints which are closest to the mapped fixed keypoint
    FeatureVectorType  fixedFeatureMapped;
    for( unsigned int d = 0; d < 3; ++d ) {
      fixedFeatureMapped.SetElement( d, fixedMappedVnl[d] );
    }

    unsigned int numberOfMappedFeatureNeighbors = 1;
    FeatureTreeType::InstanceIdentifierVectorType  featureNeighbors;
    kdTreeMovingFeatures->Search( fixedFeatureMapped, numberOfMappedFeatureNeighbors, featureNeighbors );

    for( unsigned int l = 0; l < numberOfMappedFeatureNeighbors; ++l ) {

      cdcl_feature_with_shape< 3 >::sptr  moving_feature = m_FeaturesMovingAll[0][featureNeighbors[l]];

      vnl_vector_fixed< double, 3 >  fixedNormalMappedVnl;

      // orienting the two coordinate systems (pick smallest difference between the directions)
      //vnl_vector_fixed< double, 3 >  moving_normal = moving_feature->directions_[0];
      //fixedNormalMappedVnl = fixedNormalsMappedVnl[0];
      //double smallest_difference = vcl_acos( dot_product( moving_normal, fixedNormalMappedVnl ) );
      vnl_vector_fixed< double, 3 >  moving_normal;
      double smallest_difference = vnl_math::pi;

      for( unsigned int d1 = 0; d1 < 1; ++d1 ) {
        for( unsigned int d2 = 0; d2 < 1; ++d2 ) {
      //for( unsigned int d1 = 0; d1 < fixedNormalsMappedVnl.size(); ++d1 ) {
      //  for( unsigned int d2 = 0; d2 < moving_feature->directions_.size(); ++d2 ) {
          double fixed_dot_moving = dot_product( moving_feature->directions_[d2], fixedNormalsMappedVnl[d1] );
          if( fixed_dot_moving > 1.0 ) fixed_dot_moving = 1.0;
          if( fixed_dot_moving < -1.0 ) fixed_dot_moving = -1.0;
          double difference = vcl_acos( fixed_dot_moving );

          // need to check if the direction is the tube tangent direction (if so, the difference might be flipped)
          if( ( moving_feature->shape_ == cdcl_feature_with_shape< 3 >::TUBE ) ^ ( fixed_feature->shape_ == cdcl_feature_with_shape< 3 >::TUBE ) ) {  // ^ is xor
            double flipped_difference = vcl_abs( difference - vnl_math::pi_over_2 );
            if( flipped_difference < difference ) difference = flipped_difference;
          }
          if( difference < smallest_difference ) {
            smallest_difference = difference;
            moving_normal = moving_feature->directions_[d2];
            fixedNormalMappedVnl = fixedNormalsMappedVnl[d1];
          }

          //vcl_cout << moving_feature->directions_[d1] << "  " << fixedNormalMappedVnl << "  " << difference << " ";
          //vcl_cout << d1 << "," << d2 << ": " << difference * 180.0 / vnl_math::pi << " ";

        }
      }
      //vcl_cout << endl;

      FeatureVectorType  closestFeature = kdTreeMovingFeatures->GetMeasurementVector( featureNeighbors[l] );
      vnl_vector_fixed< float, 3 >  candidate_distance = closestFeature.GetVnlVector() - fixedFeatureMapped.GetVnlVector();

      // compute histogram values for normals, binormals and distances; use descriptor match
      //double normal_difference = vcl_acos( dot_product( moving_normal, fixedNormalMappedVnl ) ) * 180.0 / vnl_math::pi; // angle might have been flipped (see above the case for tubes)
      double normal_difference = smallest_difference * 180.0 / vnl_math::pi;
      unsigned int normal_bin_number = vnl_math_rnd( normal_difference / 10.0 );

      //vcl_cout << " -> " << normal_difference << vcl_endl;

      if( candidate_distance.magnitude() < goodnessThreshold 
       && normal_difference < 20.0 ) ++goodFeatures;

      unsigned int distance_bin_number;
      if( candidate_distance.magnitude() > 3.0 ) {
        distance_bin_number = 3;
      }
      else {
        distance_bin_number = (unsigned int) ( vcl_floor( candidate_distance.magnitude() / 1.0 ) );
      }

      ++normals_distances[normal_bin_number][distance_bin_number];
    }

  }

  std::cout << "moving features: " << m_FeaturesMovingAll[0].size() << std::endl
            << "fixed features: " << m_FeaturesFixedAll[0].size() << std::endl
            << "good fixed mapped: " << goodFeatures << " percentage: " << 100.0 * double( goodFeatures ) / m_FeaturesFixedAll[0].size() << std::endl;

  // histograms of failures
  for( unsigned int a = 0; a < num_angles; ++a ) {
    vcl_cout << a*10 << "-" << (a+1)*10 << "\t";
  }
  vcl_cout << vcl_endl;

  vcl_cout << "nor_dist: " << vcl_endl;
  for( unsigned int d = 0; d < num_distances; ++d ) {
    for( unsigned int a = 0; a < num_angles; ++a ) {
      vcl_cout << vcl_setw( 4 ) << vcl_setprecision( 4 ) << vcl_setfill( '0' ) << normals_distances[a][d] * 100.0 / m_FeaturesFixedAll[0].size() << "\t";
    }
    vcl_cout << vcl_endl;
  }
 



  return 0;
}
