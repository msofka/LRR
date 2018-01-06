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

#include <cdcl/cdcl_feature.h>
#include <cdcl/cdcl_utils_VTK.h>

#include <vnl/algo/vnl_symmetric_eigensystem.h>
#include <vnl/algo/vnl_generalized_eigensystem.h>
#include <vnl/vnl_random.h>

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
// \date   September 2007


// write debugging image chips for each descriptor
#define OUTPUT_IMAGES 0


namespace {

// Kd-tree for descriptor vectors
//static const unsigned int descriptorSize = 771;
//typedef itk::Vector< float, descriptorSize >                      DescriptorVectorType;
typedef itk::Array< float >                      DescriptorVectorType;
typedef itk::Statistics::ListSample< DescriptorVectorType >       DescriptorSampleType;
typedef itk::Statistics::KdTreeGenerator< DescriptorSampleType >  DescriptorTreeGeneratorType;
typedef DescriptorTreeGeneratorType::KdTreeType                   DescriptorTreeType;


// Function object for sorting vector of descriptor indices based on distance
class DescriptorAtIndexIsGreater {
public:
  DescriptorAtIndexIsGreater( DescriptorVectorType &  FixedDescriptor,
                              DescriptorTreeType::Pointer const &  MovingDescriptorTree  )
    : m_FixedDescriptor( FixedDescriptor ),
      m_MovingDescriptorTree( MovingDescriptorTree ) {}
  bool operator()( int i, int j )
    {
      DescriptorVectorType  movingDescriptor_i = m_MovingDescriptorTree->GetMeasurementVector( i );
      vnl_vector< float >  distancei = m_FixedDescriptor - movingDescriptor_i;
      DescriptorVectorType  movingDescriptor_j = m_MovingDescriptorTree->GetMeasurementVector( j );
      vnl_vector< float >  distancej = m_FixedDescriptor - movingDescriptor_j;
      return distancei.magnitude() < distancej.magnitude();
    }
private:
  DescriptorVectorType const &  m_FixedDescriptor;
  DescriptorTreeType::Pointer  const &  m_MovingDescriptorTree;
};


//// Function object for sorting vector of descriptor indices based on ratio
//class RatioAtIndexIsGreater {
//public:
//  RatioAtIndexIsGreater( DescriptorVectorType &  FixedDescriptor,
//                         DescriptorTreeType::InstanceIdentifierVectorType const &  ClosestDescriptors,
//                         DescriptorTreeType::Pointer const &  MovingDescriptorTree  )
//    : m_FixedDescriptor( FixedDescriptor ),
//      m_ClosestDescriptors( ClosestDescriptors ),
//      m_MovingDescriptorTree( MovingDescriptorTree ) {}
//  bool operator()( int i, int j )
//    {
//      DescriptorVectorType  movingDescriptorClosest_i = m_MovingDescriptorTree->GetMeasurementVector( m_ClosestDescriptors[i] );
//      vnl_vector< float >  distanceClosest_i = m_FixedDescriptor - movingDescriptorClosest_i;
//
//      DescriptorVectorType  movingDescriptorSecondClosest_i = m_MovingDescriptorTree->GetMeasurementVector( m_ClosestDescriptors[i+1] );
//      vnl_vector< float >  distanceSecondClosest_i = m_FixedDescriptor - movingDescriptorSecondClosest_i;
//
//      double ratio_i = distanceClosest_i.magnitude() / distanceSecondClosest_i.magnitude();
//
//      DescriptorVectorType  movingDescriptorClosest_j = m_MovingDescriptorTree->GetMeasurementVector( m_ClosestDescriptors[j] );
//      vnl_vector< float >  distanceClosest_j = m_FixedDescriptor - movingDescriptorClosest_j;
//
//      DescriptorVectorType  movingDescriptorSecondClosest_j = m_MovingDescriptorTree->GetMeasurementVector( m_ClosestDescriptors[j+1] );
//      vnl_vector< float >  distanceSecondClosest_j = m_FixedDescriptor - movingDescriptorSecondClosest_j;
//
//      double ratio_j = distanceClosest_j.magnitude() / distanceSecondClosest_j.magnitude();
//
//      return ratio_i > ratio_j;
//    }
//private:
//  DescriptorVectorType const &  m_FixedDescriptor;
//  DescriptorTreeType::InstanceIdentifierVectorType const &  m_ClosestDescriptors;
//  DescriptorTreeType::Pointer  const &  m_MovingDescriptorTree;
//};


typedef std::pair< DescriptorTreeType::InstanceIdentifierVectorType, double >  distance_ratio_type;
bool DistanceRatioIsSmaller( distance_ratio_type const &  left,
                     distance_ratio_type const &  right )
{
  return left.second < right.second;
}


double RiemmanianDistance( const vnl_matrix_fixed< double, 3, 3 >  C1,
                           const vnl_matrix_fixed< double, 3, 3 >  C2 )
{
  if( C1( 0, 0 ) == 0.0 || C1( 1, 1 ) == 0.0 || C1( 2, 2 ) == 0.0
   || C2( 0, 0 ) == 0.0 || C2( 1, 1 ) == 0.0 || C2( 2, 2 ) == 0.0 ) return 0.0;
  vnl_generalized_eigensystem  geig_key( C1, C2 );
  //vcl_cout << "d0d1d2: " << geig_key.D( 0, 0 ) << " " << geig_key.D( 1, 1 ) << " " << geig_key.D( 2, 2 ) << vcl_endl;
  if( geig_key.D( 0, 0 ) == 0 || geig_key.D( 1, 1 ) == 0 || geig_key.D( 2, 2 ) == 0 ) return 0.0;
  double D0 = vcl_log( geig_key.D( 0, 0 ) );
  double D1 = vcl_log( geig_key.D( 1, 1 ) );
  double D2 = vcl_log( geig_key.D( 2, 2 ) );
  //vcl_cout << "log: " << D0 << " " << D1 << " " << D2 << vcl_endl;
  //vcl_cout << vcl_sqrt( D0*D0 + D1*D1 + D2*D2 ) << vcl_endl;
  return vcl_sqrt( D0*D0 + D1*D1 + D2*D2 );
}


double DescriptorRiemmanianDistance( vnl_vector< float > const &  FixedDescriptor,
                                     vnl_vector< float > const &  MovingDescriptor )
{
  assert( FixedDescriptor.size() == MovingDescriptor.size() );

  if( FixedDescriptor.size() % 6 != 0 ) return 0.0;
  
  double totalDistance = 0.0;
  for( unsigned int i = 0; i < FixedDescriptor.size()-6; i += 6 ) {
    vnl_matrix_fixed< double, 3, 3 >  fixedC;
    fixedC( 0, 0 ) = FixedDescriptor[ i   ];
    fixedC( 0, 1 ) = FixedDescriptor[ i+1 ];
    fixedC( 0, 2 ) = FixedDescriptor[ i+2 ];

    fixedC( 1, 1 ) = FixedDescriptor[ i+3 ];
    fixedC( 1, 2 ) = FixedDescriptor[ i+4 ];

    fixedC( 2, 2 ) = FixedDescriptor[ i+5 ];

    fixedC( 1, 0 ) = fixedC( 0, 1 );
    fixedC( 2, 0 ) = fixedC( 0, 2 );
    fixedC( 2, 1 ) = fixedC( 1, 2 );

    vnl_matrix_fixed< double, 3, 3 >  movingC;
    movingC( 0, 0 ) = MovingDescriptor[ i   ];
    movingC( 0, 1 ) = MovingDescriptor[ i+1 ];
    movingC( 0, 2 ) = MovingDescriptor[ i+2 ];

    movingC( 1, 1 ) = MovingDescriptor[ i+3 ];
    movingC( 1, 2 ) = MovingDescriptor[ i+4 ];

    movingC( 2, 2 ) = MovingDescriptor[ i+5 ];

    movingC( 1, 0 ) = movingC( 0, 1 );
    movingC( 2, 0 ) = movingC( 0, 2 );
    movingC( 2, 1 ) = movingC( 1, 2 );

    //vcl_cout << "first: " << RiemmanianDistance( fixedC, movingC ) << " second: " << RiemmanianDistance( fixedC, movingC ) << vcl_endl;
    totalDistance += RiemmanianDistance( fixedC, movingC );
  }

  return totalDistance;
}


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

// find moving descriptors closest to a fixed descriptor using Riemmanian framework
void SearchRiemmanianNeighbors( DescriptorVectorType const &  fixedDescriptor,
                                DescriptorSampleType::Pointer const &  samples,
                                unsigned int k,
                                DescriptorTreeType::InstanceIdentifierVectorType &  result )
{
  //vcl_cout << "FIXED in riemmanian: " << fixedDescriptor << vcl_endl << vcl_endl;
  vcl_vector< vcl_pair< DescriptorVectorType, double > >  distances;

  vcl_vector< DistanceID >  distIDs;
  for( DescriptorTreeType::ConstIterator it = samples->Begin(); it != samples->End(); ++it ) {
    DistanceID  distID;
    distID.vector_ = it.GetMeasurementVector();
    distID.id_ = it.GetInstanceIdentifier();
    DescriptorVectorType vect = samples->GetMeasurementVector( distID.id_ );
    //if( vect != distID.vector_ ) {
    //  std::cout << "VECTORS DIFFERENT: " << std::endl;
    //  std::cout << vect << std::endl;
    //  std::cout << distID.vector_ << std::endl;
    //}
    distID.distance_ = DescriptorRiemmanianDistance( fixedDescriptor, distID.vector_ );
    distIDs.push_back( distID );
  }


  vcl_sort( distIDs.begin(), distIDs.end(), DistanceForMeasurementIsGreater );

  for( vcl_vector< DistanceID >::size_type  i = 0; i < k; ++i ) {
    DescriptorVectorType vect = samples->GetMeasurementVector( distIDs[i].id_ );
    //if( vect != distIDs[i].vector_ ) {
    //  std::cout << "VECTORS DIFFERENT: " << std::endl;
    //  std::cout << vect << std::endl;
    //  std::cout << distIDs[i].vector_ << std::endl;
    //}
    //vcl_cout << "sort: " << i << ": " << DescriptorRiemmanianDistance( fixedDescriptor, distIDs[i].vector_ ) << "  " << distIDs[i].id_ << " basedonid: " << DescriptorRiemmanianDistance( fixedDescriptor, vect ) << " and again: " << DescriptorRiemmanianDistance( fixedDescriptor, vect ) << vcl_endl;
    //vcl_cout << "sort: " << i << ": " << distIDs[i].distance_ << "  " << distIDs[i].id_ << " basedonid: " << DescriptorRiemmanianDistance( fixedDescriptor, vect ) << vcl_endl;
  }


  //for( DescriptorTreeType::InstanceIdentifier i = 0; i < samples->Size(); ++i ) {
  //  DescriptorVectorType  movingDescriptor = samples->GetMeasurementVector( i );
  //  double distance = DescriptorRiemmanianDistance( fixedDescriptor, movingDescriptor );
  //  distances.push_back( vcl_pair< DescriptorVectorType, double >( movingDescriptor, distance ) );
  //}

  //vcl_sort( distances.begin(), distances.end(), DistanceForMeasurementIsGreater );

  result.clear();
  for( vcl_vector< DistanceID >::size_type  i = 0; i < k; ++i ) {
    //if( i == 0 ) vcl_cout << "MOVING in riemmanian " << distIDs[i].vector_ << vcl_endl << " dist: " << distIDs[i].distance_ << vcl_endl;
    DescriptorTreeType::InstanceIdentifier  id = distIDs[i].id_;
    result.push_back( id );
  }
}


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
  // e.g. indexing_shape_context3dITK.exe 241/4804/000002croppeddesc0.vtk 241/4802/000002croppeddesc0.vtk -trans 241/4802_000002cropped-4804_000002cropped.vtk -mov 241/4804/000002cropped -fix 241/4802/000002cropped
  vul_arg< const char* >        movingDescriptorFile     ( 0, "VTK moving descriptor file (of a synthetically deformed image)" );
  vul_arg< const char* >        fixedDescriptorFile      ( 0, "VTK fixed descriptor file" );
  vul_arg< const char* >        transformFile            ( "-trans", "ITK transform file (transform used to warp points)", 0 );
  vul_arg< const char* >        movingFeatureFilePrefix  ( "-fea", "VTK moving feature file prefix (of a deformed image) to test keypoint repeatability" );
  vul_arg< unsigned int >       levelInd                 ( "-ind", "level index of the feature file (-fea)" );
  vul_arg< const char* >        movingImageFile          ( "-mov", "moving image file (of a deformed image)" );
  vul_arg< const char* >        fixedImageFile           ( "-fix", "fixed image file" );

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




  // Kd-tree for moving keypoint locations
  typedef itk::Vector< float, 3 >                                 KeypointVectorType;
  typedef itk::Statistics::ListSample< KeypointVectorType >       KeypointSampleType;
  typedef itk::Statistics::KdTreeGenerator< KeypointSampleType >  KeypointTreeGeneratorType;
  typedef KeypointTreeGeneratorType::KdTreeType                   KeypointTreeType;



  // setup indexing

  DescriptorSampleType::Pointer  sampleMovingDescriptors = DescriptorSampleType::New();

  KeypointSampleType::Pointer    sampleMovingKeypoints = KeypointSampleType::New();

  vcl_vector< cdcl_keypoint< 3 >::sptr >                movingKeypoints;
  // read moving descriptors
  {
    vcl_vector< vnl_vector< float > >                    movingDescriptors;
    cdcl_read_keypoint_descriptors_VTK( movingKeypoints, movingDescriptors, movingDescriptorFile() );


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

  // moving keypoint descriptor vectors
  DescriptorTreeGeneratorType::Pointer  descriptorTreeGenerator = DescriptorTreeGeneratorType::New();

  descriptorTreeGenerator->SetSample( sampleMovingDescriptors );
  descriptorTreeGenerator->SetBucketSize( 16 );
  descriptorTreeGenerator->Update();

  DescriptorTreeType::Pointer  kdTreeMovingDescriptors = descriptorTreeGenerator->GetOutput();

  vcl_cout << "Kd tree size: " << kdTreeMovingDescriptors->Size() << vcl_endl;


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
  cdcl_read_keypoint_descriptors_VTK( fixedKeypoints, fixedDescriptors, fixedDescriptorFile() );


  KeypointSampleType::Pointer    sampleFixedKeypoints = KeypointSampleType::New();

  // read fixed keypoints
  {
    // moving keypoint locations
    //
    KeypointVectorType  fi_k;
    for( unsigned int n = 0 ; n < fixedKeypoints.size(); ++n ) {
      fi_k[0] = fixedKeypoints[n]->location_[0];
      fi_k[1] = fixedKeypoints[n]->location_[1];
      fi_k[2] = fixedKeypoints[n]->location_[2];
      //mv.Set_vnl_vector( m_FixedKeypoints[n]->location_.as_vector() );
      sampleFixedKeypoints->PushBack( fi_k );
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


  unsigned int numberOfKeypointsTested = 1000;

  srand( (unsigned)time( NULL ) );
  vnl_random rng;
  rng.reseed(333248);
  // Usually, you will want to generate a number in a specific range
  int RANGE_MIN = 0;
  int RANGE_MAX = fixedDescriptors.size()-1;


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

  std::cout << "Keypoint Distance\tDescriptor Difference\t" << std::endl;
  unsigned int goodAtFirst = 0;
  unsigned int goodAtAll = 0;
  unsigned int goodIn10 = 0;
  unsigned int goodKeypoints = 0;
  unsigned int atOrOutsideBoundary = 0;

  unsigned int numberOfDescriptorNeighbors = std::min( 100u, sampleMovingDescriptors->Size() );
  float goodnessThreshold = 9.0f;
  float keypointGoodnessThreshold = 9.0f;

  for( unsigned int t = 0; t < numberOfKeypointsTested; ++t ) {

    KeypointVectorType  queryPoint;
    for( unsigned int d = 0;  d < 3; ++d ) {
      // generate a number in a specific range
      // subtract radius of the shape context
      double context_radius = 30;
      double RANGE_MIN = minFixedPoint[d] + context_radius - 5.0; // 5 mm tolerance
      double RANGE_MAX = maxFixedPoint[d] - context_radius - 5.0;

      //double randNum = ((double) rand() / (double) RAND_MAX) * (RANGE_MAX - RANGE_MIN) + RANGE_MIN;
      double randNum = rng.drand32( RANGE_MIN, RANGE_MAX );

      queryPoint[d] = randNum;
    }

    std::cout << "Clicked: " << queryPoint << " ";

    // get the closest keypoint to the clicked point
    const unsigned int numberOfKeypointNeighbors = 1;
    KeypointTreeType::InstanceIdentifierVectorType  fixedKeypointNeighbors;
    kdTreeFixedKeypoints->Search( queryPoint, numberOfKeypointNeighbors, fixedKeypointNeighbors );

    for( unsigned int l = 0; l < numberOfKeypointNeighbors; ++l ) {

      std::cout << "n: " << l << "  " << fixedKeypointNeighbors[l] << std::endl;
      // Descriptors closest to the moving descriptor
      DescriptorTreeType::InstanceIdentifierVectorType  descriptorNeighbors;
      
      vnl_vector< float >  fixed_descriptor = fixedDescriptors[fixedKeypointNeighbors[l]];;
      DescriptorVectorType  fixedDescriptor( fixed_descriptor.size() );
      //fixedDescriptor.SetData( fixed_descriptor.data_block() );
        for( unsigned int d = 0; d < fixed_descriptor.size(); ++d ) {
          //mv.SetSize( movingDescriptors[n].size() );
          //mv.SetData( movingDescriptors[n].data_block() );
          fixedDescriptor.SetElement( d, fixed_descriptor[d] );
        }

  #define USE_RIEMMANIAN 0

  #if USE_RIEMMANIAN
      SearchRiemmanianNeighbors( fixedDescriptor, sampleMovingDescriptors, numberOfDescriptorNeighbors, descriptorNeighbors );
  #else
      // find descriptor in the moving image closest to the descriptor in the fixed image
      kdTreeMovingDescriptors->Search( fixedDescriptor, numberOfDescriptorNeighbors, descriptorNeighbors ) ;

      // sort the moving descriptors based on their distance from the fixed descriptor
      DescriptorAtIndexIsGreater  descriptorCompare( fixedDescriptor, kdTreeMovingDescriptors );

  #define SORT_ON_DISTANCE 1

  #if SORT_ON_DISTANCE
      std::sort( descriptorNeighbors.begin(), descriptorNeighbors.end(), descriptorCompare );

      //for( std::vector< DescriptorTreeType::InstanceIdentifierVectorType >::size_type  s = 0; s < descriptorNeighbors.size(); ++s ) {
      //  vnl_vector< float >  distance = fixedDescriptor - sampleMovingDescriptors->GetMeasurementVector( descriptorNeighbors[s] );
      //  std::cout << "Distance: " << distance.magnitude() << std::endl; 
      //}
  #else // SORT_ON_RATIO
      // this is not the fastest implementation, since we have to sort twice and the distances are also computed twice
      //
      // sort based on distances first
      DescriptorAtIndexIsGreater  descriptorClosest( fixedDescriptor, kdTreeMovingDescriptors );
      DescriptorTreeType::InstanceIdentifierVectorType  closestDescriptors;
      closestDescriptors.resize( descriptorNeighbors.size() );
      std::copy( descriptorNeighbors.begin(), descriptorNeighbors.end(), closestDescriptors.begin() );
      std::sort( closestDescriptors.begin(), closestDescriptors.end(), descriptorClosest );

      // compute ratios and augment the original with the ratios
      typedef std::pair< DescriptorTreeType::InstanceIdentifier, float >  distance_ratio_type;
      std::vector< distance_ratio_type >  closestOnRatioDescriptors;
      closestOnRatioDescriptors.resize( closestDescriptors.size()-1 );
      for( std::vector< distance_ratio_type >::size_type  s = 0; s < closestDescriptors.size()-1; ++s ) {
        vnl_vector< float >  distance1 = fixedDescriptor - sampleMovingDescriptors->GetMeasurementVector( closestDescriptors[s] );
        vnl_vector< float >  distance2 = fixedDescriptor - sampleMovingDescriptors->GetMeasurementVector( closestDescriptors[s+1] );
        float ratio = distance1.magnitude() / distance2.magnitude();
        //std::cout << "DISTANCE: " << distance1.magnitude() << std::endl;
        closestOnRatioDescriptors[s] = distance_ratio_type( closestDescriptors[s], ratio );
      }
      // sort on ratios
      std::sort( closestOnRatioDescriptors.begin(), closestOnRatioDescriptors.end(), DistanceRatioIsSmaller );  // don't sort the last one (there is no point further away that that one)

      // copy into the original vector
      descriptorNeighbors.clear();
      descriptorNeighbors.resize( closestOnRatioDescriptors.size() );
      for( std::vector< DescriptorTreeType::InstanceIdentifierVectorType >::size_type  s = 0; s < closestOnRatioDescriptors.size(); ++s ) {
        descriptorNeighbors[s] = closestOnRatioDescriptors[s].first;
        //std::cout << "RATIO: " << closestOnRatioDescriptors[s].first << "   " << closestOnRatioDescriptors[s].second << std::endl;
      }
  #endif

  #endif

      vnl_vector_fixed< double, 3 >  fixedNormalMappedVnl;
      vnl_vector_fixed< double, 3 >  fixedBinormalMappedVnl;
      
      // moving image was obtained by deforming the fixed image, so deform the fixed keypoint with the same transform
      // and compare the two points
      //BSplineTransformType::OutputPointType  movingMapped = bsplineTransform->TransformPoint( moving_keypoint->location_.data_block() );
      cdcl_keypoint< 3 >::sptr  fixed_keypoint  = fixedKeypoints[fixedKeypointNeighbors[l]];
      vnl_vector< double >  fixed_normal   = fixed_keypoint->normal_;
      vnl_vector< double >  fixed_binormal = fixed_keypoint->binormal_;
      vnl_vector_fixed< double, 3 >  fixedMappedVnl;
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
          // It can happen that the keypoint location is right at the boundary of the image. Do nothing in this case.
          if( !(deformationField->GetLargestPossibleRegion().IsInside( fixedIndex ) ) ) {
            std::cout << "Warning: Keypoint index is " << fixedIndex << " but the deformation field size is: " << deformationField->GetLargestPossibleRegion() << std::endl;
            ++atOrOutsideBoundary;
            continue;
          }
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
          if( !(deformationField->GetLargestPossibleRegion().IsInside( fixedIndexNor ) ) ) {
            std::cout << "Warning: Keypoint normal ahead index is " << fixedIndexNor << " but the deformation field size is: " << deformationField->GetLargestPossibleRegion() << std::endl;
            ++atOrOutsideBoundary;
            continue;
          }
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
          if( !(deformationField->GetLargestPossibleRegion().IsInside( fixedIndexBin ) ) ) {
            std::cout << "Warning: Keypoint binormal ahead index is " << fixedIndexBin << " but the deformation field size is: " << deformationField->GetLargestPossibleRegion() << std::endl;
            ++atOrOutsideBoundary;
            continue;
          }
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

      // after this loop moving_descriptor and fixed_descriptor will contain the closest (best) match
      // this loop is only for printout of the other matches
      vnl_vector< float >  moving_descriptor;
      vnl_vector_fixed< double, 3 >  distance;
      double normal_difference = 0.0;
      double binormal_difference = 0.0;
      bool goodFound = false;
      bool goodFoundIn10 = false;
      cdcl_keypoint< 3 >::sptr  moving_keypoint;
      for( unsigned int l = 0; l < numberOfDescriptorNeighbors; ++l ) {

  #if USE_RIEMMANIAN
        //vcl_cout << "neigh: " << descriptorNeighbors[numberOfDescriptorNeighbors-l-1] << vcl_endl;
        moving_keypoint = movingKeypoints[descriptorNeighbors[numberOfDescriptorNeighbors-l-1]];
        moving_descriptor = sampleMovingDescriptors->GetMeasurementVector( descriptorNeighbors[numberOfDescriptorNeighbors-l-1] );
    
        //vcl_cout << "FIXED in compare: " << fixedDescriptor << vcl_endl << vcl_endl;
        //vcl_cout << "MOVING in compare: " << moving_descriptor << vcl_endl << vcl_endl;

        //vcl_cout << "afterprintout: " << DescriptorRiemmanianDistance( fixedDescriptor, moving_descriptor );

  #else
        moving_keypoint = movingKeypoints[descriptorNeighbors[numberOfDescriptorNeighbors-l-1]];
        moving_descriptor = kdTreeMovingDescriptors->GetMeasurementVector( descriptorNeighbors[numberOfDescriptorNeighbors-l-1] );
  #endif
        distance = fixedMappedVnl - moving_keypoint->location_;

        vnl_vector< float >  difference = moving_descriptor - fixed_descriptor;

        vnl_vector< double >  moving_normal   = moving_keypoint->normal_;
        vnl_vector< double >  moving_binormal = moving_keypoint->binormal_;

        // compute histogram values for normals, binormals and distances; use the match
        normal_difference = vcl_acos( dot_product( moving_normal, fixedNormalMappedVnl ) ) * 180.0 / vnl_math::pi;
        binormal_difference = vcl_acos( dot_product( moving_binormal, fixedBinormalMappedVnl ) ) * 180.0 / vnl_math::pi;

        std::cout << numberOfDescriptorNeighbors-l-1 << ":\t";
        std::cout << distance.magnitude() << "\t";
        std::cout << difference.magnitude() << "\t";
        std::cout << DescriptorRiemmanianDistance( fixedDescriptor, moving_descriptor ) << std::endl;
        bool satisfies_threshold = distance.magnitude() < goodnessThreshold
                                && normal_difference < 20.0
                                && binormal_difference < 20.0;
        if( satisfies_threshold ) goodFound = true;
        if( satisfies_threshold && ( l > (numberOfDescriptorNeighbors-11) ) ) goodFoundIn10 = true;
      }

      // now distance variable contains distance to the closest descriptor
      if( distance.magnitude() < goodnessThreshold
          && normal_difference < 20.0
          && binormal_difference < 20.0 ) ++goodAtFirst;
      if( goodFound ) ++goodAtAll;
      if( goodFoundIn10 ) ++goodIn10;


      unsigned int normal_bin_number = vnl_math_rnd( normal_difference / 10.0 );
      unsigned int binormal_bin_number = vnl_math_rnd( binormal_difference / 10.0 );

      unsigned int distance_bin_number;
      if( distance.magnitude() > 9.0 ) {
        distance_bin_number = 3;
      }
      else {
        distance_bin_number = (unsigned int) ( vcl_floor( distance.magnitude() / 3.0 ) );
      }

      ++normals_distances[normal_bin_number][distance_bin_number];
      ++binormals_distances[binormal_bin_number][distance_bin_number];

      if( volumesRead ) {
        vcl_ostringstream  movingImageFile;
        movingImageFile << "moving_image" << vcl_setw( 6 ) << vcl_setfill( '0' ) << t << ".png";      

        WriteROIImage<MovingImageType>( movingImage, moving_keypoint, movingImageFile.str() );

        vcl_ostringstream  fixedImageFile;
        fixedImageFile << "fixed_image" << vcl_setw( 6 ) << vcl_setfill( '0' ) << t << ".png";      

        WriteROIImage<MovingImageType>( fixedImage, fixed_keypoint, fixedImageFile.str() );
      }


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

      std::cout << "Fixed: nor: " << fixed_normal << " binor: " << fixed_binormal << std::endl;

      unsigned int numberOfMappedKeypointNeighbors = 1;
      KeypointTreeType::InstanceIdentifierVectorType  mappedKeypointNeighbors;
      kdTreeMovingKeypoints->Search( fixedKeypointMapped, numberOfMappedKeypointNeighbors, mappedKeypointNeighbors );
      for( unsigned int m = 0; m < numberOfMappedKeypointNeighbors; ++m ) {
        KeypointVectorType  closestKeypoint = kdTreeMovingKeypoints->GetMeasurementVector( mappedKeypointNeighbors[m] );
        vnl_vector_fixed< float, 3 >  candidate_distance = closestKeypoint.GetVnlVector() - fixedKeypointMapped.GetVnlVector();
        DescriptorVectorType  closestDescriptor = kdTreeMovingDescriptors->GetMeasurementVector( mappedKeypointNeighbors[m] );
        vnl_vector< float >  candidate_difference = closestDescriptor - fixedDescriptor;

        std::cout << "Keypoint " << m << ": " << candidate_distance.magnitude() << "\t" << candidate_difference.magnitude() << " ";

        cdcl_keypoint<3>::sptr  moving_keypoint_closest = movingKeypoints[mappedKeypointNeighbors[m]];

        vnl_vector< double >  moving_normal   = moving_keypoint_closest->normal_;
        vnl_vector< double >  moving_binormal = moving_keypoint_closest->binormal_;

        if( volumesRead ) {
          vcl_ostringstream  movingImageFile;
          movingImageFile << "moving_image_closest" << vcl_setw( 6 ) << vcl_setfill( '0' ) << t << ".png";      

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

        std::cout << "Fixed: nor: " << vcl_setprecision( 3 ) << fixed_normal << " binor: " << fixed_binormal << std::endl;


        // compute histogram values for normals, binormals and distances; use closest keypoint
        double normal_difference = vcl_acos( dot_product( moving_normal, fixedNormalMappedVnl ) ) * 180.0 / vnl_math::pi;
        unsigned int normal_bin_number = vcl_floor( normal_difference / 10.0 );

        double binormal_difference = vcl_acos( dot_product( moving_binormal, fixedBinormalMappedVnl ) ) * 180.0 / vnl_math::pi;
        unsigned int binormal_bin_number = vcl_floor( binormal_difference / 10.0 );

        if( candidate_distance.magnitude() < keypointGoodnessThreshold
          && normal_difference < 20
          && binormal_difference < 20 ) ++goodKeypoints;

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
        fileNameMoving << "desc_candidate" << vcl_setw( 6 ) << vcl_setfill( '0' ) << k << ".png";      
        //vcl_cout << "Saving " << fileNameMoving.str() << vcl_endl;

        writer->SetFileName( fileNameMoving.str().c_str() );
        writer->Update();
  #endif
      }
      std::cout << std::endl;


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
      fileNameMoving << "desc_moving" << vcl_setw( 6 ) << vcl_setfill( '0' ) << k << ".png";
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
      fileNameFixed << "desc_fixed" << vcl_setw( 6 ) << vcl_setfill( '0' ) << k << ".png";
      //vcl_cout << "Saving " << fileNameFixed.str() << vcl_endl;

      writer->SetFileName( fileNameFixed.str().c_str() );
      writer->Update();
  #endif

    }
  }
  std::cout << "moving keypoints: " << movingKeypoints.size() << std::endl
            << "fixed keypoints: " << fixedKeypoints.size() << std::endl
            << "keypoints tested: " << numberOfKeypointsTested << std::endl
            << "descriptor neighbors for each: " << numberOfDescriptorNeighbors << std::endl
            << "good matches on the first hit: " << goodAtFirst << std::endl
            << "good matches in the first 10: " << goodIn10 << std::endl
            << "good matches in " << numberOfDescriptorNeighbors << " of all: " << goodAtAll << std::endl
            << "good keypoints (tolerance " << keypointGoodnessThreshold << ", 20, 20): " << goodKeypoints << std::endl
            << "at or outside boundary: " << atOrOutsideBoundary << std::endl;

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


  return 0;
}
