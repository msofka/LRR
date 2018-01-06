#include <vcl_utility.h>
#include <vcl_algorithm.h>
#include <vcl_iomanip.h>
#include <vcl_iostream.h>
#include <vcl_iterator.h>
#include <vcl_cstdlib.h>
#include <vcl_set.h>

#include <vul/vul_arg.h>

#include <vul/vul_timer.h>
#include <vnl/algo/vnl_svd.h>
#include <vnl/vnl_inverse.h>
#include <vnl/vnl_transpose.h>
#include <vnl/vnl_random.h>
#include <vnl/vnl_cross.h>

#include <rsdl/rsdl_dist.h>
#include <vcl_memory.h>

#include <cdcl/cdcl_utils.h>
#include <cdcl/cdcl_trans_rigid3d.h>
#include <cdcl/cdcl_trans_affine.h>
#include <cdcl/cdcl_feature_with_shape.h>

#include "itkBoundingBox.h"
#include "itkDanielssonDistanceMapImageFilter.h"

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "vtkSmartPointer.h"
#include "vtkXMLPolyDataReader.h"

#include <rsdl/rsdl_kd_tree.h>

#include <cdcl/io/vtkPolyDataToFeaturesWithShapeFilter.h>


#undef VTK_FOUND
// VTK is used to dump data for analysis in external programs.
#ifdef VTK_FOUND
#include <cdcl/cdcl_utils_VTK.h>
#endif

//:
// \file
// \brief  Compute and test distance map for rapid matching.
// \author Michal Sofka
// \date   May 2008



// Compute distance map for fast matching.
void compute_distance_map( vcl_vector< vcl_vector< cdcl_feature_with_shape< 3 >::sptr > > const &  fixed_,
                           itk::Image< unsigned int, 3 >::Pointer &  voronoi_map_ )
{

  const unsigned int res_ = 0;

  //  From ITK Examples.
  //  Since the output will contain distances measured in pixels, the
  //  pixel type should be able to represent at least the width of the image,
  //  or said in $N-D$ terms, the maximum extension along all the dimensions.
  //  The input and output image types are now defined using their respective
  //  pixel type and dimension.

  typedef  unsigned int  OutputPixelType;
  typedef itk::Image< OutputPixelType,  3 >   OutputImageType;
  //typedef itk::RGBPixel< unsigned char > RGBPixelType;
  //typedef itk::Image<RGBPixelType, 2>                           OutputImageType;

  //typedef  cdcl_feature_with_shape< 3 >::sptr  InputPixelType;
  typedef  unsigned int  InputPixelType;
  typedef itk::Image< InputPixelType,  3 >   InputImageType;

  // Compute the bounding box of all points
  //
  typedef  itk::BoundingBox< unsigned long, 3, double >  BoundingBoxType;
  typedef BoundingBoxType::PointType  PointType;

  BoundingBoxType::Pointer  fixedBoundingBox = BoundingBoxType::New();
  BoundingBoxType::PointsContainerPointer  fixedPoints = BoundingBoxType::PointsContainer::New();
  
  for( vcl_vector< cdcl_feature_with_shape< 3 >::sptr >::size_type  i = 0; i < fixed_[res_].size(); ++i ) {
    PointType  point;
    point[0] = fixed_[res_][i]->location_[0];
    point[1] = fixed_[res_][i]->location_[1];
    point[2] = fixed_[res_][i]->location_[2];
    fixedPoints->push_back( point );
  }
  fixedBoundingBox->SetPoints( fixedPoints );

  PointType  minFixedPoint = fixedBoundingBox->GetMinimum();
  PointType  maxFixedPoint = fixedBoundingBox->GetMaximum();


  // Allocate the input image, which will store feature indices
  //
  InputImageType::SpacingType  inputSpacing;
  inputSpacing.Fill( 1.0 );

  InputImageType::SizeType  inputSize;
  for( unsigned int d = 0; d < 3; ++d ) inputSize[d] = vcl_ceil( maxFixedPoint[d] - minFixedPoint[d] ) / inputSpacing[d];

  InputImageType::IndexType  inputStart;
  inputStart.Fill( 0 );

  InputImageType::RegionType  inputRegion;
  inputRegion.SetIndex( inputStart );
  inputRegion.SetSize( inputSize );

  InputImageType::Pointer  featureImage = InputImageType::New();  
  featureImage->SetOrigin( minFixedPoint );
  featureImage->SetRegions( inputRegion );
  featureImage->Allocate();
  featureImage->FillBuffer( 0 );


  // Fill the input image with feature locations
  //
  for( vcl_vector< cdcl_feature_with_shape< 3 >::sptr >::size_type  i = 0; i < fixed_[res_].size(); ++i ) {
    cdcl_feature_with_shape< 3 >::sptr  curr_feature = fixed_[res_][i];

    itk::Point< double, 3 >  location;
    location[0] = curr_feature->location_[0];
    location[1] = curr_feature->location_[1];
    location[2] = curr_feature->location_[2];

    InputImageType::IndexType  index;
    featureImage->TransformPhysicalPointToIndex( location, index );

    //featureImage->SetPixel( index, curr_feature );
    featureImage->SetPixel( index, i );
  }



  typedef itk::DanielssonDistanceMapImageFilter<
               InputImageType, OutputImageType >  FilterType;
  FilterType::Pointer filter = FilterType::New();

  filter->SetInput( featureImage );

  filter->Update();

  //// note we needed to modify the original distance map since we want the Voronoi Map to be the same type as the input
  //InputImageType::Pointer  voronoiMap = filter->GetVoronoiMap();

  voronoi_map_ = filter->GetVoronoiMap();
}


void ReadFeatures( vcl_string  image_dir,
                   vcl_vector< vcl_vector< cdcl_feature_with_shape< 3 >::sptr > > &  feature_vect )
{
  // Read features
  //
  vtkSmartPointer< vtkXMLPolyDataReader >  featureReader = vtkSmartPointer< vtkXMLPolyDataReader >::New();
  vcl_string  filename = image_dir + "_00.vtk";
  featureReader->SetFileName( filename.c_str() );
  featureReader->Update();  // need this update (bug in the converting filter)

  // create the converting filter
  vtkSmartPointer< vtkPolyDataToFeaturesWithShapeFilter >  polyDataToFeaturesFilter = vtkSmartPointer< vtkPolyDataToFeaturesWithShapeFilter >::New();
  polyDataToFeaturesFilter->SetInput( featureReader->GetOutput() );
  polyDataToFeaturesFilter->Update();
  vtkSmartPointer< vtkFeatureWithShapeAttributeSet >  featureAttributeSetVTK = polyDataToFeaturesFilter->GetOutput();

  typedef vtkFeatureWithShapeAttributeSet::FeatureSetType  FeatureSetType;
  FeatureSetType &  features = featureAttributeSetVTK->GetPoints();
  feature_vect.push_back( features );
}


void CompareMapTree( itk::Image< unsigned int, 3 >::Pointer &  voronoi_map_,
                     vcl_vector< vcl_vector< cdcl_feature_with_shape< 3 >::sptr > > &  fixed_,
                     vcl_vector< vcl_vector< cdcl_feature_with_shape< 3 >::sptr > > &  moving_ )
{
  typedef cdcl_feature_with_shape< 3 >::sptr  feature_sptr_type;
  typedef itk::Image< unsigned int, 3 >  VoronoiMapType;

  // Build kd-trees
  //
  // set the coarsest resolution
  unsigned int res_ = 0;

  const unsigned int dim = 3;
  const unsigned int dof = 12;

  // Kd-tree formed from fixed feature points, vector of fine-to-coarse resolutions.
  vcl_vector< vbl_smart_ptr< rsdl_kd_tree > >  kd_tree_;

  // build a set of kd_tree's, one for each resolution level of the fixed image
  for( vcl_vector< vcl_vector< feature_sptr_type > >::size_type n = 0; n < fixed_.size(); ++n ) {
    // build kd_tree from fixed points
    vcl_vector< rsdl_point > points;
    rsdl_point  pt( dim, 0 );
    for( vcl_vector< feature_sptr_type >::size_type j = 0; j < fixed_[n].size(); ++j ) {
      pt.set_cartesian( fixed_[n][j]->location_ );
      points.push_back( pt );
    }

    vcl_cout << "kdtree size: " << points.size() << vcl_endl;
    kd_tree_.push_back( new rsdl_kd_tree( points ) );
  }

  
  VoronoiMapType::SizeType  voronoiSize = voronoi_map_->GetLargestPossibleRegion().GetSize();

  // initialize transform to identity
  vnl_matrix_fixed< double, dim, dim >  A( 0.0 );
  A( 0, 0 ) = 1.0;
  A( 1, 1 ) = 1.0;
  A( 2, 2 ) = 1.0;

  vnl_vector_fixed< double, dim >  t( 0.0 );
  vnl_vector_fixed< double, dim >  center_moving( 0.0 );

  cdcl_trans_affine<dim,dof>::sptr  trans_ = new cdcl_trans_affine<dim, dof>( A, t, center_moving );

  // query point for kd-tree
  rsdl_point  query_pt( dim, 0 );
  // use n nearest neighbors
  int n = 1;
  vcl_vector< rsdl_point >  points( n );
  vcl_vector< int >  indices( n );

  // we will pick number_points_ number of moving points randomly
  vcl_random_shuffle( moving_[res_].begin(), moving_[res_].end() );

  // find matches based on Euclidean distance, and store them
  //
  unsigned int number_points_ = 1000000;
  number_points_ = vcl_min( number_points_, (unsigned int) moving_[res_].size() );
  vcl_cout << "Using: " << number_points_ << " moving points." << vcl_endl;

  vul_timer  timer_tree;
  timer_tree.mark();
  // go through all moving points
//  for( typename vcl_vector< feature_sptr_type >::size_type  i = 0; i < moving_[res_].size(); ++i ) {
  for( vcl_vector< feature_sptr_type >::size_type  i = 0; i < number_points_; ++i ) {
    // moving point p
    feature_sptr_type  p = moving_[res_][i];
    vnl_vector_fixed< double, dim >  mapped_moving = trans_->map_loc( p->location_ );

    query_pt.set_cartesian( mapped_moving );
    // query kd-tree and find set of nearest points as potential matches
    kd_tree_[res_]->n_nearest( query_pt, n, points, indices );

    feature_sptr_type  q = fixed_[res_][indices[0]];

    //vcl_cout << "TREE: " << q->location_ << "  " << vcl_endl; 
    
  }

  long time_tree = timer_tree.real();

  vul_timer  timer_map;
  timer_map.mark();
  for( vcl_vector< feature_sptr_type >::size_type  i = 0; i < number_points_; ++i ) {
    // moving point p
    feature_sptr_type  p = moving_[res_][i];
    vnl_vector_fixed< double, dim >  mapped_moving = trans_->map_loc( p->location_ );

    itk::Point< double, 3 >  mapped_moving_itk;
    mapped_moving_itk[0] = mapped_moving[0];
    mapped_moving_itk[1] = mapped_moving[1];
    mapped_moving_itk[2] = mapped_moving[2];
    VoronoiMapType::IndexType  mapped_moving_index;
    voronoi_map_->TransformPhysicalPointToIndex( mapped_moving_itk, mapped_moving_index );

    for( unsigned int d = 0; d < 3; ++d ) {
      if( mapped_moving_index[d] < 0 ) mapped_moving_index[d] = 0;
      if( mapped_moving_index[d] >= voronoiSize[d] ) mapped_moving_index[d] = VoronoiMapType::IndexValueType( voronoiSize[d] - 1 );
    }
    unsigned int fixed_ind = voronoi_map_->GetPixel( mapped_moving_index );

    feature_sptr_type  q_itk = fixed_[res_][fixed_ind];


    //vcl_cout << "MAP: " << q_itk->location_ << vcl_endl; 
    
  }

  long time_map = timer_map.real();


  vcl_cout << "MAP  TREE: " << time_map << "  " << time_tree << vcl_endl; 
}


int
main( int argc, char* argv[] )
{
  vul_arg< std::string >        fixedImageDir            ( 0, "Fixed Dicom volume directory" );
  vul_arg< std::string >        movingImageDir           ( "-moving", "Moving Dicom volume directory", "" );

  vul_arg_parse( argc, argv );

  typedef cdcl_feature_with_shape< 3 >::sptr  feature_sptr_type;

  vcl_vector< vcl_vector< feature_sptr_type > >  fixed_;
  ReadFeatures( fixedImageDir(), fixed_ );

  typedef itk::Image< unsigned int, 3 >  VoronoiMapType;
  VoronoiMapType::Pointer  fixed_voronoi_map_;

  typedef   itk::Image< unsigned int, 3 >  UnsignedIntImageType;

#if 1
  // Compute distance map
  compute_distance_map( fixed_, fixed_voronoi_map_ );
#else
  typedef   itk::ImageFileReader< UnsignedIntImageType >  MapFileReaderType;
  MapFileReaderType::Pointer  mapFileReader = MapFileReaderType::New();
  std::string  fileName = fixedImageDir() + "voronoi.mhd";
  mapFileReader->SetFileName( fileName );
  mapFileReader->Update();
  fixed_voronoi_map_ = mapFileReader->GetOutput();
#endif

  if( movingImageDir.set() && movingImageDir() != "" ) {
    vcl_vector< vcl_vector< feature_sptr_type > >  moving_;
    ReadFeatures( movingImageDir(), moving_ );

    CompareMapTree( fixed_voronoi_map_, fixed_, moving_ );
  }

  typedef   itk::ImageFileWriter< UnsignedIntImageType >  MapFileWriterType;
  MapFileWriterType::Pointer  mapFileWriter = MapFileWriterType::New();
  std::string  fileNameWrite = fixedImageDir() + "voronoi.mhd";
  mapFileWriter->SetFileName( fileNameWrite );
  mapFileWriter->SetInput( fixed_voronoi_map_ );
  mapFileWriter->Update();


}
