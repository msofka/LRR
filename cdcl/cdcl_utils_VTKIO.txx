#ifndef cdcl_utils_VTKIO_txx_
#define cdcl_utils_VTKIO_txx_

#include "cdcl_feature.h"
#include "cdcl_keypoint.h"
#include "cdcl_match.h"
#include "cdcl_trans.h"
#include "cdcl_utils_VTKIO.h"

#include <vcl_string.h>

#include "vtkPolyData.h"
#include "vtkSmartPointer.h"



#include "vtkPointData.h"
#include "vtkDoubleArray.h"
#include "vtkFloatArray.h"

//:
// \file
// \brief  Writing matches, points, and covariances (axes of the ellipse or ellipsoid).
// \author Michal Sofka
// \date   May 2006


#define USE_FLOAT 1

// Write keypoints with their descriptors in a VTK format
template < unsigned int dim >
void
cdcl_keypoints_to_poly_data_VTK( vcl_vector< vbl_smart_ptr< cdcl_keypoint< dim > > >  const &  keypoints,
                                 vtkSmartPointer< vtkPolyData >                             &  poly_data )
{
  typedef typename cdcl_feature< dim >::sptr      feature_sptr_type;

#if USE_FLOAT
  vtkSmartPointer< vtkPoints >      vtkpoints   = vtkSmartPointer< vtkPoints >::New();
  vtkSmartPointer< vtkFloatArray >  covariances = vtkSmartPointer< vtkFloatArray >::New();
  vtkSmartPointer< vtkFloatArray >  strengths   = vtkSmartPointer< vtkFloatArray >::New();
  vtkSmartPointer< vtkFloatArray >  normals     = vtkSmartPointer< vtkFloatArray >::New();
#else
  vtkSmartPointer< vtkPoints >       vtkpoints   = vtkSmartPointer< vtkPoints >::New();
  vtkSmartPointer< vtkDoubleArray >  covariances = vtkSmartPointer< vtkDoubleArray >::New();
  vtkSmartPointer< vtkDoubleArray >  strengths   = vtkSmartPointer< vtkDoubleArray >::New();
  vtkSmartPointer< vtkDoubleArray >  normals     = vtkSmartPointer< vtkDoubleArray >::New();
#endif

  normals->SetName( "normals" );
  normals->SetNumberOfComponents( keypoints[0]->normal_.size() );

  covariances->SetName( "covariances" );
  covariances->SetNumberOfComponents( dim*dim );

  strengths->SetName( "strengths" );
  strengths->SetNumberOfComponents( 1 );

  // write data points
  //
  for( typename vcl_vector< feature_sptr_type >::size_type  i = 0; i < keypoints.size(); ++i ) {

#if USE_FLOAT
    vnl_vector_fixed< float, dim >  location_float;
    location_float[0] = float( keypoints[i]->location_[0] );
    location_float[1] = float( keypoints[i]->location_[1] );
    location_float[2] = float( keypoints[i]->location_[2] );

    vnl_matrix_fixed< float, dim, dim >  covariance_float;
    for( unsigned int d1 = 0; d1 < dim; ++d1 )
      for( unsigned int d2 = 0; d2 < dim; ++d2 ) {
        covariance_float( d1, d2 ) = float( keypoints[i]->covariance_( d1, d2 ) );
      }

    vnl_vector_fixed< float, dim >  normal_float;
    normal_float[0] = float( keypoints[i]->normal_[0] );
    normal_float[1] = float( keypoints[i]->normal_[1] );
    normal_float[2] = float( keypoints[i]->normal_[2] );

    vtkpoints->InsertNextPoint( location_float.data_block() );
    covariances->InsertNextTupleValue( covariance_float.data_block() );
    normals->InsertNextTupleValue( normal_float.data_block() );
#else
    vtkpoints->InsertNextPoint( keypoints[i]->location_.data_block() );
    covariances->InsertNextTupleValue( keypoints[i]->covariance_.data_block() );
    normals->InsertNextTupleValue( keypoints[i]->normal_.data_block() );
#endif

    //vnl_matrix_fixed< double, dim, dim >  outer_product;
    //outer_product = vnl_inverse( keypoints[i]->covariance_ );

    ////// eigenvectors and eigenvalues of the covariance matrix
    ////vnl_matrix< double >  V( dim, dim );
    ////vnl_vector< double >  D( dim );
    ////D.squared_magnitude();

    ////vnl_symmetric_eigensystem_compute( outer_product, V, D );

    //double strength = vnl_trace( outer_product );
    ////vcl_cout << "strength: " << strength << vcl_endl;
    //strengths->InsertNextValue( strength );
    strengths->InsertNextValue( keypoints[i]->strength_ );
  }

  poly_data->SetPoints( vtkpoints );
  poly_data->GetPointData()->SetScalars( strengths );
  poly_data->GetPointData()->SetVectors( normals );
  poly_data->GetPointData()->SetTensors( covariances );

}


// Write keypoints with their descriptors in a VTK format
VCL_DEFINE_SPECIALIZATION
void
cdcl_keypoints_to_poly_data_VTK( vcl_vector< vbl_smart_ptr< cdcl_keypoint< 3 > > >  const &  keypoints,
                                 vtkSmartPointer< vtkPolyData >                           &  poly_data )
{
  const unsigned int dim = 3;
  typedef cdcl_feature< dim >::sptr      feature_sptr_type;

  vtkSmartPointer< vtkPoints >       vtkpoints   = vtkSmartPointer< vtkPoints >::New();

#if USE_FLOAT
  vtkSmartPointer< vtkFloatArray >  covariances = vtkSmartPointer< vtkFloatArray >::New();
  vtkSmartPointer< vtkFloatArray >  strengths   = vtkSmartPointer< vtkFloatArray >::New();
  vtkSmartPointer< vtkFloatArray >  normals     = vtkSmartPointer< vtkFloatArray >::New();
  vtkSmartPointer< vtkFloatArray >  binormals   = vtkSmartPointer< vtkFloatArray >::New();
#else
  vtkSmartPointer< vtkDoubleArray >  covariances = vtkSmartPointer< vtkDoubleArray >::New();
  vtkSmartPointer< vtkDoubleArray >  strengths   = vtkSmartPointer< vtkDoubleArray >::New();
  vtkSmartPointer< vtkDoubleArray >  normals     = vtkSmartPointer< vtkDoubleArray >::New();
  vtkSmartPointer< vtkDoubleArray >  binormals   = vtkSmartPointer< vtkDoubleArray >::New();
#endif

  normals->SetName( "normals" );
  normals->SetNumberOfComponents( keypoints[0]->normal_.size() );

  binormals->SetName( "binormals" );
  binormals->SetNumberOfComponents( keypoints[0]->binormal_.size() );
 
  covariances->SetName( "covariances" );
  covariances->SetNumberOfComponents( dim*dim );

  strengths->SetName( "strengths" );
  strengths->SetNumberOfComponents( 1 );

  // write data points
  //
  for( vcl_vector< feature_sptr_type >::size_type  i = 0; i < keypoints.size(); ++i ) {

#if USE_FLOAT
    vnl_vector_fixed< float, dim >  location_float;
    location_float[0] = float( keypoints[i]->location_[0] );
    location_float[1] = float( keypoints[i]->location_[1] );
    location_float[2] = float( keypoints[i]->location_[2] );

    vnl_matrix_fixed< float, dim, dim >  covariance_float;
    for( unsigned int d1 = 0; d1 < dim; ++d1 )
      for( unsigned int d2 = 0; d2 < dim; ++d2 ) {
        covariance_float( d1, d2 ) = float( keypoints[i]->covariance_( d1, d2 ) );
      }

    vnl_vector_fixed< float, dim >  normal_float;
    normal_float[0] = float( keypoints[i]->normal_[0] );
    normal_float[1] = float( keypoints[i]->normal_[1] );
    normal_float[2] = float( keypoints[i]->normal_[2] );

    vnl_vector_fixed< float, dim >  binormal_float;
    binormal_float[0] = float( keypoints[i]->binormal_[0] );
    binormal_float[1] = float( keypoints[i]->binormal_[1] );
    binormal_float[2] = float( keypoints[i]->binormal_[2] );

    vtkpoints->InsertNextPoint( location_float.data_block() );
    covariances->InsertNextTupleValue( covariance_float.data_block() );
    normals->InsertNextTupleValue( normal_float.data_block() );

    strengths->InsertNextValue( keypoints[i]->strength_ );

    binormals->InsertNextTupleValue( binormal_float.data_block() );

#else
    vtkpoints->InsertNextPoint( keypoints[i]->location_.data_block() );
    covariances->InsertNextTupleValue( keypoints[i]->covariance_.data_block() );

    //vnl_matrix_fixed< double, dim, dim >  outer_product;
    //outer_product = vnl_inverse( keypoints[i]->covariance_ );

    ////// eigenvectors and eigenvalues of the covariance matrix
    ////vnl_matrix< double >  V( dim, dim );
    ////vnl_vector< double >  D( dim );
    ////D.squared_magnitude();

    ////vnl_symmetric_eigensystem_compute( outer_product, V, D );

    //double strength = vnl_trace( outer_product );
    ////vcl_cout << "strength: " << strength << vcl_endl;
    //strengths->InsertNextValue( strength );
    strengths->InsertNextValue( keypoints[i]->strength_ );

    normals->InsertNextTupleValue( keypoints[i]->normal_.data_block() );

    binormals->InsertNextTupleValue( keypoints[i]->binormal_.data_block() );
#endif
  }

  poly_data->SetPoints( vtkpoints );
  poly_data->GetPointData()->SetScalars( strengths );
  poly_data->GetPointData()->SetActiveVectors( "normals" );
  poly_data->GetPointData()->SetVectors( normals );
  poly_data->GetPointData()->SetActiveVectors( "binormals" );
  poly_data->GetPointData()->SetVectors( binormals );
  poly_data->GetPointData()->SetTensors( covariances );

}


// Write keypoints with their descriptors in a VTK format
template < unsigned int dim >
void
cdcl_poly_data_to_keypoints_VTK( vtkSmartPointer< vtkPolyData >                       const &  poly_data,
                                 vcl_vector< vbl_smart_ptr< cdcl_keypoint< dim > > >        &  keypoints )
// THIS IS 2D, 3D is below
{
  // get points and covariances
  vtkSmartPointer< vtkPoints >       vtkpoints      = poly_data->GetPoints();
  vtkSmartPointer< vtkDataArray >    covariances    = poly_data->GetPointData()->GetTensors();
  vtkSmartPointer< vtkDoubleArray >  strengths      = vtkDoubleArray::SafeDownCast( poly_data->GetPointData()->GetScalars( "strengths" ) );
  vtkSmartPointer< vtkDoubleArray >  normals        = vtkDoubleArray::SafeDownCast( poly_data->GetPointData()->GetVectors( "normals" ) );
  vtkSmartPointer< vtkDoubleArray >  binormals      = vtkDoubleArray::SafeDownCast( poly_data->GetPointData()->GetVectors( "binormals" ) );


  // store the data
  keypoints.clear();
  for( vtkIdType  i = 0; i < vtkpoints->GetNumberOfPoints(); ++i ) {
    vnl_vector_fixed< double, dim >  location;
    vtkpoints->GetPoint( i, location.data_block() );

    vnl_matrix_fixed< double, dim, dim >  covariance;
    covariances->GetTuple( i, covariance.data_block() );

    typename cdcl_feature<dim>::sptr  feature = new cdcl_feature< dim >( location, covariance );

    vnl_vector_fixed< double, dim >  normal;
    normals->GetTupleValue( i, normal.data_block() );

    double strength;
    strengths->GetTuple( i, &strength );

    typename cdcl_keypoint<dim>::sptr  keypoint = new cdcl_keypoint<dim>( location, strength, covariance, normal );
    keypoints.push_back( keypoint );
  }

}


VCL_DEFINE_SPECIALIZATION
void
cdcl_poly_data_to_keypoints_VTK( vtkSmartPointer< vtkPolyData >                     const &  poly_data,
                                 vcl_vector< vbl_smart_ptr< cdcl_keypoint< 3 > > >        &  keypoints )
{
  // get points and covariances
  vtkSmartPointer< vtkPoints >       vtkpoints      = poly_data->GetPoints();
  vtkSmartPointer< vtkDataArray >    covariances    = poly_data->GetPointData()->GetTensors();

#if USE_FLOAT
  vtkSmartPointer< vtkFloatArray >   strengths      = vtkFloatArray::SafeDownCast( poly_data->GetPointData()->GetScalars( "strengths" ) );
  poly_data->GetPointData()->SetActiveVectors( "normals" );
  vtkSmartPointer< vtkFloatArray >   normals        = vtkFloatArray::SafeDownCast( poly_data->GetPointData()->GetVectors( "normals" ) );
  poly_data->GetPointData()->SetActiveVectors( "binormals" );
  vtkSmartPointer< vtkFloatArray >   binormals      = vtkFloatArray::SafeDownCast( poly_data->GetPointData()->GetVectors( "binormals" ) );
#else
  vtkSmartPointer< vtkDoubleArray >   strengths      = vtkDoubleArray::SafeDownCast( poly_data->GetPointData()->GetScalars( "strengths" ) );
  poly_data->GetPointData()->SetActiveVectors( "normals" );
  vtkSmartPointer< vtkDoubleArray >   normals        = vtkDoubleArray::SafeDownCast( poly_data->GetPointData()->GetVectors( "normals" ) );
  poly_data->GetPointData()->SetActiveVectors( "binormals" );
  vtkSmartPointer< vtkDoubleArray >   binormals      = vtkDoubleArray::SafeDownCast( poly_data->GetPointData()->GetVectors( "binormals" ) );
#endif

  if( !normals || !binormals ) {
    std::cout << "ERROR: The data coord type float/double does not match that of the keypoint reader. File: " 
              << __FILE__ << " on line " << __LINE__ << std::endl;
    return;
  }

  const unsigned int dim = 3;

  // store the data
  keypoints.clear();
  for( vtkIdType  i = 0; i < vtkpoints->GetNumberOfPoints(); ++i ) {
    vnl_vector_fixed< double, dim >  location;
    vtkpoints->GetPoint( i, location.data_block() );

    vnl_matrix_fixed< double, dim, dim >  covariance( 0.0 );
    if( covariances ) covariances->GetTuple( i, covariance.data_block() );

    //double strength;
    //strengths->GetTuple( i, &strength );
    double strength = 0.0;
    if( strengths ) strengths->GetTuple( i, &strength );
#if USE_FLOAT
    vnl_vector_fixed< float, dim >  normal_float;
    normals->GetTupleValue( i, normal_float.data_block() );

    vnl_vector_fixed< float, dim >  binormal_float;
    binormals->GetTupleValue( i, binormal_float.data_block() );

    vnl_vector_fixed< double, dim >  normal;
    normal[0] = normal_float[0];
    normal[1] = normal_float[1];
    normal[2] = normal_float[2];

    vnl_vector_fixed< double, dim >  binormal;
    binormal[0] = binormal_float[0];
    binormal[1] = binormal_float[1];
    binormal[2] = binormal_float[2];
#else
    vnl_vector_fixed< double, dim >  normal;
    normals->GetTupleValue( i, normal.data_block() );

    vnl_vector_fixed< double, dim >  binormal;
    binormals->GetTupleValue( i, binormal.data_block() );
#endif
    cdcl_keypoint<dim>::sptr  keypoint = new cdcl_keypoint<dim>( location, strength, covariance, normal, binormal );

    keypoints.push_back( keypoint );
  }

}



#define CDCL_UTILS_VTKIO_INSTANTIATE( dim )                                                         \
template                                                                                             \
void                                                                                                        \
cdcl_poly_data_to_keypoints_VTK( vtkSmartPointer< vtkPolyData >                       const &  poly_data,   \
                                 vcl_vector< vbl_smart_ptr< cdcl_keypoint< dim > > >        &  keypoints ); \
template                                                                                                    \
void                                                                                                        \
cdcl_keypoints_to_poly_data_VTK( vcl_vector< vbl_smart_ptr< cdcl_keypoint< dim > > >  const &  keypoints,   \
                                 vtkSmartPointer< vtkPolyData >                             &  poly_data );

#endif
