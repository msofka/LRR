#ifndef cdcl_utils_VTK_txx_
#define cdcl_utils_VTK_txx_

#include <sstream>
#include <iomanip>

#include <vnl/algo/vnl_symmetric_eigensystem.h>
#include <vnl/vnl_transpose.h>

#include <vnl/vnl_trace.h>
#include <vnl/vnl_inverse.h>

#include "vtkPoints.h"
#include "vtkTensor.h"
#include "vtkPolyData.h"
#include "vtkPointData.h"
#include "vtkVertex.h"
#include "vtkLine.h"
#include "vtkIntArray.h"
#include "vtkFloatArray.h"
#include "vtkDataArray.h"
#include "vtkDoubleArray.h"
#include "vtkFloatArray.h"
#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkXMLPolyDataWriter.h"
#include <vtkXMLPolyDataReader.h>
#include "vtkSmartPointer.h"

#include "cdcl_utils_VTK.h"
#include "cdcl_utils_VTKIO.h"

//:
// \file
// \brief  Writing matches, points, and covariances (axes of the ellipse or ellipsoid).
// \author Michal Sofka
// \date   May 2006


#define MATCHES 0
#define COVARIANCES 0

namespace {

inline bool is_positive_definite( vnl_matrix_fixed< double, 3, 3 > const &  covariance )
{
  // The following rule can be rewritten in terms of determinants -> less expensive
  //
  // positive definite: all major subdeterminants are positive
  // negative definite: all major subdeterminants change signs, beginning with negative

  // eigen decomposition of the covariance matrix
  double eigenValues[3];
  vnl_symmetric_eigensystem_compute_eigenvals( covariance( 0, 0 ), covariance( 0, 1 ), covariance( 0, 2 ),
                                                                   covariance( 1, 1 ), covariance( 1, 2 ),
                                                                                       covariance( 2, 2 ),
                                               eigenValues[0], eigenValues[1], eigenValues[2] );

  if( eigenValues[0] > 0 && eigenValues[1] > 0 && eigenValues[2] > 0 ) {
      //std::cout << "Warning: Current point is not a maximum (matrix not negative definite)." << std::endl;
      return true;
  }

  return false;
}

}


namespace {

inline bool is_positive_definite( vnl_matrix_fixed< double, 2, 2 > const &  covariance )
{
  // The following rule can be rewritten in terms of determinants -> less expensive
  //
  // positive definite: all major subdeterminants are positive
  // negative definite: all major subdeterminants change signs, beginning with negative

  // eigen decomposition of the covariance matrix
  vnl_matrix<double>  V;
  vnl_vector<double>  D;

  vnl_symmetric_eigensystem_compute( covariance.as_matrix(), V, D );

  if( D[0] > 0 && D[1] > 0 ) {
      //std::cout << "Warning: Current point is not a maximum (matrix not negative definite)." << std::endl;
      return true;
  }

  return false;
}

}


// Write matches in a VTK format
template < unsigned int dim, unsigned int dof >
void
cdcl_write_matches_VTK( vcl_vector< vbl_smart_ptr< cdcl_feature< dim > > > const &  moving,
                        vcl_vector< vbl_smart_ptr< cdcl_feature< dim > > > const &  fixed,
                        vcl_vector< vbl_smart_ptr< cdcl_match< dim > > >   const &  matches,
                        vbl_smart_ptr< cdcl_trans< dim, dof > >            const &  trans,
                        unsigned int                                                iteration )
{
  typedef typename cdcl_feature< dim >::sptr      feature_sptr_type;
  typedef typename cdcl_match< dim >::sptr        match_sptr_type;
  typedef typename cdcl_trans< dim, dof >::sptr   trans_sptr_type;
  typedef typename vcl_vector< match_sptr_type >  match_set_type;

  vtkSmartPointer< vtkPoints >     pts      = vtkSmartPointer< vtkPoints >::New();
  vtkSmartPointer< vtkIntArray >   set      = vtkSmartPointer< vtkIntArray >::New();
  vtkSmartPointer< vtkFloatArray > weights  = vtkSmartPointer< vtkFloatArray >::New();
  vtkSmartPointer< vtkCellArray >  lines    = vtkSmartPointer< vtkCellArray >::New();
  vtkSmartPointer< vtkCellArray >  vertices = vtkSmartPointer< vtkCellArray >::New();

  vcl_ostringstream  suffix;
  suffix << std::setw( 3 ) << std::setfill( '0' ) << iteration;

  set->SetNumberOfComponents( 1 );
  set->SetName("0moving-1fixed");

  weights->SetNumberOfComponents( 1 );
  weights->SetName( "Weights" );

  // id's for q (fixed) and p (moving) points
  unsigned int nq = 0;
  unsigned int np = 0;
  unsigned int n = 0;

  // write matches
  //
  for( typename vcl_vector< match_sptr_type >::size_type  i = 0; i < matches.size(); ++i ) {
    // moving point p
    feature_sptr_type  p = matches[i]->from_;

    // mapped moving point
    vnl_vector_fixed< double, dim >  const &  mapped_moving = trans->map_loc( p->location_ );

    //pts->InsertNextPoint( moving_mapped.data_block() );

    // np is current index, start inserting fixed points from index nq
    nq = np + 1;

    // go through all matches (fixed points)
    for( typename vcl_vector< match_sptr_type >::size_type  j = 0; j < matches[i]->to_.size(); ++j, ++nq ) {
      // fixed point q
      feature_sptr_type  q = matches[i]->to_[j];
   
      pts->InsertNextPoint( mapped_moving.data_block() );
      pts->InsertNextPoint( q->location_.data_block() );
      vtkIdType lineIds[2] = {n, n+1};

      //vtkIdType lineIds[2] = {np, nq};
      lines->InsertNextCell( 2, lineIds );  // number of Ids and an array
      // the previous is significantly faster than construction below
      //vtkSmartPointer< vtkLine > line = vtkSmartPointer< vtkLine >::New();
      //line->GetPointIds()->SetId( 0, n );
      //line->GetPointIds()->SetId( 1, n+1 );

      // 0 = moving, 1 = fixed
      set->InsertNextValue( 0 );
      set->InsertNextValue( 1 );

      weights->InsertNextValue( matches[i]->w_[j] );

      n += 2;

    }
    // next, moving point will be inserted, so prepare its index
    np = nq + 1;
  }

  // mixed cells must be inserted in the order of vertices (vtkVertex and vtkPolyVertex), lines (vtkLine and vtkPolyLine), polygons (vtkTriangle, vtkQuad, vtkPolygon), and triangle strips (vtkTriangleStrip)
  // when adding scalars for lines and using vertices and lines, we need add scalars also for vertices (fill them with arbitrary value)
  vtkSmartPointer< vtkPolyData > pdata = vtkSmartPointer< vtkPolyData >::New();
  pdata->SetPoints( pts );
  pdata->GetPointData()->SetScalars( set );
  pdata->SetLines( lines );
  pdata->GetCellData()->SetScalars( weights );

  vtkSmartPointer< vtkXMLPolyDataWriter > writer = vtkSmartPointer< vtkXMLPolyDataWriter >::New();

  if( MATCHES ) {
    writer->SetInput( pdata );
    vcl_string filename = "matches" + suffix.str() + ".vtp";
    writer->SetFileName( filename.c_str() );
    writer->Write();
  }


  vtkSmartPointer< vtkPoints >     cov_pts   = vtkSmartPointer< vtkPoints >::New();


  // id's for q (fixed) and p (moving) points
  n = 0;

  // write covariances
  //
  // Jacobian w.r.t. p (moving point)
  vnl_matrix_fixed< double, dim, dim > const &  Jp = trans->jacobian_wrt_loc();
  for( typename vcl_vector< match_sptr_type >::size_type  i = 0; i < matches.size(); ++i ) {
    // moving point p
    feature_sptr_type  p = matches[i]->from_;

    // mapped moving point
    vnl_vector_fixed< double, dim >  const &  mapped_moving = trans->map_loc( p->location_ );

    // Jacobian w.r.t. transformation parameters
    vnl_matrix_fixed< double, dim, dof > const &  Jth = trans->jacobian_wrt_par( p->location_ );

    // eigenvectors and eigenvalues of the covariance matrix
    vnl_matrix< double >  V( dim, dim );
    vnl_vector< double >  D( dim );

    // go through all matches (fixed points)
    for( typename vcl_vector< match_sptr_type >::size_type  j = 0; j < matches[i]->to_.size(); ++j, ++nq ) {
      // fixed point q
      feature_sptr_type  q = matches[i]->to_[j];
  
      // covariance of a correspondence match
      vnl_matrix_fixed< double, dim, dim >  Cij = Jth * trans->get_covariance() * Jth.transpose() +
                                                  Jp * p->covariance_ * Jp.transpose() + 
                                                  q->covariance_;

      // eigen decomposition of the covariance matrix
      const vnl_matrix< double >  A = Cij;
      vnl_symmetric_eigensystem_compute( Cij, V, D );
    
      vnl_vector_fixed< double, dim >  mid = q->location_ + ( ( mapped_moving - q->location_ ) * 0.5 );

      for( unsigned int d = 0; d < dim; ++d ) {
        // square root of eigenvalues
        double e = vcl_sqrt( D(d) );

        // end points of covariance crosses
        vnl_vector_fixed< double, dim >  P1 = mid + e * V.get_column(d) * 0.5;
        vnl_vector_fixed< double, dim >  P2 = mid - e * V.get_column(d) * 0.5;

        cov_pts->InsertNextPoint( P1.data_block() );
        cov_pts->InsertNextPoint( P2.data_block() );

        vtkIdType lineIds[2] = {n, n+1};
        lines->InsertNextCell( 2, lineIds );  // number of Ids and an array
        n += 2;
      }

    }
  }

  // mixed cells must be inserted in the order of vertices (vtkVertex and vtkPolyVertex), lines (vtkLine and vtkPolyLine), polygons (vtkTriangle, vtkQuad, vtkPolygon), and triangle strips (vtkTriangleStrip)
  // when adding scalars for lines and using vertices and lines, we need add scalars also for vertices (fill them with arbitrary value)
  vtkSmartPointer< vtkPolyData > cov_data = vtkSmartPointer< vtkPolyData >::New();
  cov_data->SetPoints( cov_pts );
//  cov_data->GetPointData()->SetScalars( set );
  cov_data->SetLines( lines );
//  cov_data->GetCellData()->SetScalars( weights );

  if( COVARIANCES ) {
    writer->SetInput( cov_data );
    vcl_string filename = "covariances" + suffix.str() + ".vtp";
    writer->SetFileName( filename.c_str() );
    writer->Write();
  }


  // write fixed data points
  //
  vtkSmartPointer< vtkPoints > fixed_pts = vtkSmartPointer< vtkPoints >::New();
  for( typename vcl_vector< feature_sptr_type >::size_type  i = 0; i < fixed.size(); ++i ) {
    fixed_pts->InsertNextPoint( fixed[i]->location_.data_block() );
    vtkIdType vertexId[1] = { i };
    vertices->InsertNextCell( 1, vertexId );
  }

  vtkSmartPointer< vtkPolyData > fixed_data = vtkSmartPointer< vtkPolyData >::New();
  fixed_data->SetPoints( fixed_pts );
  //fixed_data->SetVerts( vertices );

  writer->SetInput( fixed_data );
  vcl_string filename = "fixed" + suffix.str() + ".vtp";
  writer->SetFileName( filename.c_str() );
  writer->Write();



  // write moving data points
  //
  vtkSmartPointer< vtkPoints > moving_pts = vtkSmartPointer< vtkPoints >::New();
  for( typename vcl_vector< feature_sptr_type >::size_type  i = 0; i < moving.size(); ++i ) {
    // mapped moving point
    vnl_vector_fixed< double, dim >  const &  moving_mapped = trans->map_loc( moving[i]->location_ );
    moving_pts->InsertNextPoint( moving_mapped.data_block() );
    vtkIdType vertexId[1] = { i };
    vertices->InsertNextCell( 1, vertexId );
  }

  vtkSmartPointer< vtkPolyData > moving_data = vtkSmartPointer< vtkPolyData >::New();
  moving_data->SetPoints( moving_pts );
  //moving_data->SetVerts( vertices );

  writer->SetInput( moving_data );
  filename = "moving" + suffix.str() + ".vtp";
  writer->SetFileName( filename.c_str() );
  writer->Write();

}


// Write keypoints with their descriptors in a VTK format
template < unsigned int dim >
void
cdcl_write_keypoint_descriptors_VTK( vcl_vector< vbl_smart_ptr< cdcl_keypoint< dim > > >  const &  keypoints,
                                     vcl_vector< vnl_vector< float > >                    const &  descriptors,
                                     vcl_string                                           const    filename )
{
  if( keypoints.size() != descriptors.size() ) {
    vcl_cout << "Error: vectors of keypoints has size: " << keypoints.size() 
             << " and vector of descriptors have size: " << descriptors.size() << vcl_endl;
    return;
  }

  vtkSmartPointer< vtkXMLPolyDataWriter > writer = vtkSmartPointer< vtkXMLPolyDataWriter >::New();

  vtkSmartPointer< vtkPolyData >  poly_data = vtkSmartPointer< vtkPolyData >::New();

  cdcl_keypoints_to_poly_data_VTK( keypoints, poly_data );

  typedef typename cdcl_feature< dim >::sptr      feature_sptr_type;

  vtkSmartPointer< vtkFloatArray >  descriptorData = vtkSmartPointer< vtkFloatArray >::New();
  descriptorData->SetName( "descriptors" );
  descriptorData->SetNumberOfComponents( descriptors[0].size() );

  // write data points
  //
  for( typename vcl_vector< vnl_vector< float > >::size_type  i = 0; i < descriptors.size(); ++i ) {
    descriptorData->InsertNextTupleValue( descriptors[i].data_block() );
  }
  poly_data->GetPointData()->SetScalars( descriptorData );

  writer->SetInput( poly_data );
  writer->SetFileName( filename.c_str() );
  writer->Write();
}


// Write keypoints with their descriptors in a VTK format
template < unsigned int dim >
bool
cdcl_read_keypoint_descriptors_VTK( vcl_vector< vbl_smart_ptr< cdcl_keypoint< dim > > >        &  keypoints,
                                    vcl_vector< vnl_vector< float > >                          &  descriptors,
                                    vcl_string                                           const    filename )
{
  typedef typename cdcl_feature< dim >::sptr                                 feature_sptr_type;

  // read poly data
  vtkSmartPointer< vtkXMLPolyDataReader >  polyDataReader = vtkXMLPolyDataReader::New();

  polyDataReader->SetFileName( filename.c_str() );
  int file_exists = polyDataReader->CanReadFile( filename.c_str() );
  if( file_exists ) {
    vcl_cout << "Loading " << filename << vcl_endl;
  }
  else {
    vcl_cout << "Error: cannot read points: " << filename.c_str() << vcl_endl;
    return false;
  }
  polyDataReader->Update();
  
  vtkSmartPointer< vtkPolyData >  poly_data = polyDataReader->GetOutput();

  cdcl_poly_data_to_keypoints_VTK( poly_data, keypoints );

  // read descriptors
  //
  // get points and covariances
  vtkSmartPointer< vtkPoints >       vtkpoints      = poly_data->GetPoints();
  vtkSmartPointer< vtkFloatArray >   descriptorData = vtkFloatArray::SafeDownCast( poly_data->GetPointData()->GetScalars( "descriptors" ) );

  // store the data
  for( vtkIdType  i = 0; i < vtkpoints->GetNumberOfPoints(); ++i ) {

    vnl_vector< float >  descriptor( descriptorData->GetNumberOfComponents() );
    descriptorData->GetTupleValue( i, descriptor.data_block() );
    descriptors.push_back( descriptor );
  }

  return true;
}


// Convert cdcl features to VTK polydata for display purposes
template < unsigned int dim >
void
cdcl_features_to_poly_data( vcl_vector< vbl_smart_ptr< cdcl_feature< dim > > > const &  points,
                            vtkSmartPointer< vtkPolyData >                           &  poly_data )
{
  typedef typename cdcl_feature< dim >::sptr      feature_sptr_type;

  vtkSmartPointer< vtkPoints >       vtkpoints   = vtkSmartPointer< vtkPoints >::New();
  vtkSmartPointer< vtkDoubleArray >  covariances = vtkSmartPointer< vtkDoubleArray >::New();
  vtkSmartPointer< vtkDoubleArray >  strengths   = vtkSmartPointer< vtkDoubleArray >::New();

  covariances->SetName( "covariances" );
  covariances->SetNumberOfComponents( dim*dim );

  strengths->SetName( "strengths" );
  strengths->SetNumberOfComponents( 1 );

  // write data points
  //
  for( typename vcl_vector< feature_sptr_type >::size_type  i = 0; i < points.size(); ++i ) {
    vtkpoints->InsertNextPoint( points[i]->location_.data_block() );
    covariances->InsertNextTupleValue( points[i]->covariance_.data_block() );

    //vnl_matrix_fixed< double, dim, dim >  outer_product;
    //outer_product = vnl_inverse( points[i]->covariance_ );

    ////// eigenvectors and eigenvalues of the covariance matrix
    ////vnl_matrix< double >  V( dim, dim );
    ////vnl_vector< double >  D( dim );
    ////D.squared_magnitude();

    ////vnl_symmetric_eigensystem_compute( outer_product, V, D );

    //double strength = vnl_trace( outer_product );
    ////vcl_cout << "strength: " << strength << vcl_endl;
    //strengths->InsertNextValue( strength );
    strengths->InsertNextValue( points[i]->strength_ );
  }

  if( !poly_data ) poly_data = vtkSmartPointer< vtkPolyData >::New();
  else poly_data->Reset();

  poly_data->SetPoints( vtkpoints );
  poly_data->GetPointData()->SetScalars( strengths );
  poly_data->GetPointData()->SetTensors( covariances );
}


// Convert cdcl features to VTK polydata for display purposes
template < unsigned int dim >
void
cdcl_features_to_poly_data( vcl_vector< vbl_smart_ptr< cdcl_feature_with_shape< dim > > > const &  points,
                            vtkSmartPointer< vtkPolyData >                                      &  poly_data )
{
  typedef typename cdcl_feature_with_shape< dim >::sptr      feature_sptr_type;

  vtkSmartPointer< vtkPoints >       vtkpoints   = vtkSmartPointer< vtkPoints >::New();
  vtkSmartPointer< vtkDoubleArray >  covariances = vtkSmartPointer< vtkDoubleArray >::New();
  vtkSmartPointer< vtkDoubleArray >  strengths   = vtkSmartPointer< vtkDoubleArray >::New();

  covariances->SetName( "covariances" );
  covariances->SetNumberOfComponents( dim*dim );

  strengths->SetName( "strengths" );
  strengths->SetNumberOfComponents( 1 );

  // write data points
  //
  for( typename vcl_vector< feature_sptr_type >::size_type  i = 0; i < points.size(); ++i ) {
    vtkpoints->InsertNextPoint( points[i]->location_.data_block() );
    covariances->InsertNextTupleValue( points[i]->covariance_.data_block() );

    //vnl_matrix_fixed< double, dim, dim >  outer_product;
    //outer_product = vnl_inverse( points[i]->covariance_ );

    ////// eigenvectors and eigenvalues of the covariance matrix
    ////vnl_matrix< double >  V( dim, dim );
    ////vnl_vector< double >  D( dim );
    ////D.squared_magnitude();

    ////vnl_symmetric_eigensystem_compute( outer_product, V, D );

    //double strength = vnl_trace( outer_product );
    ////vcl_cout << "strength: " << strength << vcl_endl;
    //strengths->InsertNextValue( strength );
    strengths->InsertNextValue( points[i]->strength_ );
  }

  if( !poly_data ) poly_data = vtkSmartPointer< vtkPolyData >::New();
  else poly_data->Reset();

  poly_data->SetPoints( vtkpoints );
  poly_data->GetPointData()->SetScalars( strengths );
  poly_data->GetPointData()->SetTensors( covariances );
}


// Convert cdcl features to VTK polydata for display purposes
template < unsigned int dim >
void
cdcl_features_to_poly_data( vcl_vector< vbl_smart_ptr< cdcl_feature_ICP< dim > > > const &  points,
                            vtkSmartPointer< vtkPolyData >                               &  poly_data )
{
  typedef typename cdcl_feature_ICP< dim >::sptr      feature_sptr_type;

  vtkSmartPointer< vtkPoints >       vtkpoints   = vtkSmartPointer< vtkPoints >::New();
  //vtkSmartPointer< vtkDoubleArray >  covariances = vtkSmartPointer< vtkDoubleArray >::New();
  vtkSmartPointer< vtkDoubleArray >  strengths   = vtkSmartPointer< vtkDoubleArray >::New();

  //covariances->SetName( "covariances" );
  //covariances->SetNumberOfComponents( dim*dim );

  strengths->SetName( "strengths" );
  strengths->SetNumberOfComponents( 1 );

  // write data points
  //
  for( typename vcl_vector< feature_sptr_type >::size_type  i = 0; i < points.size(); ++i ) {
    vtkpoints->InsertNextPoint( points[i]->location_.data_block() );
    //covariances->InsertNextTupleValue( points[i]->covariance_.data_block() );

    //vnl_matrix_fixed< double, dim, dim >  outer_product;
    //outer_product = vnl_inverse( points[i]->covariance_ );

    ////// eigenvectors and eigenvalues of the covariance matrix
    ////vnl_matrix< double >  V( dim, dim );
    ////vnl_vector< double >  D( dim );
    ////D.squared_magnitude();

    ////vnl_symmetric_eigensystem_compute( outer_product, V, D );

    //double strength = vnl_trace( outer_product );
    ////vcl_cout << "strength: " << strength << vcl_endl;
    //strengths->InsertNextValue( strength );
    strengths->InsertNextValue( points[i]->strength_ );
  }

  if( !poly_data ) poly_data = vtkSmartPointer< vtkPolyData >::New();
  else poly_data->Reset();

  poly_data->SetPoints( vtkpoints );
  poly_data->GetPointData()->SetScalars( strengths );
  //poly_data->GetPointData()->SetTensors( covariances );
}



#define CDCL_UTILS_VTK_INSTANTIATE_DOF( dim, dof )                                              \
template                                                                                        \
void                                                                                            \
cdcl_write_matches_VTK( vcl_vector< vbl_smart_ptr< cdcl_feature< dim > > > const &  moving,     \
                        vcl_vector< vbl_smart_ptr< cdcl_feature< dim > > > const &  fixed,      \
                        vcl_vector< vbl_smart_ptr< cdcl_match< dim > > >   const &  matches,    \
                        vbl_smart_ptr< cdcl_trans< dim, dof > >            const &  trans,      \
                        unsigned int                                                iteration );

#define CDCL_UTILS_VTK_INSTANTIATE( dim )                                                       \
template                                                                                             \
void                                                                                                 \
cdcl_features_to_poly_data( vcl_vector< vbl_smart_ptr< cdcl_feature< dim > > > const &  points,      \
                            vtkSmartPointer< vtkPolyData >                           &  poly_data ); \
template                                                                                             \
void                                                                                                 \
cdcl_features_to_poly_data( vcl_vector< vbl_smart_ptr< cdcl_feature_with_shape< dim > > > const &  points,      \
                            vtkSmartPointer< vtkPolyData >                                      &  poly_data ); \
template                                                                                             \
void                                                                                                 \
cdcl_features_to_poly_data( vcl_vector< vbl_smart_ptr< cdcl_feature_ICP< dim > > > const &  points,      \
                            vtkSmartPointer< vtkPolyData >                                      &  poly_data ); \
template                                                                                             \
bool                                                                                                 \
cdcl_read_keypoint_descriptors_VTK( vcl_vector< vbl_smart_ptr< cdcl_keypoint< dim > > >       &  keypoints,     \
                                    vcl_vector< vnl_vector< float > >                         &  descriptors,    \
                                    vcl_string                                          const    filename );    \
template                                                                                             \
void                                                                                                 \
cdcl_write_keypoint_descriptors_VTK( vcl_vector< vbl_smart_ptr< cdcl_keypoint< dim > > >  const &  keypoints,    \
                          vcl_vector< vnl_vector< float > >                  const &  descriptors,   \
                          vcl_string                                          const    filename );



// those utils that have specialization for 3d, cannot be instantiated with dim == 3
// so define a separate macro and use it only for 2d instantiations
#define CDCL_UTILS_IO_VTK_INSTANTIATE( dim )                                                         \
template                                                                                             \
void                                                                                                        \
cdcl_poly_data_to_keypoints_VTK( vtkSmartPointer< vtkPolyData >                       const &  poly_data,   \
                                 vcl_vector< vbl_smart_ptr< cdcl_keypoint< dim > > >        &  keypoints ); \
//template                                                                                                    \
//void                                                                                                        \
//cdcl_keypoints_to_poly_data_VTK( vcl_vector< vbl_smart_ptr< cdcl_keypoint< dim > > >  const &  keypoints,   \
//                                 vtkSmartPointer< vtkPolyData >                             &  poly_data );
//

// include instantiations here because we have template specialization
// the specialization is treated in the same way as a non template and would cause multiply defined functions
// (warning on VCC, error on gcc)
//CDCL_UTILS_VTK_INSTANTIATE_DOF( 2, 4 );
//CDCL_UTILS_VTK_INSTANTIATE_DOF( 2, 6 );
//
//CDCL_UTILS_IO_VTK_INSTANTIATE( 2 );
//
//CDCL_UTILS_VTK_INSTANTIATE_DOF( 3, 6 );
//CDCL_UTILS_VTK_INSTANTIATE_DOF( 3, 12 );
//
//CDCL_UTILS_VTK_INSTANTIATE( 2 );
//CDCL_UTILS_VTK_INSTANTIATE( 3 );



#endif
