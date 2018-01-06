#include <vcl_iomanip.h>
#include <vcl_sstream.h>
#include <vcl_algorithm.h>
#include <vnl/algo/vnl_symmetric_eigensystem.h>

#include "cdcl_extract_data.h"
#include <cdcl/cdcl_utils_VTK.h>

#include "vtkImageData.h"
#include "vtkSmartPointer.h"
#include "vtkOpenGLRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkOpenGLRenderWindow.h"
#include "vtkTensorGlyphScaled.h"
#include "vtkWindowToImageFilter.h"
#include "vtkActor.h"
#include "vtkPolyDataMapper.h"
#include "vtkJPEGWriter.h"
#include "vtkGlyphSource2D.h"
#include "vtkProperty.h"
#include "vtkGraphicsFactory.h"
#include "vtkImagingFactory.h"
#include "vtkCamera.h"
#include "vtkPNGWriter.h"
#include "vtkImageActor.h"
#include "vtkWindowLevelLookupTable.h"
#include "vtkImageMapToColors.h"
#include "vtkRenderWindowInteractor.h"

#include "vtkCylinderSource.h"
#include "vtkGlyph3D.h"
#include "vtkPNGReader.h"
#include "vtkSphereSource.h"
#include "vtkExtractGeometry.h"
#include "vtkExtractPolyDataGeometry.h"
#include "vtkBox.h"
#include "vtkImageViewer2.h"
#include "vtkPolyDataWriter.h"
#include "vtkXMLPolyDataWriter.h"

#include "itkImage.h"
#include "itkImageFileWriter.h"
#include "itkRescaleIntensityImageFilter.h"


//#include "vtkMesaRenderer.h"
//#include "vtkMesaRenderWindow.h"
//#include "vtkMesaActor.h"
//#include "vtkMesaPolyDataMapper.h" 


//:
// \file
// \brief  Displaying points and matches with a display callback function.
//         The dependence of the estimation classes on the display is avoided
//         by giving a choice of supplying callback or not.
//         The other advantage is that data members do not neeed to be exposed
//         in order to be displayed.
//         In case callback is not available, estimation can still run
//         without displying progress.
//         Previously, callback generated display (when dependent on vgui)
//         but now it only prepares values (poly data) to render with VTK.
// \author Michal Sofka
// \date   Aug 2007


namespace {

// Convert cdcl features to VTK polydata for display purposes
template < unsigned int dim, class match_type >
void
cdcl_matches_to_poly_data( vcl_vector< vbl_smart_ptr< match_type > > const &  matches,
                           vtkSmartPointer< vtkPolyData >                  &  poly_data )
{
  typedef typename match_type::feature_sptr_type  feature_sptr_type;
  typedef typename match_type::sptr               match_sptr_type;

  vtkSmartPointer< vtkPoints >       vtkpoints   = vtkSmartPointer< vtkPoints >::New();
  vtkSmartPointer< vtkDoubleArray >  weights     = vtkSmartPointer< vtkDoubleArray >::New();

  weights->SetName( "weights" );
  weights->SetNumberOfComponents( 1 );

  // write data points
  //
  for( typename vcl_vector< match_sptr_type >::size_type  i = 0; i < matches.size(); ++i ) {

    for( typename vcl_vector< match_sptr_type >::size_type  j = 0; j < matches[i]->to_.size(); ++j ) {

      // fixed point q
      feature_sptr_type  q = matches[i]->to_[j];

      vtkpoints->InsertNextPoint( q->location_.data_block() );

      weights->InsertNextValue( matches[i]->w_[j] );

    }
  }

  if( !poly_data ) poly_data = vtkSmartPointer< vtkPolyData >::New();
  else poly_data->Reset();

  poly_data->SetPoints( vtkpoints );
  poly_data->GetPointData()->SetActiveScalars( "weights" );
  poly_data->GetPointData()->SetScalars( weights );
}


template < unsigned int dim, class match_type >
void
save_matches_as_polydata( vcl_vector< vbl_smart_ptr< match_type > > const &  matches,
                          vcl_string  file_name )
{
  vtkSmartPointer< vtkPolyData >  matches_poly_data = vtkSmartPointer< vtkPolyData >::New();

  cdcl_matches_to_poly_data< dim, match_type >( matches, matches_poly_data );

  // write the features
  vtkSmartPointer< vtkXMLPolyDataWriter >  featureWriter = vtkSmartPointer< vtkXMLPolyDataWriter >::New();

  featureWriter->SetInput( matches_poly_data );
  featureWriter->SetFileName( file_name.c_str() );
  featureWriter->Write();
}

}


template < unsigned int dim, unsigned int dof, class feature_typeT >
cdcl_extract_data_callback3d< dim, dof, feature_typeT >::cdcl_extract_data_callback3d()
{
  m_PolyDataCovariances = vtkSmartPointer< vtkPolyData >::New();
  m_PolyDataMovingROIMapped = vtkSmartPointer< vtkPolyData >::New();

  //if( !PolyDataCovariances ) PolyDataCovariances = vtkSmartPointer< vtkPolyData >::New();
  //else PolyDataCovariances->Reset();

  //m_PolyDataCovariances = PolyDataCovariances;


  //if( !PolyDataMovingROIMapped ) PolyDataMovingROIMapped = vtkSmartPointer< vtkPolyData >::New();
  //else PolyDataMovingROIMapped->Reset();

  //m_PolyDataMovingROIMapped = PolyDataMovingROIMapped;

}


template < unsigned int dim, unsigned int dof, class feature_typeT >
void
cdcl_extract_data_callback3d< dim, dof, feature_typeT >::initialize_moving_mapped_features( vcl_vector< feature_sptr_type > const &  moving,
                                                                                            trans_sptr_type                 const &  trans )
{
  // Jacobian w.r.t. p (moving point)
  vnl_matrix_fixed< double, dim, dim > const &  Jp = trans->jacobian_wrt_loc();

  vcl_vector< feature_sptr_type >  one_level_mapped;
  // go through all matches
  for( typename vcl_vector< feature_sptr_type >::size_type  i = 0; i < moving.size(); ++i ) {
    // moving point p
    feature_sptr_type  p = moving[i];
    vnl_vector_fixed< double, dim >  mapped_moving = trans->map_loc( p->location_ );

#define HAS_COVARIANCE 0

#if HAS_COVARIANCE
    vnl_matrix_fixed< double, dim, dim >  mapped_covariance = Jp * p->covariance_ * Jp.transpose();
    one_level_mapped.push_back( new feature_type( mapped_moving, mapped_covariance ) );
#else
    feature_sptr_type  new_feature = new feature_type;
    new_feature->location_ = mapped_moving;
    one_level_mapped.push_back( new_feature );
#endif
  }

  cdcl_features_to_poly_data( one_level_mapped, m_PolyDataMovingROIMapped );

}



template < unsigned int dim, unsigned int dof, class feature_typeT >
void
cdcl_extract_data_callback3d< dim, dof, feature_typeT >::display_points( vcl_vector< feature_sptr_type > const &  moving,
                                                          vcl_vector< feature_sptr_type > const &  fixed,
                                                          vcl_vector< match_sptr_type >   const &  matches,
                                                          trans_sptr_type                 const &  trans,
                                                          unsigned int                             iteration )
{

  // feature set for creating covariances polydata
  vcl_vector< feature_sptr_type >  covariances;

  // Jacobian w.r.t. p (moving point)
  vnl_matrix_fixed< double, dim, dim > const &  Jp = trans->jacobian_wrt_loc();

  vcl_vector< feature_sptr_type >  one_level_mapped;
  // go through all matches
  for( typename vcl_vector< match_sptr_type >::size_type  i = 0; i < matches.size(); ++i ) {
    // moving point p
    feature_sptr_type  p = matches[i]->from_;
    vnl_vector_fixed< double, dim >  mapped_moving = trans->map_loc( p->location_ );

#if !HAS_COVARIANCE
    feature_sptr_type  new_feature = new feature_type;
    new_feature->location_ = mapped_moving;
    one_level_mapped.push_back( new_feature );

#else

    // Jacobian w.r.t. transformation parameters
    vnl_matrix_fixed< double, dim, dof > const &  Jth = trans->jacobian_wrt_par( p->location_ );

    double max_weight = *vcl_max_element( matches[i]->w_.begin(), matches[i]->w_.end() );

    // eigenvectors and eigenvalues of the covariance matrix
    vnl_matrix< double >  V( dim, dim );
    vnl_vector< double >  D( dim );

    // go through all matches (fixed points)
    for( typename vcl_vector< match_sptr_type >::size_type  j = 0; j < matches[i]->to_.size(); ++j ) {
     if( matches[i]->w_[j] != max_weight ) continue;

      // fixed point q
      feature_sptr_type  q = matches[i]->to_[j];


      // covariance of a correspondence match
      vnl_matrix_fixed< double, dim, dim >  Cij = this->transfer_covar( Jth, trans ) +
                                                  mapped_covariance +
                                                  q->covariance_;

      vnl_vector_fixed< double, dim >  mid = q->location_ + ( ( mapped_moving - q->location_ ) * 0.5 );

      covariances.push_back( new feature_type( mid, Cij ) );
    }
#endif
  }

#if HAS_COVARIANCE
  cdcl_features_to_poly_data( covariances, m_PolyDataCovariances );
#endif

  cdcl_features_to_poly_data( one_level_mapped, m_PolyDataMovingROIMapped );

//vcl_ostringstream  matchesName;
//matchesName << "matches" << vcl_setw( 6 ) << vcl_setfill( '0' ) << iteration << ".vtk";
//
//save_matches_as_polydata< dim, cdcl_match< dim, feature_type > >( matches, matchesName.str() );
}


template < unsigned int dim, unsigned int dof, class feature_typeT >
void
cdcl_extract_data_callback3d< dim, dof, feature_typeT >::finalize_display( unsigned int  iteration )
{
}



// Instantiations are here.
#define CDCL_EXTRACT_DATA_CALLBACK3D_INSTANTIATE( dim, dof, feature_typeT )  \
template                                                 \
class cdcl_extract_data_callback3d< dim, dof, feature_typeT >;

