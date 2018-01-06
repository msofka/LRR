#include <vcl_iomanip.h>
#include <vcl_sstream.h>
#include <vcl_algorithm.h>
#include <vgui/vgui.h>
#include <vgui/vgui_find.h>
#include <vgui/vgui_utils.h>
#include <vgui/vgui_viewer2D_tableau.h>
#include <vgui/vgui_viewer2D_tableau_sptr.h>
#include <vnl/algo/vnl_symmetric_eigensystem.h>

#include "cdcl_display.h"

//:
// \file
// \brief  Displaying points and matches with a display callback function.
//         The dependence of the estimation classes on vgui is avoided
//         by giving a choice of supplying callback or not.
//         In case vgui is not available, estimation can still run
//         without displying progress.
//         Function finalize_display must be called implicitly.
//         The main functionality is to let cdcl_estimation run display_points
//         and once this is done for each region, run finalize_display at the end.
//         The 2d and 3d version cannot be united with a single parent
//         because vgui_easy2d and vgui_easy3d do not have a single parent.
// \author Michal Sofka
// \date   May 2006

#define POINTS 1
#define CORRESPONDENCES 1
#define COVARIANCES 1

#define DARK_BACKGROUND 0

template < unsigned int dim, unsigned int dof >
void
cdcl_display_callback2d< dim, dof >::display_points( vcl_vector< feature_sptr_type > const &  moving,
                                                     vcl_vector< feature_sptr_type > const &  fixed,
                                                     vcl_vector< match_sptr_type >   const &  matches,
                                                     trans_sptr_type                 const &  trans,
                                                     unsigned int                             iteration )
{
  // compute minimum and maximum x and y coordinates
  double min_x = fixed[0]->location_[0];
  double max_x = fixed[0]->location_[0];
  double min_y = fixed[0]->location_[1];
  double max_y = fixed[0]->location_[1];
  for( typename vcl_vector< match_sptr_type >::size_type  i = 0; i < fixed.size(); ++i ) {
    min_x = vcl_min( min_x, fixed[i]->location_[0] );
    min_y = vcl_min( min_y, fixed[i]->location_[1] );
    max_x = vcl_max( max_x, fixed[i]->location_[0] );
    max_y = vcl_max( max_y, fixed[i]->location_[1] );
  }

  double width  = max_x - min_x;
  double height = max_y - min_y;

  vcl_cout << "scale_for_disp: " << (0.5*vcl_sqrt(width*width + height*height)) << vcl_endl;
  double scale_for_disp = 10;//  300.0 / (0.5*vcl_sqrt(width*width + height*height));
  //double scale_for_disp = 300.0 / (0.5*vcl_sqrt(width*width + height*height));
//scale_for_disp = 800.0;
  this->res_ = scale_for_disp;

  // Jacobian w.r.t. p (moving point)
  vnl_matrix_fixed< double, dim, dim > const &  Jp = trans->jacobian_wrt_loc();


  double radius = 6.0;

  // go through all moving points, transform and draw them
  for( typename vcl_vector< feature_sptr_type >::size_type  i = 0; i < moving.size(); ++i ) {
    // moving point p
    feature_sptr_type  p = moving[i];
    vnl_vector_fixed< double, dim >  loc = trans->map_loc( p->location_ );

    // draw moving point
    #if DARK_BACKGROUND
      easy2d_->set_foreground( 1.0, 1.0, 0.0 );  // dark background
    #else
      easy2d_->set_foreground( 1.0, 0.0, 0.0 );  // bright background
    #endif
    
    easy2d_->set_line_width( 4.0 );
    easy2d_->add_circle( loc[0] * res_,
                         loc[1] * res_,
                         radius );
  }

  // go through all fixed points and draw them
  for( typename vcl_vector< feature_sptr_type >::size_type  j = 0; j < fixed.size(); ++j ) {
    // fixed point q
    feature_sptr_type  q = fixed[j];

    // draw fixed point
    #if DARK_BACKGROUND
      easy2d_->set_foreground( 1.0, 0.0, 0.0 );  // dark background
    #else
      easy2d_->set_foreground( 0.0, 0.0, 1.0 );  // bright background
    #endif

    easy2d_->set_line_width( 4.0 );
    double half_line = radius;
    float x0 = q->location_[0] * res_ + half_line;
    float y0 = q->location_[1] * res_ + half_line;
    float x1 = q->location_[0] * res_ - half_line;
    float y1 = q->location_[1] * res_ - half_line;
    easy2d_->add_line( x0, y0, x1, y1 );
    x0 = q->location_[0] * res_ + half_line;
    y0 = q->location_[1] * res_ - half_line;
    x1 = q->location_[0] * res_ - half_line;
    y1 = q->location_[1] * res_ + half_line;
    easy2d_->add_line( x0, y0, x1, y1 );

    //easy2d_->add_circle( q->location_[0] * res_,
    //                     q->location_[1] * res_,
    //                     radius );
  }


  // go through all matches and draw correspondences
  double min_weight = 0.0;
  double max_weight = 1.0;
  for( typename vcl_vector< match_sptr_type >::size_type  i = 0; i < matches.size(); ++i ) {
    // moving point p
    feature_sptr_type  p = matches[i]->from_;
    vnl_vector_fixed< double, dim >  mapped_moving = trans->map_loc( p->location_ );

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

      const float color_coef = ( matches[i]->w_[j] - min_weight ) / ( max_weight - min_weight );

#if CORRESPONDENCES 
      // draw correspondence
      easy2d_->set_line_width( 2.0 );
      #if DARK_BACKGROUND
        easy2d_->set_foreground( 0.0, 1.0, 0.0 );  // dark background
      #else
        easy2d_->set_foreground( 0.6f, 0.0f, 1.0f-color_coef );  // bright background
      #endif
      easy2d_->add_line( mapped_moving[0] * res_, mapped_moving[1] * res_,
                         q->location_[0] * res_, q->location_[1] * res_ );
#endif


#if COVARIANCES

      // covariance of a correspondence match
      vnl_matrix_fixed< double, dim, dim >  Cij = this->transfer_covar( Jth, trans ) +
                                                  Jp * p->covariance_ * vnl_transpose( Jp ) +
                                                  q->covariance_;

      const vnl_matrix< double >  A = Cij;
      vnl_symmetric_eigensystem_compute( A, V, D );
    
      // select sqrt of eigenvalues as ellipse width and height
      // easy2d ellipse width and height correspond to minor and major axes
      float height = vcl_sqrt( D( 1 ) );
      float width  = vcl_sqrt( D( 0 ) );

      // compute ellipse orientation
      float phi = atan2( V( 0, 1 ), V( 1, 1 ) );

      vnl_vector_fixed< double, 2 >  mid = q->location_ + ( ( mapped_moving - q->location_ ) * 0.5 );

      // draw covariance ellipse
      easy2d_->set_line_width( 2.0 );
      #if DARK_BACKGROUND
        easy2d_->set_foreground( 0.0, 1.0, 1.0 );  // dark background
      #else
        easy2d_->set_foreground( 0.0f, 0.7f, 0.0f );  // bright background
      #endif

      easy2d_->add_ellipse( mid[0] * res_, mid[1] * res_, width * res_, height * res_, phi );



  #if 0
      // draw covariance cross (inner part is feature covariance, outer part is parameter covariance)
      // the lengths are determined by a projection on the eigenvectors computed from Cij
      //
      // covariance of a correspondence match
      vnl_matrix_fixed< double, dim, dim >  Cth = Jth * trans->get_covariance() * vnl_transpose( Jth );
      vnl_matrix_fixed< double, dim, dim >  Cpq = Jp * p->covariance_ * vnl_transpose( Jp ) + 
                                                  q->covariance_;

      // square root of eigenvalues
      vnl_matrix_fixed< double, dim, 1 >  column;
      column.set_column( 0, V.get_column(0) );
      vnl_matrix_fixed< double, 1, 1 > eth0 = vnl_transpose( column ) * Cth * column;
      vnl_matrix_fixed< double, 1, 1 > epq0 = vnl_transpose( column ) * Cpq * column;
      column.set_column( 0, V.get_column(1) );
      vnl_matrix_fixed< double, 1, 1 > eth1 = vnl_transpose( column ) * Cth * column;
      vnl_matrix_fixed< double, 1, 1 > epq1 = vnl_transpose( column ) * Cpq * column;

      double ev0 = vcl_sqrt( epq0( 0, 0 ) );
      double ev1 = vcl_sqrt( epq1( 0, 0 ) );

      // end points of covariance crosses
      vnl_vector_fixed< double, dim >  P1 = mid * res_ + ev0 * V.get_column(0) * res_;
      vnl_vector_fixed< double, dim >  P2 = mid * res_ - ev0 * V.get_column(0) * res_;
      vnl_vector_fixed< double, dim >  Q1 = mid * res_ + ev1 * V.get_column(1) * res_;
      vnl_vector_fixed< double, dim >  Q2 = mid * res_ - ev1 * V.get_column(1) * res_;

      // draw feature covariance cross
      easy2d_->set_foreground( 0.0, 1.0, 1.0 );
      easy2d_->add_line( P1[0], P1[1], P2[0], P2[1] );
      easy2d_->add_line( Q1[0], Q1[1], Q2[0], Q2[1] );

      ev0 = vcl_sqrt( eth0( 0, 0 ) );
      ev1 = vcl_sqrt( eth1( 0, 0 ) );

      vnl_vector_fixed< double, dim >  U1 = P1 + ev0 * V.get_column(0) * res_;
      vnl_vector_fixed< double, dim >  U2 = P2 - ev0 * V.get_column(0) * res_;
      vnl_vector_fixed< double, dim >  V1 = Q1 + ev1 * V.get_column(1) * res_;
      vnl_vector_fixed< double, dim >  V2 = Q2 - ev1 * V.get_column(1) * res_;

      // draw parameter covariance cross
      easy2d_->set_foreground( 1.0, 0.0, 1.0 );
      easy2d_->add_line( P1[0], P1[1], U1[0], U1[1] );
      easy2d_->add_line( P2[0], P2[1], U2[0], U2[1] );
      easy2d_->add_line( Q1[0], Q1[1], V1[0], V1[1] );
      easy2d_->add_line( Q2[0], Q2[1], V2[0], V2[1] );
  #endif



      //// covariance of the moving point location
      //vnl_symmetric_eigensystem_compute( Jp * p->covariance_ * vnl_transpose( Jp ), V, D );
    
      //// select sqrt of eigenvalues as ellipse width and height
      //height = 2*vcl_sqrt( D( 1 ) );
      //width  = 2*vcl_sqrt( D( 0 ) );

      //// compute ellipse orientation
      //phi = atan2( V( 0, 1 ), V( 1, 1 ) );

      //mid = mapped_moving;

      //// draw covariance ellipse
      //easy2d_->set_foreground( 0.0, 1.0, 0.0 );
      //easy2d_->add_ellipse( mid[0] * res_, mid[1] * res_, width * res_, height * res_, phi );

#endif
    }
  }

}


template < unsigned int dim, unsigned int dof >
void
cdcl_display_callback2d< dim, dof >::finalize_display( unsigned int  iteration )
{
  vcl_ostringstream  matches_name;
  matches_name << "iteration" << vcl_setw( 4 ) << vcl_setfill( '0' ) << iteration;
  //easy2d_->print_psfile( matches_name.str() + ".ps", 1, true, 0, 0 );
  //vcl_cout << "Printed " << matches_name.str() << vcl_endl;

  // center on the current fixed region
  vgui_viewer2D_tableau_sptr  viewer;
  viewer.vertical_cast( vgui_find_above_by_type_name( easy2d_, "vgui_viewer2D_tableau" ) );
  
  // use for cdc_refinement
  //viewer->token.scaleX =  3.0;
  //viewer->token.scaleY =  3.0;
  //viewer->token.offsetX = viewer->token.scaleX*(-min_x+20);
  //viewer->token.offsetY = viewer->token.scaleY*(-min_y+20);

  // use for estimate and convergence on H
  viewer->token.scaleX =  0.2f;//1.0;//scale_for_disp;
  viewer->token.scaleY =  0.2f;//1.0;//scale_for_disp;
  //viewer->token.scaleX =  1.0f;///scale_for_disp;
  //viewer->token.scaleY =  1.0f;///scale_for_disp;
  viewer->token.offsetX = 0.0f;
  viewer->token.offsetY = 0.0f;
//viewer->token.offsetX = 300.0f;
//viewer->token.offsetY = 300.0f;
  //viewer->token.offsetX = viewer->token.scaleX*(-min_x*scale_for_disp + 50);
  //viewer->token.offsetY = viewer->token.scaleY*(-min_y*scale_for_disp + 50);

  viewer->post_redraw();

  easy2d_->post_redraw();
  vgui::run_till_idle();
  //vgui::run();
  vcl_string  buffer_fname = matches_name.str() + ".png";
  vgui_utils::dump_colour_buffer( buffer_fname.c_str() );

  easy2d_->clear();
}


template < unsigned int dim, unsigned int dof >
void
cdcl_display_callback3d< dim, dof >::display_points( vcl_vector< feature_sptr_type > const &  moving,
                                                     vcl_vector< feature_sptr_type > const &  fixed,
                                                     vcl_vector< match_sptr_type >   const &  matches,
                                                     trans_sptr_type                 const &  trans,
                                                     unsigned int                             iteration )
{
  // Jacobian w.r.t. p (moving point)
  vnl_matrix_fixed< double, dim, dim > const &  Jp = trans->jacobian_wrt_loc();

  double res = 200;

  // go through all matches and draw correspondences
  double min_weight = 0.0;
  double max_weight = 1.0;
  for( typename vcl_vector< match_sptr_type >::size_type  i = 0; i < matches.size(); ++i ) {
    // moving point p
    feature_sptr_type  p = matches[i]->from_;
    vnl_vector_fixed< double, dim >  mapped_moving = trans->map_loc( p->location_ );

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

      const double color_coef = ( matches[i]->w_[j] - min_weight ) / ( max_weight - min_weight );

#if CORRESPONDENCES
      // draw correspondence
      easy3d_->set_foreground( 0.5*color_coef, 0.0, 0.5-0.5*color_coef );
      easy3d_->add_line( mapped_moving[0] * res, mapped_moving[1] * res, mapped_moving[2] * res,
                         q->location_[0] * res, q->location_[1] * res, q->location_[2] * res );
#endif

#if COVARIANCES
      // covariance of a correspondence match
      vnl_matrix_fixed< double, dim, dim >  Cth = Jth * trans->get_covariance() * vnl_transpose( Jth );
      vnl_matrix_fixed< double, dim, dim >  Cpq = Jp * p->covariance_ * vnl_transpose( Jp ) + 
                                                  q->covariance_;
      vnl_matrix_fixed< double, dim, dim >  Cij = Cth + Cpq;
                                                  

      vnl_vector_fixed< double, dim >  mid = q->location_ + ( ( mapped_moving - q->location_ ) * 0.5 );

  #if 0
      // draw covariance cross (single cross for each covariance)
      
      // eigen decomposition of the covariance matrix
      const vnl_matrix< double >  A = Cij;
      vnl_symmetric_eigensystem_compute( A, V, D );
    
      // square root of eigenvalues
      double e0 = vcl_sqrt( D(0) );
      double e1 = vcl_sqrt( D(1) );
      double e2 = vcl_sqrt( D(2) );

      // end points of covariance crosses
      vnl_vector_fixed< double, dim >  P1 = mid * res + e0 * V.get_column(0) * res;
      vnl_vector_fixed< double, dim >  P2 = mid * res - e0 * V.get_column(0) * res;
      vnl_vector_fixed< double, dim >  Q1 = mid * res + e1 * V.get_column(1) * res;
      vnl_vector_fixed< double, dim >  Q2 = mid * res - e1 * V.get_column(1) * res;
      vnl_vector_fixed< double, dim >  R1 = mid * res + e2 * V.get_column(2) * res;
      vnl_vector_fixed< double, dim >  R2 = mid * res - e2 * V.get_column(2) * res;

      // draw covariance cross
      easy3d_->set_foreground( 1.0, 1.0, 1.0 );
      easy3d_->add_line( P1[0], P1[1], P1[2], P2[0], P2[1], P2[2] );
      easy3d_->add_line( Q1[0], Q1[1], Q1[2], Q2[0], Q2[1], Q2[2] );
      easy3d_->add_line( R1[0], R1[1], R1[2], R2[0], R2[1], R2[2] );

  #else

      // draw covariance cross (inner part is feature covariance, outer part is parameter covariance)
      // the lengths are determined by a projection on the eigenvectors computed from Cij
      //
      // eigen decomposition of the covariance matrix
      const vnl_matrix< double >  A = Cij;
      vnl_symmetric_eigensystem_compute( A, V, D );
    
      // square root of eigenvalues
      vnl_matrix_fixed< double, dim, 1 >  column;
      column.set_column( 0, V.get_column(0) );
      vnl_matrix_fixed< double, 1, 1 > eth0 = vnl_transpose( column ) * Cth * column;
      vnl_matrix_fixed< double, 1, 1 > epq0 = vnl_transpose( column ) * Cpq * column;
      column.set_column( 0, V.get_column(1) );
      vnl_matrix_fixed< double, 1, 1 > eth1 = vnl_transpose( column ) * Cth * column;
      vnl_matrix_fixed< double, 1, 1 > epq1 = vnl_transpose( column ) * Cpq * column;
      column.set_column( 0, V.get_column(2) );
      vnl_matrix_fixed< double, 1, 1 > eth2 = vnl_transpose( column ) * Cth * column;
      vnl_matrix_fixed< double, 1, 1 > epq2 = vnl_transpose( column ) * Cpq * column;

      double ev0 = vcl_sqrt( epq0( 0, 0 ) );
      double ev1 = vcl_sqrt( epq1( 0, 0 ) );
      double ev2 = vcl_sqrt( epq2( 0, 0 ) );

      // end points of covariance crosses
      vnl_vector_fixed< double, dim >  P1 = mid * res + ev0 * V.get_column(0) * res;
      vnl_vector_fixed< double, dim >  P2 = mid * res - ev0 * V.get_column(0) * res;
      vnl_vector_fixed< double, dim >  Q1 = mid * res + ev1 * V.get_column(1) * res;
      vnl_vector_fixed< double, dim >  Q2 = mid * res - ev1 * V.get_column(1) * res;
      vnl_vector_fixed< double, dim >  R1 = mid * res + ev2 * V.get_column(2) * res;
      vnl_vector_fixed< double, dim >  R2 = mid * res - ev2 * V.get_column(2) * res;

      // draw feature covariance cross
      easy3d_->set_foreground( 0.0, 1.0, 1.0 );
      easy3d_->add_line( P1[0], P1[1], P1[2], P2[0], P2[1], P2[2] );
      easy3d_->add_line( Q1[0], Q1[1], Q1[2], Q2[0], Q2[1], Q2[2] );
      easy3d_->add_line( R1[0], R1[1], R1[2], R2[0], R2[1], R2[2] );

      ev0 = vcl_sqrt( eth0( 0, 0 ) );
      ev1 = vcl_sqrt( eth1( 0, 0 ) );
      ev2 = vcl_sqrt( eth2( 0, 0 ) );

      vnl_vector_fixed< double, dim >  U1 = P1 + ev0 * V.get_column(0) * res;
      vnl_vector_fixed< double, dim >  U2 = P2 - ev0 * V.get_column(0) * res;
      vnl_vector_fixed< double, dim >  V1 = Q1 + ev1 * V.get_column(1) * res;
      vnl_vector_fixed< double, dim >  V2 = Q2 - ev1 * V.get_column(1) * res;
      vnl_vector_fixed< double, dim >  W1 = R1 + ev2 * V.get_column(2) * res;
      vnl_vector_fixed< double, dim >  W2 = R2 - ev2 * V.get_column(2) * res;

      // draw parameter covariance cross
      easy3d_->set_foreground( 1.0, 1.0, 0.0 );
      easy3d_->add_line( P1[0], P1[1], P1[2], U1[0], U1[1], U1[2] );
      easy3d_->add_line( P2[0], P2[1], P2[2], U2[0], U2[1], U2[2] );
      easy3d_->add_line( Q1[0], Q1[1], Q1[2], V1[0], V1[1], V1[2] );
      easy3d_->add_line( Q2[0], Q2[1], Q2[2], V2[0], V2[1], V2[2] );
      easy3d_->add_line( R1[0], R1[1], R1[2], W1[0], W1[1], W1[2] );
      easy3d_->add_line( R2[0], R2[1], R2[2], W2[0], W2[1], W2[2] );

  #endif

#endif

    }
  }


  double radius = 5.0;

#if POINTS
  // go through all moving points, transform and draw them
  for( typename vcl_vector< feature_sptr_type >::size_type  i = 0; i < moving.size(); ++i ) {
    // moving point p
    feature_sptr_type  p = moving[i];
    vnl_vector_fixed< double, dim >  loc = trans->map_loc( p->location_ );

    // draw moving point
    easy3d_->set_foreground( 0.0, 1.0, 0.0 );
    easy3d_->add_point( loc[0] * res,
                        loc[1] * res,
                        loc[2] * res );
  }

  // go through all fixed points and draw them
  for( typename vcl_vector< feature_sptr_type >::size_type  j = 0; j < fixed.size(); ++j ) {
    // fixed point q
    feature_sptr_type  q = fixed[j];

    // draw fixed point
    easy3d_->set_foreground( 1.0, 0.0, 1.0 );
    easy3d_->add_point( q->location_[0] * res,
                        q->location_[1] * res,
                        q->location_[2] * res );
  }
  easy3d_->set_point_radius( radius );
#endif

}


template < unsigned int dim, unsigned int dof >
void
cdcl_display_callback3d< dim, dof >::finalize_display( unsigned int  iteration )
{
  vcl_ostringstream  matches_name;
  matches_name << "iteration" << vcl_setw( 3 ) << vcl_setfill( '0' ) << iteration;

  easy3d_->post_redraw();
  //vgui::run();
  vgui::run_till_idle();
  // examine if some flag being set to indicate it should stop

  vcl_string  buffer_fname = matches_name.str() + ".png";
  vgui_utils::dump_colour_buffer( buffer_fname.c_str() );

  easy3d_->clear();
}


// Instantiations are here.
#define CDCL_DISPLAY_CALLBACK2D_INSTANTIATE( dim, dof )  \
template                                                 \
class cdcl_display_callback2d< dim, dof >;

#define CDCL_DISPLAY_CALLBACK3D_INSTANTIATE( dim, dof )  \
template                                                 \
class cdcl_display_callback3d< dim, dof >;

#define CDCL_DISPLAY_CALLBACK2D_TRANSFER_INSTANTIATE( dim, dof )  \
template                                                          \
class cdcl_display_callback2d_transfer< dim, dof >;

#define CDCL_DISPLAY_CALLBACK3D_TRANSFER_INSTANTIATE( dim, dof )  \
template                                                          \
class cdcl_display_callback3d_transfer< dim, dof >;
