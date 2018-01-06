#ifndef cdcl_display_h_
#define cdcl_display_h_

#include <cdcl/cdcl_feature.h>
#include <cdcl/cdcl_match.h>
#include <cdcl/cdcl_trans.h>
#include <cdcl/cdcl_estimation_transfer.h>

#include <vnl/vnl_transpose.h>

#include <vgui/vgui_easy2D_tableau.h>
#include <vgui/vgui_easy3D_tableau.h>

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


// Callback which displays points is passed to the registration object.
// This way the registration object does not depend on vgui.
template < unsigned int dim, unsigned int dof >
class cdcl_display_callback2d
{
public:
  cdcl_display_callback2d( vgui_easy2D_tableau_sptr &  easy2d )
    : easy2d_( easy2d ) {};

  typedef typename cdcl_feature< dim >::sptr      feature_sptr_type;
  typedef typename cdcl_match< dim >::sptr        match_sptr_type;
  typedef typename cdcl_trans< dim, dof >::sptr   trans_sptr_type;

  // Display 2d points.
  virtual void display_points( vcl_vector< feature_sptr_type > const &  moving,
                               vcl_vector< feature_sptr_type > const &  fixed,
                               vcl_vector< match_sptr_type >   const &  matches,
                               trans_sptr_type                 const &  trans,
                               unsigned int                             iteration );

  // Callback which calls member display_points to draw points on tableau.
  // It is done this way because 1) we need to know the caller, 2) callback must be static.
  static void display_points( void * caller,   
                              vcl_vector< feature_sptr_type > const &  moving,
                              vcl_vector< feature_sptr_type > const &  fixed,
                              vcl_vector< match_sptr_type >   const &  matches,
                              trans_sptr_type                 const &  transform,
                              unsigned int                             iteration )
  {
    cdcl_display_callback2d* mySelf = (cdcl_display_callback2d*)  caller;
    mySelf->display_points( moving, fixed, matches, transform, iteration );
    mySelf->finalize_display( iteration );
  }

  // Finalize display callback tableau (post redraw and save).
  virtual void finalize_display( unsigned int  iteration );

protected:
  // Transfer error covariance.
  inline virtual vnl_matrix_fixed< double, dim, dim >  transfer_covar( vnl_matrix_fixed< double, dim, dof > const &  Jth,
                                                                       trans_sptr_type                      const &  trans ) {
    return Jth * trans->get_covariance() * vnl_transpose( Jth );
  }

  vgui_easy2D_tableau_sptr  easy2d_;
    
  // Resolution. Scaling factor by which all points will be multiplied.
  double res_;
};


// Callback which displays points is passed to the registration object.
// This way the registration object does not depend on vgui.
// Transfer error covariances are used.
template < unsigned int dim, unsigned int dof >
class cdcl_display_callback2d_transfer : public cdcl_display_callback2d< dim, dof >
{
public:
  cdcl_display_callback2d_transfer( vgui_easy2D_tableau_sptr &  easy2d )
    : cdcl_display_callback2d< dim, dof >( easy2d ) {};

  typedef cdcl_display_callback2d< dim, dof >             superclass_type;
  typedef typename superclass_type::feature_sptr_type                      feature_sptr_type;
  typedef typename superclass_type::match_sptr_type                        match_sptr_type;
  typedef typename superclass_type::trans_sptr_type                        trans_sptr_type;
  //typedef cdcl_estimation_transfer_piecewise<dim, dof>                     estimation_type;
  typedef cdcl_estimation_transfer<dim, dof>                               estimation_type;
  //typedef typename cdcl_estimation_ICP<dim, dof> >                         estimation_type;
  typedef typename estimation_type::sptr                                   estimation_sptr_type;

  // Display 2d points.
  virtual void display_points( vcl_vector< feature_sptr_type > const &  moving,
                               vcl_vector< feature_sptr_type > const &  fixed,
                               vcl_vector< match_sptr_type >   const &  matches,
                               trans_sptr_type                 const &  trans,
                               unsigned int                             iteration )
  { superclass_type::display_points( moving, fixed, matches, trans, iteration ); };

  // Callback which calls member display_points to draw points on tableau.
  // It is done this way because 1) we need to know the caller, 2) callback must be static.
  static void display_points( void * caller,   
                              vcl_vector< feature_sptr_type > const &  moving,
                              vcl_vector< feature_sptr_type > const &  fixed,
                              vcl_vector< match_sptr_type >   const &  matches,
                              trans_sptr_type                 const &  transform,
                              unsigned int                             iteration )
  {
    cdcl_display_callback2d_transfer* mySelf = (cdcl_display_callback2d_transfer*)  caller;
    mySelf->display_points( moving, fixed, matches, transform, iteration );
    mySelf->finalize_display( iteration );
  }

protected:
  // Transfer error covariance.
  inline virtual vnl_matrix_fixed< double, dim, dim >  transfer_covar( vnl_matrix_fixed< double, dim, dof > const &  Jth,
                                                                       trans_sptr_type                      const &  trans ) {
    return trans->get_covarianceJ();
  }

};


// Callback which displays points is passed to the registration object.
// This way the registration object does not depend on vgui.
template < unsigned int dim, unsigned int dof >
class cdcl_display_callback3d
{
public:
  cdcl_display_callback3d( vgui_easy3D_tableau_sptr &  easy3d )
    : easy3d_( easy3d ) {};

  typedef typename cdcl_feature< dim >::sptr      feature_sptr_type;
  typedef typename cdcl_match< dim >::sptr        match_sptr_type;
  typedef typename cdcl_trans< dim, dof >::sptr   trans_sptr_type;

  // Display 3d points.
  virtual void display_points( vcl_vector< feature_sptr_type > const &  moving,
                               vcl_vector< feature_sptr_type > const &  fixed,
                               vcl_vector< match_sptr_type >   const &  matches,
                               trans_sptr_type                 const &  trans,
                               unsigned int                             iteration );

  // Callback which calls member display_points to draw points on tableau.
  // It is done this way because 1) we need to know the caller, 2) callback must be static.
  static void display_points( void * caller,   
                              vcl_vector< feature_sptr_type > const &  moving,
                              vcl_vector< feature_sptr_type > const &  fixed,
                              vcl_vector< match_sptr_type >   const &  matches,
                              trans_sptr_type                 const &  transform,
                              unsigned int                             iteration )
  {
    cdcl_display_callback3d* mySelf = (cdcl_display_callback3d*)  caller;
    mySelf->display_points( moving, fixed, matches, transform, iteration );
    mySelf->finalize_display( iteration );
  }
  
  // Finalize display callback tableau.
  virtual void finalize_display( unsigned int  iteration );

protected:
  // Transfer error covariance.
  inline virtual vnl_matrix_fixed< double, dim, dim >  transfer_covar( vnl_matrix_fixed< double, dim, dof > const &  Jth,
                                                                       trans_sptr_type                      const &  trans ) {
    return Jth * trans->get_covariance() * vnl_transpose( Jth );
  }

  vgui_easy3D_tableau_sptr  easy3d_;

  // Resolution. Scaling factor by which all points will be multiplied.
  double res_;
};


// Callback which displays points is passed to the registration object.
// This way the registration object does not depend on vgui.
// Transfer error covariances are used.
template < unsigned int dim, unsigned int dof >
class cdcl_display_callback3d_transfer : public cdcl_display_callback3d< dim, dof >
{
public:
  cdcl_display_callback3d_transfer( vgui_easy3D_tableau_sptr &  easy3d )
    : cdcl_display_callback3d< dim, dof >( easy3d ) {};

  typedef cdcl_display_callback3d< dim, dof >             superclass_type;
  typedef typename superclass_type::feature_sptr_type     feature_sptr_type;
  typedef typename superclass_type::match_sptr_type       match_sptr_type;
  typedef typename superclass_type::trans_sptr_type       trans_sptr_type;
  typedef cdcl_estimation_transfer<dim, dof>              estimation_type;

  // Display 3d points.
  virtual void display_points( vcl_vector< feature_sptr_type > const &  moving,
                               vcl_vector< feature_sptr_type > const &  fixed,
                               vcl_vector< match_sptr_type >   const &  matches,
                               trans_sptr_type                 const &  trans,
                               unsigned int                             iteration )
  { superclass_type::display_points( moving, fixed, matches, trans, iteration ); };

  // Callback which calls member display_points to draw points on tableau.
  // It is done this way because 1) we need to know the caller, 2) callback must be static.
  static void display_points( void * caller,   
                              vcl_vector< feature_sptr_type > const &  moving,
                              vcl_vector< feature_sptr_type > const &  fixed,
                              vcl_vector< match_sptr_type >   const &  matches,
                              trans_sptr_type                 const &  transform,
                              unsigned int                             iteration )
  {
    cdcl_display_callback3d_transfer* mySelf = (cdcl_display_callback3d_transfer*)  caller;
    mySelf->display_points( moving, fixed, matches, transform, iteration );
    mySelf->finalize_display( iteration );
  }
  
protected:
  // Transfer error covariance.
  inline virtual vnl_matrix_fixed< double, dim, dim >  transfer_covar( vnl_matrix_fixed< double, dim, dof > const &  Jth,
                                                                       trans_sptr_type                      const &  trans ) {
    return trans->get_covarianceJ();
  }

};


#endif
