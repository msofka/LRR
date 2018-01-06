#include <vcl_iomanip.h>
#include <vcl_iostream.h>
#include <vcl_set.h>

#include <vul/vul_timer.h>

#include "cdcl_estimation_symmetric.h"
#include "cdcl_obj_fun_par_sym.h"

#include "cdcl_utils.h"
#include "cdcl_lbfgs.h"

//:
// \file
// \brief  Symmetric estimation.
//         The idea is to have forward and backward estimation object.
//         We can then run matching and covariance estimation for forward and backward case
//         and provide parameter estimation function for combined symmetric case.
// \author Michal Sofka
// \date   May 2006


 //Construct estimation object given moving and fixed feature sets,
 //and an initial transformation. All movind and fixed points are used.
template < unsigned int dim, unsigned int dof, class estimation_type >
cdcl_estimation_symmetric<dim, dof, estimation_type>::cdcl_estimation_symmetric( vcl_vector< typename cdcl_estimation_symmetric<dim, dof, estimation_type>::feature_sptr_type >  const &  moving, 
                                                                vcl_vector< typename cdcl_estimation_symmetric<dim, dof, estimation_type>::feature_sptr_type >  const &  fixed,
                                                                typename cdcl_estimation_symmetric<dim, dof, estimation_type>::trans_sptr_type                  const &  trans,
                                                                typename cdcl_estimation_symmetric<dim, dof, estimation_type>::trans_sptr_type                  const &  trans_inverse )
  : normalize_matches_( false )
{
  forward_ = new estimation_type( moving, fixed, trans );
  backward_ = new estimation_type( fixed, moving, trans_inverse );
};


// Construct estimation object given moving and fixed feature sets,
// an initial transformation and moving region corners.
template < unsigned int dim, unsigned int dof, class estimation_type >
cdcl_estimation_symmetric<dim, dof, estimation_type>::cdcl_estimation_symmetric( vcl_vector< typename cdcl_estimation_symmetric<dim, dof, estimation_type>::feature_sptr_type >  const &  moving, 
                                                                vcl_vector< typename cdcl_estimation_symmetric<dim, dof, estimation_type>::feature_sptr_type >  const &  fixed,
                                                                typename cdcl_estimation_symmetric<dim, dof, estimation_type>::trans_sptr_type                  const &  trans,
                                                                typename cdcl_estimation_symmetric<dim, dof, estimation_type>::trans_sptr_type                  const &  trans_inverse,
                                                                vnl_vector_fixed< double, dim >  const &  moving_x0,
                                                                vnl_vector_fixed< double, dim >  const &  moving_x1,
                                                                vnl_vector_fixed< double, dim >  const &  fixed_x0,
                                                                vnl_vector_fixed< double, dim >  const &  fixed_x1 )
{
  typedef vcl_vector< feature_sptr_type >  feature_vector_type;
  // only consider points that are in the current region
  //   
  // extract fixed features from the current region only
  // create new copies so that normalization of each dataset individually is possible
  feature_vector_type  fixed_inside;
  for( typename feature_vector_type::size_type  i = 0; i < fixed.size(); ++i ) {
    vnl_vector_fixed< double, dim >  curr = fixed[i]->location_;
    if( fixed_x0[0] < curr[0] && curr[0] < fixed_x1[0]
     && fixed_x0[1] < curr[1] && curr[1] < fixed_x1[1] ) fixed_inside.push_back( new feature_type( *(fixed[i]) ) );
  }

  // extract moving features from the current region only
  // create new copies so that normalization of each dataset individually is possible
  feature_vector_type  moving_inside;
  for( typename feature_vector_type::size_type  i = 0; i < moving.size(); ++i ) {
    vnl_vector_fixed< double, dim >  curr = moving[i]->location_;
    if( moving_x0[0] < curr[0] && curr[0] < moving_x1[0]
     && moving_x0[1] < curr[1] && curr[1] < moving_x1[1] ) moving_inside.push_back( new feature_type( *(moving[i]) ) );
  }

  forward_ = new estimation_type( moving_inside, fixed, trans );
  backward_ = new estimation_type( fixed_inside, moving, trans_inverse );
};


#define PROGRESS_OUTPUT( procedure, timer )  vcl_cout << "----> Computing " << procedure << " ... " << vcl_endl; timer.mark()
#define TIMER_OUTPUT( procedure, timer )  vcl_cout << "----> Computing " << procedure << " done in " << timer.real()/1000 << "." << timer.real()%1000 << " sec." << vcl_endl


template < unsigned int dim, unsigned int dof, class estimation_type >
void cdcl_estimation_symmetric<dim, dof, estimation_type>::initialize()
{
  forward_->initialize();
  backward_->initialize();
}


// Run one EM iteration of estimation (parameters, weights, covariance, weights).
// Return true when converged.
template < unsigned int dim, unsigned int dof, class estimation_type >
bool cdcl_estimation_symmetric<dim, dof, estimation_type>::one_iteration( void *  caller, 
                                                         display_callback_type  display_points,
                                                         unsigned int  iteration )
{
  vul_timer timer;

  unsigned int startEM = 5;

  vcl_cout << "******************** Iteration " << iteration << "********************" << vcl_endl;

  bool resolution_switch_allowed = iteration > startEM;
  PROGRESS_OUTPUT( "weights", timer );
  forward_->compute_weights( resolution_switch_allowed );
  backward_->compute_weights( resolution_switch_allowed );
  TIMER_OUTPUT( "weights", timer );
  
  if( caller ) display_points( caller, forward_->moving_[forward_->res_], forward_->fixed_[forward_->res_], forward_->matches_, forward_->trans_, iteration );
  if( caller ) display_points( caller, backward_->moving_[backward_->res_], backward_->fixed_[backward_->res_], backward_->matches_, backward_->trans_, 9000+iteration );

  bool parameters_converged = false;
  if( iteration > startEM ) {
    PROGRESS_OUTPUT( "parameters", timer );
    parameters_converged = this->estimate_parameters();
    TIMER_OUTPUT( "parameters", timer );

    PROGRESS_OUTPUT( "weights", timer );
// WRONG FOR RIGID: WEIGHTS CANNOT BE COMPUTED HERE BECAUSE WE WOULD NEED TO RECOMPOSE TRANSFORM
    forward_->compute_weights( resolution_switch_allowed );
    backward_->compute_weights( resolution_switch_allowed );
// WRONG FOR RIGID: WEIGHTS CANNOT BE COMPUTED HERE BECAUSE WE WOULD NEED TO RECOMPOSE TRANSFORM
    TIMER_OUTPUT( "weights", timer );
  }

  PROGRESS_OUTPUT( "covariance", timer );
  bool covariance_converged = this->estimate_covariance();
  TIMER_OUTPUT( "covariance", timer );

  vcl_cout << "Forward:" << vcl_endl;
  forward_->get_transform()->print( vcl_cout );
  vcl_cout << "Backward:" << vcl_endl;
  backward_->get_transform()->print( vcl_cout );

  return parameters_converged && covariance_converged;
                              //&& forward_->get_transform()->get_covariance().frobenius_norm() < 1e-3
                              //&& backward_->get_transform()->get_covariance().frobenius_norm() < 1e-3;
}


// Run estimation.
template < unsigned int dim, unsigned int dof, class estimation_type >
bool cdcl_estimation_symmetric<dim, dof, estimation_type>::run( void *  caller, 
                                               display_callback_type  display_points )
{
  this->initialize();

  bool converged = false;
  unsigned int max_iterations = 70;
  unsigned int iteration = 0;

  while( !converged && iteration < max_iterations ) {

    converged = this->one_iteration( caller, display_points, iteration );

    ++iteration;
  }

  return converged;
}


// Estimate parameters.
template < unsigned int dim, unsigned int dof, class estimation_type >
bool cdcl_estimation_symmetric<dim, dof, estimation_type>::estimate_parameters()
{
  // Normalize before estimation.
  //
  // Normalization parameters.
  vnl_vector_fixed< double, dim >  center_moving( 0.0 ), center_fixed( 0.0 );
  double avg_rad_moving( 0.0 ), avg_rad_fixed( 0.0 );

  if( normalize_matches_ ) {
    // convert transformation and matches back to normalized coordinates
    forward_->cdcl_normalize_matches( center_moving, avg_rad_moving, center_fixed, avg_rad_fixed );
    forward_->trans_->normalize( avg_rad_moving, avg_rad_fixed, center_moving, center_fixed );

    // normalize backward transformation and matches with flipped radii and centers
    this->cdcl_normalize_matches_known( backward_->matches_, center_fixed, avg_rad_fixed , center_moving, avg_rad_moving );
    backward_->trans_->normalize( avg_rad_fixed, avg_rad_moving, center_fixed, center_moving );
  }

  // objective function
  cdcl_obj_fun_par_sym< dim, dof, typename estimation_type::obj_fun_par_type >  f_par_sym( forward_->matches_, forward_->trans_, backward_->matches_, backward_->trans_ );

  vnl_vector< double >  x = forward_->trans_->get_parameterization();  // can afford to copy the vector
  vnl_vector< double >  x_old = x;

  // estimator
  cdcl_lbfgs  minimizer( f_par_sym );
  minimizer.set_max_function_evals( 15 );


  //minimizer.set_f_tolerance( 0.0001 );
  //minimizer.set_x_tolerance( 1e-3 );
  minimizer.set_verbose( true );
  minimizer.set_check_derivatives( false );
  minimizer.set_trace( true );
  //minimizer.default_step_length = 0.1;
  bool converged = minimizer.minimize( x );
  
  if( ( !converged && minimizer.get_failure_code() < 0 ) || x.size() == 0 ) {
    vcl_cerr << "Error during minimization, failure code: " << minimizer.get_failure_code() << " . Keeping the old value." << vcl_endl;
    forward_->trans_->set_parameterization( x_old );
    //exit( 1 );
 }
  else {
    forward_->trans_->set_parameterization( x );
  }

  // minimizer start and end errors are computed as RMS
  bool parameters_converged = ( minimizer.get_start_error() - minimizer.get_end_error() ) < 1e-5;
  vcl_cout << "converged: " << converged << " parameters_converged: " << parameters_converged << " start error: " << minimizer.get_start_error() << " end error: " << minimizer.get_end_error() << " diff: " << minimizer.get_start_error() - minimizer.get_end_error() << vcl_endl;

  // update the backward transformation, make sure to keep the same covariance
  vnl_matrix_fixed< double, dof, dof >  backward_covar = backward_->trans_->get_covariance();
  vnl_matrix_fixed< double, dim, dim >  backward_covarJ = backward_->trans_->get_covarianceJ();
  backward_->set_transform( forward_->trans_->inverse() );
  backward_->trans_->set_covariance( backward_covar );
  backward_->trans_->set_covarianceJ( backward_covarJ );

  if( normalize_matches_ ) {
    // convert transformation and matches back to unnormalized coordinates
    forward_->trans_->unnormalize( avg_rad_moving, avg_rad_fixed, center_moving, center_fixed );
    forward_->cdcl_unnormalize_matches( center_moving, avg_rad_moving, center_fixed, avg_rad_fixed );

    // unnormalize backward transformation and matches with flipped radii and centers
    backward_->trans_->unnormalize( avg_rad_fixed, avg_rad_moving, center_fixed, center_moving );
    backward_->cdcl_unnormalize_matches( center_fixed, avg_rad_fixed, center_moving, avg_rad_moving );
  }

  return parameters_converged;
}


// Estimate covariance.
template < unsigned int dim, unsigned int dof, class estimation_type >
bool cdcl_estimation_symmetric<dim, dof, estimation_type>::estimate_covariance()
{
  bool forward_covariance_converged = forward_->estimate_covariance();
  bool backward_covariance_converged = backward_->estimate_covariance();
  return forward_covariance_converged && backward_covariance_converged;
}


// Normalize moving and fixed points which the matches were formed with
// Use supplied known centers and radii
template < unsigned int dim, unsigned int dof, class estimation_type >
void cdcl_estimation_symmetric<dim, dof, estimation_type>::cdcl_normalize_matches_known( vcl_vector< match_sptr_type >        & matches,
                                   vnl_vector_fixed< typename feature_type::coord_type, dim >                   const & center_moving, 
                                   typename feature_type::coord_type                                            const & avg_radius_moving, 
                                   vnl_vector_fixed< typename feature_type::coord_type, dim >                   const & center_fixed, 
                                   typename feature_type::coord_type                                            const & avg_radius_fixed )
{
  // set of fixed points that got matched to moving points
  // utilize unique container concept, i.e. no two points are identical
  // we need to make sure that points don't get normalized twice
  vcl_set< feature_sptr_type >  fixed_matched;

  // go through all moving points
  for( typename vcl_vector< match_sptr_type >::size_type  i = 0; i < matches.size(); ++i ) {
    // moving point p
    feature_sptr_type  p = matches[i]->from_;

    // convert locations to normalized coordinate system
    p->location_ -= center_moving;
    p->location_ /= avg_radius_moving;

    // convert covariances to normalized coordinate system
    p->covariance_ /= (avg_radius_moving*avg_radius_moving);

    // go through all matches (fixed points)
    for( typename vcl_vector< match_sptr_type >::size_type  j = 0; j < matches[i]->to_.size(); ++j ) {
      // fixed point q
      feature_sptr_type  q = matches[i]->to_[j];
      // if it was inserted (i.e. not already there)
      if( fixed_matched.insert( q ).second ) {
        // convert locations to normalized coordinate system
        q->location_ -= center_fixed;
        q->location_ /= avg_radius_fixed;

        // convert covariances to normalized coordinate system
        q->covariance_ /= (avg_radius_fixed*avg_radius_fixed);
      }
    }
  }

}


#define CDCL_ESTIMATION_SYMMETRIC_INSTANTIATE( dim, dof, estimation_type ) \
template                                                  \
class cdcl_estimation_symmetric< dim, dof, estimation_type >;
