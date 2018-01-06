#include <vcl_iomanip.h>
#include <vcl_iostream.h>
#include <vcl_algorithm.h>
#include <vcl_utility.h>
#include <vcl_set.h>

#include <vul/vul_timer.h>

#include <vnl/vnl_inverse.h>
#include <vnl/vnl_transpose.h>
#include <vnl/vnl_trace.h>
#include <vnl/algo/vnl_svd.h>
#include <vnl/vnl_random.h>


#include "cdcl_estimation_symmetric_ICP_matching_all.h"

#include "cdcl_utils.h"
#include "cdcl_lbfgs.h"
#include "cdcl_trans_rigid3d.h"

//:
// \file
// \brief  Symmetric estimation.
//         The idea is to have forward and backward estimation object.
//         Constraints from inverse transform are included in the estimation.
//         This is the same as cdcl_estimation_symmetric_ICP but with 
//         cdcl_estimation_symmetric_ICP_matching_all as the estimation type.
// \author Michal Sofka
// \date   Feb 2008


// Construct estimation object given moving and fixed feature sets,
// and an initial transformation. All movind and fixed points are used.
template < unsigned int dim, unsigned int dof >
cdcl_estimation_symmetric_ICP_matching_all<dim, dof>::cdcl_estimation_symmetric_ICP_matching_all( vcl_vector< feature_sptr_type >  const &  moving, 
                                                                vcl_vector< feature_sptr_type >  const &  fixed,
                                                                trans_sptr_type                  const &  trans,
                                                                trans_sptr_type                  const &  trans_inverse )
  : normalize_matches_( true ),
    oscillation_count_( 0 ),
    error_difference_( 0.0 ),
    weighted_error_( 0.0 ),
    number_matches_( 1000 ),
    finest_level_( false )
{
  forward_->weight_by_strength( true );
  forward_->normalize_matches();
  backward_->weight_by_strength( true );
  backward_->normalize_matches();
};


// Construct estimation object given moving and fixed feature sets,
// an initial transformation and moving region corners.
template < unsigned int dim, unsigned int dof >
cdcl_estimation_symmetric_ICP_matching_all<dim, dof>::cdcl_estimation_symmetric_ICP_matching_all( vcl_vector< feature_sptr_type >  const &  moving, 
                                                                vcl_vector< feature_sptr_type >  const &  fixed,
                                                                trans_sptr_type                  const &  trans,
                                                                trans_sptr_type                  const &  trans_inverse,
                                                                vnl_vector_fixed< double, dim >  const &  moving_x0,
                                                                vnl_vector_fixed< double, dim >  const &  moving_x1,
                                                                vnl_vector_fixed< double, dim >  const &  fixed_x0,
                                                                vnl_vector_fixed< double, dim >  const &  fixed_x1 )
  : normalize_matches_( true ),
    oscillation_count_( 0 ),
    error_difference_( 0.0 ),
    weighted_error_( 0.0 ),
    number_matches_( 1000 ),
    finest_level_( false )
{
  typedef vcl_vector< feature_sptr_type >  feature_vector_type;
  // only consider points that are in the current region
  //   
  // extract fixed features from the current region only
  // create new copies so that normalization of each dataset individually is possible
  feature_vector_type  fixed_inside;
  for( typename feature_vector_type::size_type  i = 0; i < fixed.size(); ++i ) {
    vnl_vector_fixed< double, dim >  curr = fixed[i]->location_;
    bool inside = true;
    for( unsigned int d = 0; d < dim; ++d ) if( fixed_x0[d] > curr[d] || curr[d] > fixed_x1[d] ) inside = false;
    if( inside ) fixed_inside.push_back( new feature_type( *(fixed[i]) ) );
  }
  vcl_cout << "Features in fixed window: " << fixed_inside.size() << vcl_endl;

  // extract moving features from the current region only
  // create new copies so that normalization of each dataset individually is possible
  feature_vector_type  moving_inside;
  for( typename feature_vector_type::size_type  i = 0; i < moving.size(); ++i ) {
    vnl_vector_fixed< double, dim >  curr = moving[i]->location_;
    bool inside = true;
    for( unsigned int d = 0; d < dim; ++d ) if( moving_x0[d] > curr[d] || curr[d] > moving_x1[d] ) inside = false;
    if( inside ) moving_inside.push_back( new feature_type( *(moving[i]) ) );
  }
  vcl_cout << "Features in moving window: " << moving_inside.size() << vcl_endl;



  vnl_vector_fixed< double, dim >  fixed_size = fixed_x1 - fixed_x0;
  vnl_vector_fixed< double, dim >  fixed_x0_larger = fixed_x0 - 0.5 * fixed_size;
  vnl_vector_fixed< double, dim >  fixed_x1_larger = fixed_x1 + 0.5 * fixed_size;

  vnl_vector_fixed< double, dim >  moving_size = moving_x1 - fixed_x0;
  vnl_vector_fixed< double, dim >  moving_x0_larger = moving_x0 - 0.5 * moving_size;
  vnl_vector_fixed< double, dim >  moving_x1_larger = moving_x1 + 0.5 * moving_size;

  // only consider points that are in the current region
  //   
  // extract fixed features from the current region only
  // create new copies so that normalization of each dataset individually is possible
  feature_vector_type  fixed_inside_larger;
  for( typename feature_vector_type::size_type  i = 0; i < fixed.size(); ++i ) {
    vnl_vector_fixed< double, dim >  curr = fixed[i]->location_;
    bool inside = true;
    for( unsigned int d = 0; d < dim; ++d ) if( fixed_x0_larger[d] > curr[d] || curr[d] > fixed_x1_larger[d] ) inside = false;
    if( inside ) fixed_inside_larger.push_back( new feature_type( *(fixed[i]) ) );
  }

  // extract moving features from the current region only
  // create new copies so that normalization of each dataset individually is possible
  feature_vector_type  moving_inside_larger;
  for( typename feature_vector_type::size_type  i = 0; i < moving.size(); ++i ) {
    vnl_vector_fixed< double, dim >  curr = moving[i]->location_;
    bool inside = true;
    for( unsigned int d = 0; d < dim; ++d ) if( moving_x0_larger[d] > curr[d] || curr[d] > moving_x1_larger[d] ) inside = false;
    if( inside ) moving_inside_larger.push_back( new feature_type( *(moving[i]) ) );
  }



  forward_ = new estimation_type( moving_inside, fixed_inside_larger, trans );
  backward_ = new estimation_type( fixed_inside, moving_inside_larger, trans_inverse );

  forward_->weight_by_strength( true );
  forward_->normalize_matches();
  backward_->weight_by_strength( true );
  backward_->normalize_matches();
};


// Construct estimation object given moving and fixed feature sets,
// an initial transformation and moving region corners.
template < unsigned int dim, unsigned int dof >
cdcl_estimation_symmetric_ICP_matching_all<dim, dof>::cdcl_estimation_symmetric_ICP_matching_all( vcl_vector< feature_sptr_type >  const &  moving, 
                                                                vcl_vector< feature_sptr_type >  const &  fixed,
                                                                trans_sptr_type                  const &  trans,
                                                                trans_sptr_type                  const &  trans_inverse,
                                                                vcl_vector< feature_sptr_type >  const &  moving_inside, 
                                                                vcl_vector< feature_sptr_type >  const &  fixed_inside )
  : normalize_matches_( true ),
    oscillation_count_( 0 ),
    error_difference_( 0.0 ),
    weighted_error_( 0.0 ),
    number_matches_( 1000 ),
    finest_level_( false )
{
  forward_ = new estimation_type( moving_inside, fixed, trans );
  backward_ = new estimation_type( fixed_inside, moving, trans_inverse );

  forward_->weight_by_strength( true );
  forward_->normalize_matches();
  backward_->weight_by_strength( true );
  backward_->normalize_matches();
};


// Set the current moving and fixed regions. Only features withing the regions are used in matching.
template < unsigned int dim, unsigned int dof >
void cdcl_estimation_symmetric_ICP_matching_all<dim, dof>::set_regions( vcl_vector< feature_sptr_type >  const &  moving, 
                                                           vcl_vector< feature_sptr_type >  const &  fixed,
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
    bool inside = true;
    for( unsigned int d = 0; d < dim; ++d ) if( fixed_x0[d] > curr[d] || curr[d] > fixed_x1[d] ) inside = false;
    if( inside ) fixed_inside.push_back( new feature_type( *(fixed[i]) ) );
  }

  // extract moving features from the current region only
  // create new copies so that normalization of each dataset individually is possible
  feature_vector_type  moving_inside;
  for( typename feature_vector_type::size_type  i = 0; i < moving.size(); ++i ) {
    vnl_vector_fixed< double, dim >  curr = moving[i]->location_;
    bool inside = true;
    for( unsigned int d = 0; d < dim; ++d ) if( moving_x0[d] > curr[d] || curr[d] > moving_x1[d] ) inside = false;
    if( inside ) moving_inside.push_back( new feature_type( *(moving[i]) ) );
  }

  this->set_features( moving_inside, fixed_inside );
}


// Set the current moving and fixed regions. Only features withing the regions are used in matching.
template < unsigned int dim, unsigned int dof >
void cdcl_estimation_symmetric_ICP_matching_all<dim, dof>::set_features( vcl_vector< feature_sptr_type >  const &  moving_inside, 
                                                                         vcl_vector< feature_sptr_type >  const &  fixed_inside )
{
  typedef vcl_vector< feature_sptr_type >  feature_vector_type;

  vcl_vector< feature_vector_type >  multi_moving_inside;
  multi_moving_inside.push_back( moving_inside );

  vcl_vector< feature_vector_type >  multi_fixed_inside;
  multi_fixed_inside.push_back( fixed_inside );

  forward_->set_moving_data( multi_moving_inside );
  backward_->set_moving_data( multi_fixed_inside );

  // oscillation count and error_difference are no longer valid
  oscillation_count_ = 0;
  error_difference_ = 0.0;

  number_matches_ = 1000;
  finest_level_ = false;
}


// Set initial transform.
template < unsigned int dim, unsigned int dof >
void cdcl_estimation_symmetric_ICP_matching_all<dim, dof>::set_transform( trans_sptr_type const &  trans )
{
  forward_->set_transform( trans );
  backward_->set_transform( trans->inverse() );

  // oscillation count and error_difference are no longer valid
  oscillation_count_ = 0;
  error_difference_ = 0.0;

  number_matches_ = 1000;
  finest_level_ = false;
};


template < unsigned int dim, unsigned int dof >
void cdcl_estimation_symmetric_ICP_matching_all<dim, dof>::initialize()
{
  this->forward_->initialize();
  this->backward_->initialize();
}


template < unsigned int dim, unsigned int dof >
void cdcl_estimation_symmetric_ICP_matching_all<dim, dof>::find_closest_euclidean()
{
  forward_->find_closest_euclidean();
  backward_->find_closest_euclidean();
}


#define PROGRESS_OUTPUT( procedure, timer )  vcl_cout << "----> Computing " << procedure << " ... " << vcl_endl; timer.mark()
#define TIMER_OUTPUT( procedure, timer )  vcl_cout << "----> Computing " << procedure << " done in " << timer.real()/1000 << "." << timer.real()%1000 << " sec." << vcl_endl


// Run one EM iteration of estimation (parameters, weights, covariance, weights).
// Return true when converged.
template < unsigned int dim, unsigned int dof >
bool cdcl_estimation_symmetric_ICP_matching_all<dim, dof>::one_iteration( void *  caller, 
                                                         display_callback_type  display_points,
                                                         unsigned int  iteration )
{
  vul_timer timer;

  vcl_cout << "******************** Iteration " << iteration << "********************" << vcl_endl;

  double weighted_error_before = forward_->weighted_error();

  vcl_cout << "forward: " << forward_->moving_[forward_->res_].size() << "  " << forward_->fixed_[forward_->res_].size() << vcl_endl;
  vcl_cout << "backward: " << backward_->moving_[backward_->res_].size() << "  " << backward_->fixed_[backward_->res_].size() << vcl_endl;

  PROGRESS_OUTPUT( "weights", timer );
  this->find_closest_euclidean();
  TIMER_OUTPUT( "weights", timer );
  
  //if( caller ) display_points( caller, forward_->moving_[forward_->res_], forward_->fixed_[forward_->res_], forward_->matches_, forward_->trans_, iteration );
  // since the callback stores the points to be saved later, we would need another instance to display the other direction
  //if( caller ) display_points( caller, backward_->moving_[backward_->res_], backward_->fixed_[backward_->res_], backward_->matches_, backward_->trans_, 9000+iteration );

  PROGRESS_OUTPUT( "parameters", timer );
  this->estimate_LS( true );
  TIMER_OUTPUT( "parameters", timer );

  if( caller ) display_points( caller, forward_->moving_[forward_->res_], forward_->fixed_[forward_->res_], forward_->matches_, forward_->trans_, iteration );

  forward_->trans_->print( vcl_cout );
  backward_->trans_->print( vcl_cout );


  // convergence testing
  //
  double weighted_error_after = forward_->weighted_error();

  double curr_error_difference = weighted_error_after - weighted_error_before;
  
  // rate of change
  double diff = (weighted_error_after - weighted_error_before) / weighted_error_after ;

  vcl_cout << "Converged: after before rate_diff: " << weighted_error_after << " " << weighted_error_before << " " << diff << vcl_endl;

  bool converged = vcl_abs( diff ) < 1e-4;

  weighted_error_ = weighted_error_after;

  // if we are not at the finest level but started oscillating
  if( !finest_level_ && ( oscillation_count_ > 1 || ( vcl_abs( diff ) < 0.01 ) ) ) converged = true;  // not on the finest level, so allow switching (<- lower requirement on oscillation count)
  else
  if( !converged ) {
    // look for oscillation
    // There are two kinds of oscillation:           
    // 1. The error increases and decreases in turn
    // 2. Once it is good enough, the error increases slightly

    // first situation
    if( weighted_error_after > 1.5 ) {  // not good enough for termination
      if( curr_error_difference * error_difference_ < 0.0 ) { // change of sign
        ++oscillation_count_;
      }
      else {
        if( oscillation_count_ > 0 ) --oscillation_count_;
      }

      if( oscillated() ) converged = true;  // stuck in oscillation -> stop estimation
    }
    else { // second situation
      if( curr_error_difference > 0.0 ) { // error increases again
        ++oscillation_count_;
      }
//      if( finest_level_ ) {
        if( oscillation_count_ > 3 ) converged = true;
//      }
    }

  }
  vcl_cout << "Oscillation count: " << oscillation_count_ << vcl_endl;

  error_difference_ = curr_error_difference;

  if( converged ) {
    if( !finest_level_ ) { // increase the number of points used -> switch resolution
      converged = false;
      number_matches_ *= 2;
      bool forward_at_max = forward_->set_number_points( number_matches_ );
      bool backward_at_max = backward_->set_number_points( number_matches_ );
      finest_level_ = forward_at_max && backward_at_max;
      oscillation_count_ = 0;
      error_difference_ = 0.0;
    }
  }

  ++iteration;

  return converged && finest_level_;
}


// Add constraints from from_matches to to_matches.
template < unsigned int dim, unsigned int dof >
void cdcl_estimation_symmetric_ICP_matching_all<dim, dof>::add_constraints( vcl_vector< match_sptr_type >  const &  from_matches, vcl_vector< match_sptr_type > &  to_matches )
{
  // go through all moving points
  for( typename vcl_vector< feature_sptr_type >::size_type  i = 0; i < from_matches.size(); ++i ) {
    // moving point p
    feature_sptr_type  p = from_matches[i]->from_;
  
    vcl_vector< feature_sptr_type >  p_vec;
    feature_sptr_type  p_clone = new feature_type( *p );
    p_vec.push_back( p_clone );

    for( typename vcl_vector< feature_sptr_type >::size_type  j = 0; j < from_matches[i]->to_.size(); ++j ) {
      // fixed point q
      feature_sptr_type  q = from_matches[i]->to_[j];

      vcl_vector< double >  w_vec;
      w_vec.push_back( from_matches[i]->w_[j] );
   
      // create a new match with from / to roles reverse
      feature_sptr_type  q_clone = new feature_type( *q );
      match_sptr_type  new_match = new match_type( q_clone, p_vec, w_vec );
      to_matches.push_back( new_match );
    }
  }

}


// Estimate initial covariance matrix.
template < unsigned int dim, unsigned int dof >
void cdcl_estimation_symmetric_ICP_matching_all<dim, dof>::estimate_LS( bool estimate_parameters )
{
  trans_sptr_type  forward_global_transform;
  if( forward_->trans_->is_type( cdcl_trans_rigid3d::type_id() ) ) { // if transform is rigid with small angle approximation
    // from now on, transform is incremental
    // save current transform for recomposing
    forward_global_transform = forward_->trans_;
    forward_->trans_ = forward_->trans_->create_new(); // will create identity transformation

    forward_->trans_->set_covariance( forward_global_transform->get_covariance() );

    // for rigid transform, we are estimating incremental transform
    // therefore transform points first
    vnl_matrix_fixed< double, dim, dim > const &  Jp = forward_global_transform->jacobian_wrt_loc();
    // go through all moving points
//vcl_cout << "TRANSFORMING FOR RIGID" << vcl_endl;
    for( typename vcl_vector< match_sptr_type >::size_type  i = 0; i < forward_->matches_.size(); ++i ) {
      // moving point p
      feature_sptr_type  p = forward_->matches_[i]->from_;
      
      // reset the original from point
      // transform location and covariance
      const vnl_vector_fixed< double, dim >       location   = forward_global_transform->map_loc( p->location_ );
      vnl_matrix_fixed< double, dim, dim > covariance = Jp * p->covariance_ * Jp.transpose();
      forward_->matches_[i]->from_ = new feature_type( *p );
      forward_->matches_[i]->from_->location_ = location;
      forward_->matches_[i]->from_->covariance_ = covariance;
//vcl_cout << forward_->matches_[i]->from_->location_ << "  " << forward_->matches_[i]->to_[0]->location_ << vcl_endl;
    }
    forward_->trans_->set_parameterization( vnl_vector< double >( 6, 0.0 ) );
  }

  trans_sptr_type  backward_global_transform;
  if( backward_->trans_->is_type( cdcl_trans_rigid3d::type_id() ) ) { // if transform is rigid with small angle approximation
    // from now on, transform is incremental
    // save current transform for recomposing
    backward_global_transform = backward_->trans_;
    backward_->trans_ = backward_->trans_->create_new(); // will create identity transformation

    backward_->trans_->set_covariance( backward_global_transform->get_covariance() );

    // for rigid transform, we are estimating incremental transform
    // therefore transform points first
    vnl_matrix_fixed< double, dim, dim > const &  Jp = backward_global_transform->jacobian_wrt_loc();
    // go through all moving points
    for( typename vcl_vector< match_sptr_type >::size_type  i = 0; i < backward_->matches_.size(); ++i ) {
      // moving point p
      feature_sptr_type  p = backward_->matches_[i]->from_;
      
      // reset the original from point
      // transform location and covariance
      const vnl_vector_fixed< double, dim >       location   = backward_global_transform->map_loc( p->location_ );
      vnl_matrix_fixed< double, dim, dim > covariance = Jp * p->covariance_ * Jp.transpose();
      backward_->matches_[i]->from_ = new feature_type( *p );
      backward_->matches_[i]->from_->location_ = location;
      backward_->matches_[i]->from_->covariance_ = covariance;
    }
    backward_->trans_->set_parameterization( vnl_vector< double >( 6, 0.0 ) );
  }

  double forward_scale = forward_->estimate_scale_and_assign_weight();
  if( forward_scale < 0.005 ) forward_scale = 0.005; // in normalized coordinates

  double backward_scale = backward_->estimate_scale_and_assign_weight();
  if( backward_scale < 0.005 ) backward_scale = 0.005; // in normalized coordinates

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
    this->cdcl_normalize_matches_known( backward_->matches_, center_fixed, avg_rad_fixed, center_moving, avg_rad_moving );
    backward_->trans_->normalize( avg_rad_fixed, avg_rad_moving, center_fixed, center_moving );
  }






//  // output for debugging
//  //
//  // FORWARD
//  {
//  vnl_matrix_fixed< double, dof, dof >  AtA( 0.0 );
//  vnl_matrix_fixed< double, dof, 1 >    Atb( 0.0 );
//
//  vnl_matrix_fixed< double, dim, 1 >  q_loc;
//  double sum_w = 0.0;
//  for( typename vcl_vector< feature_sptr_type >::size_type  i = 0; i < forward_->matches_.size(); ++i ) {
//    // moving point p
//    feature_sptr_type  p = forward_->matches_[i]->from_;
//
//    vnl_vector_fixed< double, dim >  mapped_moving = forward_->trans_->map_loc( p->location_ );
//  
//    // Jacobian w.r.t. parameters
//    vnl_matrix_fixed< double, dim, dof > const &  Jth = forward_->trans_->jacobian_wrt_par( p->location_ );
//
//    vnl_matrix_fixed< double, dim, dim >  Cij = Jth * forward_->trans_->get_covariance() * Jth.transpose();
//
//    double Cij_trace = vnl_trace( Cij );
//
//    //vcl_cout << p->location_ << "  pc  " << trans_->center_moving_ << vcl_endl;
//
//    for( typename vcl_vector< feature_sptr_type >::size_type  j = 0; j < forward_->matches_[i]->to_.size(); ++j ) {
//      // fixed point q
//      feature_sptr_type  q = forward_->matches_[i]->to_[j];
//      if( forward_->trans_->is_type( cdcl_trans_rigid3d::type_id() ) ) { // if transform is rigid with small angle approximation
//        // rigid is: p' = p + R p + t
//        q_loc.set_column( 0, q->location_-p->location_ );
//      }
//      else {
//        q_loc.set_column( 0, q->location_ );
//      }
//
//      double weight = forward_->matches_[i]->w_[j];      
//
//      AtA += Jth.transpose() * q->error_projector_ * Jth * weight;
//      Atb += Jth.transpose() * q->error_projector_ * q_loc * weight;
// 
//      sum_w += weight;
//
//      //vcl_cout << Jth << vcl_endl << AtA << vcl_endl << vcl_endl;
//      
//      vnl_vector_fixed< double, dim >  const &  err_vec = mapped_moving - q->location_;
//
//      vcl_cout << p->location_ << "\t" << q->location_ << "\t" << weight << "\t" << Cij_trace << "\t" << err_vec.magnitude() << vcl_endl;
//
//    }
//  }
//  }




  // temporarily add constraints from backward to forward and from now on use only forward
  const typename vcl_vector< match_sptr_type >::size_type  end_of_forward = forward_->matches_.size();
  this->add_constraints( backward_->matches_, forward_->matches_ );

  if( forward_->weight_by_strength_ ) forward_->weight_by_strength();

forward_->weight_spatially();

  // Estimate initial matrix by setting up a least squares problem
  // A x = q
  // A^T A
  vnl_matrix_fixed< double, dof, dof >  AtA( 0.0 );
  vnl_matrix_fixed< double, dof, 1 >    Atb( 0.0 );
  
  // FORWARD
  {
  typedef typename vcl_vector< typename vcl_vector< feature_sptr_type >::size_type >  indices_vector_type;
  indices_vector_type  indices_vector( forward_->matches_.size(), 0 );
  for( typename indices_vector_type::size_type n = 0; n < forward_->matches_.size(); ++n ) {
    indices_vector[n] = n;
  }
  vcl_random_shuffle( indices_vector.begin(), indices_vector.end() );
  
  vnl_matrix_fixed< double, dim, 1 >  q_loc;
  //double sum_w = 0.0;
  // go through all moving points
  //vcl_random_shuffle( forward_->matches_.begin(), forward_->matches_.end() ); // Warning: currently can't use shuffle because we would later erase some forward matches, when deleting backward matches addition
//  if( number_matches_ >= forward_->matches_.size() ) {
//    number_matches_ = forward_->matches_.size();
//    finest_level_ = true;
//  }
  typename vcl_vector< feature_sptr_type >::size_type  matches_used = number_matches_;
  //matches_used = vcl_min( matches_used, forward_->matches_.size() );
  matches_used = forward_->matches_.size();
  //for( typename vcl_vector< feature_sptr_type >::size_type  i = end_of_forward; i < forward_->matches_.size(); ++i ) {
  for( typename vcl_vector< feature_sptr_type >::size_type  i = 0; i < matches_used; ++i ) {
  //for( typename vcl_vector< feature_sptr_type >::size_type  i = 0; i < end_of_forward; ++i ) {

    // moving point p
    //typename vcl_vector< feature_sptr_type >::size_type  ind = indices_vector[i];
    typename vcl_vector< feature_sptr_type >::size_type  ind = i;
    feature_sptr_type  p = forward_->matches_[ind]->from_;
  
    // Jacobian w.r.t. parameters
    vnl_matrix_fixed< double, dim, dof > const &  Jth = forward_->trans_->jacobian_wrt_par( p->location_ );
    //vcl_cout << p->location_ << "  pc  " << trans_->center_moving_ << vcl_endl;

    for( typename vcl_vector< feature_sptr_type >::size_type  j = 0; j < forward_->matches_[ind]->to_.size(); ++j ) {
      // fixed point q
      feature_sptr_type  q = forward_->matches_[ind]->to_[j];
      if( forward_->trans_->is_type( cdcl_trans_rigid3d::type_id() ) ) { // if transform is rigid with small angle approximation
        // rigid is: p' = p + R p + t
        q_loc.set_column( 0, q->location_-p->location_ );
      }
      else {
        q_loc.set_column( 0, q->location_ );
      }

      double weight = forward_->matches_[ind]->w_[j];

    //  AtA += Jth.transpose() * q->error_projector_ * Jth * weight;
    //  Atb += Jth.transpose() * q->error_projector_ * q_loc * weight;
 
    //  sum_w += weight;

    //  //vcl_cout << Jth << vcl_endl << AtA << vcl_endl << vcl_endl;

    //}


    //// Each has only one match

    //// fixed point q
    //feature_sptr_type  q = (*it)->to_[0];
    //
    //// removed the condition on rigid from above

    //q_loc.set_column( 0, q->location_ );


    AtA += Jth.transpose() * q->error_projector_ * Jth * weight;
    Atb += Jth.transpose() * q->error_projector_ * q_loc * weight;

    //sum_w += weight;

    //vcl_cout << Jth << vcl_endl << AtA << vcl_endl << vcl_endl;


  }


  //AtA *= ( 1.0 / sum_w );
  //Atb *= ( 1.0 / sum_w );
  }

  }

  // since forward / backward have their own scale the obj. fun is:  sum  1/2 forward + 1/2 backward
  // each has 3 different scales (for different feature shapes)
  //AtA *= ( 1.0 / 6.0 );
  //Atb *= ( 1.0 / 6.0 );


  vnl_svd< double > svd( AtA );
  vnl_matrix_fixed< double, dof, dof >  invAA = svd.inverse();

  // don't need to multiply by scale, since it is already included in the robust weights [yang:cvpr04]
  forward_->trans_->set_covariance( /*scale*scale **/ /** 2.0*/ /**/ /*scale_corners * scale_corners **/ invAA );

  // assign inverse covariance
  // note: would have to be estimated
  //backward_->trans_->set_covariance( );



  if( estimate_parameters ) {
    vnl_matrix_fixed< double, dof, 1 >  new_params = invAA * Atb;
    //vcl_cout << "PARAMS: " << new_params << vcl_endl;

    forward_->trans_->set_parameterization( new_params.get_column( 0 ) );

    //backward_->trans_ = forward_->trans_->inverse();
  }

  // erase the constraints that were added
  forward_->matches_.erase( forward_->matches_.begin() + end_of_forward, forward_->matches_.end() );

  if( normalize_matches_ ) {
    // convert transformation and matches back to unnormalized coordinates
    forward_->trans_->unnormalize( avg_rad_moving, avg_rad_fixed, center_moving, center_fixed );
    forward_->cdcl_unnormalize_matches( center_moving, avg_rad_moving, center_fixed, avg_rad_fixed );

    // unnormalize backward transformation and matches with flipped radii and centers
    backward_->trans_->unnormalize( avg_rad_fixed, avg_rad_moving, center_fixed, center_moving );
    backward_->cdcl_unnormalize_matches( center_fixed, avg_rad_fixed, center_moving, avg_rad_moving );
  }

  if( forward_->trans_->is_type( cdcl_trans_rigid3d::type_id() ) ) { // if transform is rigid with small angle approximation
    // new_transform is incremental rigid transform
    forward_global_transform->recompose_increment( forward_->trans_ );
    forward_global_transform->set_covariance( forward_->trans_->get_covariance() );
    forward_->trans_ = forward_global_transform;
    // from now on, transform is not incremental
  }

  //if( backward_->trans_->is_type( cdcl_trans_rigid3d::type_id() ) ) { // if transform is rigid with small angle approximation
  //  // new_transform is incremental rigid transform
  //  backward_global_transform->recompose_increment( backward_->trans_ );
  //  backward_global_transform->set_covariance( backward_->trans_->get_covariance() );
  //  backward_->trans_ = backward_global_transform;
  //  // from now on, transform is not incremental
  //}

  if( estimate_parameters ) {
    backward_->trans_ = forward_->trans_->inverse();
  }


  vcl_cout << "Initial forward covariance: " << vcl_endl << forward_->trans_->get_covariance() << vcl_endl << "Scale: " << forward_scale << vcl_endl;
  //vcl_cout << "Initial backward covariance: " << vcl_endl << backward_->trans_->get_covariance() << vcl_endl << "Scale: " << backward_scale << vcl_endl;

}


template < unsigned int dim, unsigned int dof >
void cdcl_estimation_symmetric_ICP_matching_all<dim, dof>::estimate_LS_backward( bool estimate_parameters )
{
  trans_sptr_type  backward_global_transform;
  if( backward_->trans_->is_type( cdcl_trans_rigid3d::type_id() ) ) { // if transform is rigid with small angle approximation
    // from now on, transform is incremental
    // save current transform for recomposing
    backward_global_transform = backward_->trans_;
    backward_->trans_ = backward_->trans_->create_new(); // will create identity transformation

    backward_->trans_->set_covariance( backward_global_transform->get_covariance() );

    // for rigid transform, we are estimating incremental transform
    // therefore transform points first
    vnl_matrix_fixed< double, dim, dim > const &  Jp = backward_global_transform->jacobian_wrt_loc();
    // go through all moving points
    for( typename vcl_vector< match_sptr_type >::size_type  i = 0; i < backward_->matches_.size(); ++i ) {
      // moving point p
      feature_sptr_type  p = backward_->matches_[i]->from_;
      
      // reset the original from point
      // transform location and covariance
      const vnl_vector_fixed< double, dim >       location   = backward_global_transform->map_loc( p->location_ );
      vnl_matrix_fixed< double, dim, dim > covariance = Jp * p->covariance_ * Jp.transpose();
      backward_->matches_[i]->from_ = new feature_type( *p );
      backward_->matches_[i]->from_->location_ = location;
      backward_->matches_[i]->from_->covariance_ = covariance;
    }
    backward_->trans_->set_parameterization( vnl_vector< double >( 6, 0.0 ) );
  }


  double forward_scale = forward_->estimate_scale_and_assign_weight();
  if( forward_scale < 0.005 ) forward_scale = 0.005; // in normalized coordinates

  double backward_scale = backward_->estimate_scale_and_assign_weight();
  if( backward_scale < 0.005 ) backward_scale = 0.005; // in normalized coordinates

  // temporarily add constraints from forward to backward and from now on use only backward
  const typename vcl_vector< match_sptr_type >::size_type  end_of_backward = backward_->matches_.size();
  this->add_constraints( forward_->matches_, backward_->matches_ );


  // Normalize before estimation.
  //
  // Normalization parameters.
  vnl_vector_fixed< double, dim >  center_moving( 0.0 ), center_fixed( 0.0 );
  double avg_rad_moving( 0.0 ), avg_rad_fixed( 0.0 );

  if( normalize_matches_ ) {
    // convert transformation and matches back to normalized coordinates
    backward_->cdcl_normalize_matches( center_moving, avg_rad_moving, center_fixed, avg_rad_fixed );
    backward_->trans_->normalize( avg_rad_moving, avg_rad_fixed, center_moving, center_fixed );

    // normalize forward transformation and matches with flipped radii and centers
    //this->cdcl_normalize_matches_known( forward_->matches_, center_fixed, avg_rad_fixed, center_moving, avg_rad_moving );
    //forward_->trans_->normalize( avg_rad_fixed, avg_rad_moving, center_fixed, center_moving );
  }

  if( backward_->weight_by_strength_ ) backward_->weight_by_strength();

  // Estimate initial matrix by setting up a least squares problem
  // A x = q
  // A^T A
  vnl_matrix_fixed< double, dof, dof >  AtA( 0.0 );
  vnl_matrix_fixed< double, dof, 1 >    Atb( 0.0 );
  
  // backward
  {
  vnl_matrix_fixed< double, dim, 1 >  q_loc;
  double sum_w = 0.0;
  // go through all moving points
  //for( typename vcl_vector< feature_sptr_type >::size_type  i = end_of_backward; i < backward_->matches_.size(); ++i ) {
  for( typename vcl_vector< feature_sptr_type >::size_type  i = 0; i < backward_->matches_.size(); ++i ) {
    // moving point p
    feature_sptr_type  p = backward_->matches_[i]->from_;
  
    // Jacobian w.r.t. parameters
    vnl_matrix_fixed< double, dim, dof > const &  Jth = backward_->trans_->jacobian_wrt_par( p->location_ );
    //vcl_cout << p->location_ << "  pc  " << trans_->center_moving_ << vcl_endl;

    for( typename vcl_vector< feature_sptr_type >::size_type  j = 0; j < backward_->matches_[i]->to_.size(); ++j ) {
      // fixed point q
      feature_sptr_type  q = backward_->matches_[i]->to_[j];
      if( backward_->trans_->is_type( cdcl_trans_rigid3d::type_id() ) ) { // if transform is rigid with small angle approximation
        // rigid is: p' = p + R p + t
        q_loc.set_column( 0, q->location_-p->location_ );
      }
      else {
        q_loc.set_column( 0, q->location_ );
      }

      double weight = backward_->matches_[i]->w_[j];

      AtA += Jth.transpose() * q->error_projector_ * Jth * weight;
      Atb += Jth.transpose() * q->error_projector_ * q_loc * weight;
 
      sum_w += weight;

      //vcl_cout << Jth << vcl_endl << AtA << vcl_endl << vcl_endl;

    }
  }
  //AtA *= ( 1.0 / sum_w );
  //Atb *= ( 1.0 / sum_w );
  }

  // since backward / forward have their own scale the obj. fun is:  sum  1/2 backward + 1/2 forward
  // each has 3 different scales (for different feature shapes)
  //AtA *= ( 1.0 / 6.0 );
  //Atb *= ( 1.0 / 6.0 );


  vnl_svd< double > svd( AtA );
  vnl_matrix_fixed< double, dof, dof >  invAA = svd.inverse();

  // don't need to multiply by scale, since it is already included in the robust weights [yang:cvpr04]
  backward_->trans_->set_covariance( /*scale*scale **/ /** 2.0*/ /**/ /*scale_corners * scale_corners **/ invAA );

  // assign inverse covariance
  // note: would have to be estimated
  //forward_->trans_->set_covariance( );

  if( estimate_parameters ) {
    vnl_matrix_fixed< double, dof, 1 >  new_params = invAA * Atb;
    //vcl_cout << "PARAMS: " << new_params << vcl_endl;

    backward_->trans_->set_parameterization( new_params.get_column( 0 ) );

    //forward_->trans_ = backward_->trans_->inverse();
  }

  // erase the constraints that were added
  backward_->matches_.erase( backward_->matches_.begin() + end_of_backward, backward_->matches_.end() );

  if( normalize_matches_ ) {
    // convert transformation and matches back to unnormalized coordinates
    backward_->trans_->unnormalize( avg_rad_moving, avg_rad_fixed, center_moving, center_fixed );
    backward_->cdcl_unnormalize_matches( center_moving, avg_rad_moving, center_fixed, avg_rad_fixed );

    // unnormalize forward transformation and matches with flipped radii and centers
    //forward_->trans_->unnormalize( avg_rad_fixed, avg_rad_moving, center_fixed, center_moving );
    //forward_->cdcl_unnormalize_matches( center_fixed, avg_rad_fixed, center_moving, avg_rad_moving );
  }

  if( estimate_parameters ) {
    forward_->trans_ = backward_->trans_->inverse();
  }

  vcl_cout << "Initial backward covariance: " << vcl_endl << backward_->trans_->get_covariance() << vcl_endl << "Scale: " << backward_scale << vcl_endl;
  //vcl_cout << "Initial forward covariance: " << vcl_endl << forward_->trans_->get_covariance() << vcl_endl << "Scale: " << forward_scale << vcl_endl;

  if( backward_->trans_->is_type( cdcl_trans_rigid3d::type_id() ) ) { // if transform is rigid with small angle approximation
    // new_transform is incremental rigid transform
    backward_global_transform->recompose_increment( backward_->trans_ );
    backward_global_transform->set_covariance( backward_->trans_->get_covariance() );
    backward_->trans_ = backward_global_transform;
    // from now on, transform is not incremental
  }
}


// Normalize moving and fixed points which the matches were formed with
// Use supplied known centers and radii
template < unsigned int dim, unsigned int dof >
void cdcl_estimation_symmetric_ICP_matching_all<dim, dof>::cdcl_normalize_matches_known( vcl_vector< match_sptr_type >        & matches,
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



#define CDCL_ESTIMATION_SYMMETRIC_ICP_MATCHING_ALL_INSTANTIATE( dim, dof ) \
template                                                  \
class cdcl_estimation_symmetric_ICP_matching_all< dim, dof >;
