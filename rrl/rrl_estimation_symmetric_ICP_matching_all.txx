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


#include "rrl_estimation_symmetric_ICP_matching_all.h"

#include <cdcl/cdcl_utils.h>
#include <cdcl/cdcl_lbfgs.h>
#include <cdcl/cdcl_trans_rigid3d.h>

//:
// \file
// \brief  Symmetric estimation.
//         The idea is to have forward and backward estimation object.
//         Constraints from inverse transform are included in the estimation.
//         This is the same as cdcl_estimation_symmetric_ICP but with 
//         rrl_estimation_symmetric_ICP_matching_all as the estimation type.
// \author Michal Sofka
// \date   Feb 2008

//#define DEBUG( statement ) statement
#define DEBUG( statement ) 0



// computes
//    AtA += Jth.transpose() * p->error_projector_ * Jth * weight;
//    Atb += Jth.transpose() * p->error_projector_ * q_loc * weight;
// the first multiplication and transpose is done together
// jacobian has many zeros, so some multiplications are not necessary
// the second multiplication is done in a for loop
// only upper triangle of AtA is filled
#define COMPUTE_ATA( error_proj, point_location, AtA, Atb )               \
    Jth_point_proj( 0, 0 ) = Jth( 0, 0 ) * error_proj( 0, 0 );            \
    Jth_point_proj( 0, 1 ) = Jth( 0, 0 ) * error_proj( 0, 1 );            \
    Jth_point_proj( 0, 2 ) = Jth( 0, 0 ) * error_proj( 0, 2 );            \
                                                                          \
    Jth_point_proj( 1, 0 ) = Jth( 0, 1 ) * error_proj( 0, 0 );            \
    Jth_point_proj( 1, 1 ) = Jth( 0, 1 ) * error_proj( 0, 1 );            \
    Jth_point_proj( 1, 2 ) = Jth( 0, 1 ) * error_proj( 0, 2 );            \
                                                                          \
    Jth_point_proj( 2, 0 ) = Jth( 0, 2 ) * error_proj( 0, 0 );            \
    Jth_point_proj( 2, 1 ) = Jth( 0, 2 ) * error_proj( 0, 1 );            \
    Jth_point_proj( 2, 2 ) = Jth( 0, 2 ) * error_proj( 0, 2 );            \
                                                                          \
    Jth_point_proj( 3, 0 ) = Jth( 1, 3 ) * error_proj( 1, 0 );            \
    Jth_point_proj( 3, 1 ) = Jth( 1, 3 ) * error_proj( 1, 1 );            \
    Jth_point_proj( 3, 2 ) = Jth( 1, 3 ) * error_proj( 1, 2 );            \
                                                                          \
    Jth_point_proj( 4, 0 ) = Jth( 1, 4 ) * error_proj( 1, 0 );            \
    Jth_point_proj( 4, 1 ) = Jth( 1, 4 ) * error_proj( 1, 1 );            \
    Jth_point_proj( 4, 2 ) = Jth( 1, 4 ) * error_proj( 1, 2 );            \
                                                                          \
    Jth_point_proj( 5, 0 ) = Jth( 1, 5 ) * error_proj( 1, 0 );            \
    Jth_point_proj( 5, 1 ) = Jth( 1, 5 ) * error_proj( 1, 1 );            \
    Jth_point_proj( 5, 2 ) = Jth( 1, 5 ) * error_proj( 1, 2 );            \
                                                                          \
    Jth_point_proj( 6, 0 ) = Jth( 2, 6 ) * error_proj( 2, 0 );            \
    Jth_point_proj( 6, 1 ) = Jth( 2, 6 ) * error_proj( 2, 1 );            \
    Jth_point_proj( 6, 2 ) = Jth( 2, 6 ) * error_proj( 2, 2 );            \
                                                                          \
    Jth_point_proj( 7, 0 ) = Jth( 2, 7 ) * error_proj( 2, 0 );            \
    Jth_point_proj( 7, 1 ) = Jth( 2, 7 ) * error_proj( 2, 1 );            \
    Jth_point_proj( 7, 2 ) = Jth( 2, 7 ) * error_proj( 2, 2 );            \
                                                                          \
    Jth_point_proj( 8, 0 ) = Jth( 2, 8 ) * error_proj( 2, 0 );            \
    Jth_point_proj( 8, 1 ) = Jth( 2, 8 ) * error_proj( 2, 1 );            \
    Jth_point_proj( 8, 2 ) = Jth( 2, 8 ) * error_proj( 2, 2 );            \
                                                                          \
    Jth_point_proj( 9, 0 ) = error_proj( 0, 0 );                          \
    Jth_point_proj( 9, 1 ) = error_proj( 0, 1 );                          \
    Jth_point_proj( 9, 2 ) = error_proj( 0, 2 );                          \
                                                                          \
    Jth_point_proj( 10, 0 ) = error_proj( 1, 0 );                         \
    Jth_point_proj( 10, 1 ) = error_proj( 1, 1 );                         \
    Jth_point_proj( 10, 2 ) = error_proj( 1, 2 );                         \
                                                                          \
    Jth_point_proj( 11, 0 ) = error_proj( 2, 0 );                         \
    Jth_point_proj( 11, 1 ) = error_proj( 2, 1 );                         \
    Jth_point_proj( 11, 2 ) = error_proj( 2, 2 );                         \
                                                                          \
    for (unsigned int i = 0; i < dof; ++i)                                \
      for (unsigned int j = i; j < dof; ++j) {                            \
        double accum = Jth_point_proj[i][0] * Jth[0][j];                  \
        for (unsigned int k = 1; k < dim; ++k)                            \
          accum += Jth_point_proj[i][k] * Jth[k][j];                      \
        AtA[i][j] += accum;                                               \
      }                                                                   \
    for (unsigned int i = 0; i < dof; ++i) {                              \
      double accum = Jth_point_proj[i][0] * q->location_[0];              \
      for (unsigned int k = 1; k < dim; ++k)                              \
        accum += Jth_point_proj[i][k] * point_location[k];                \
      Atb[i][0] += accum;                                                 \
    }



// Construct estimation object given moving and fixed feature sets,
// and an initial transformation. All movind and fixed points are used.
template < unsigned int dim, unsigned int dof >
rrl_estimation_symmetric_ICP_matching_all<dim, dof>::rrl_estimation_symmetric_ICP_matching_all( vcl_vector< feature_sptr_type >  const &  moving, 
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

  m_Threader = ThreaderType::New();
};


// Construct estimation object given moving and fixed feature sets,
// an initial transformation and moving region corners.
template < unsigned int dim, unsigned int dof >
rrl_estimation_symmetric_ICP_matching_all<dim, dof>::rrl_estimation_symmetric_ICP_matching_all( vcl_vector< feature_sptr_type >  const &  moving, 
                                                                vcl_vector< feature_sptr_type >  const &  fixed,
                                                                trans_sptr_type                  const &  trans,
                                                                trans_sptr_type                  const &  trans_inverse,
                                                                location_type  const &  moving_x0,
                                                                location_type  const &  moving_x1,
                                                                location_type  const &  fixed_x0,
                                                                location_type  const &  fixed_x1,
																																bool constructed )
  : normalize_matches_( true ),
    oscillation_count_( 0 ),
    error_difference_( 0.0 ),
    weighted_error_( 0.0 ),
    number_matches_( 1000 ),
    finest_level_( false )
{
	constructed = true;

  typedef vcl_vector< feature_sptr_type >  feature_vector_type;
  // only consider points that are in the current region
  //   
  // extract fixed features from the current region only
  // create new copies so that normalization of each dataset individually is possible
  feature_vector_type  fixed_inside;
  for( typename feature_vector_type::size_type  i = 0; i < fixed.size(); ++i ) {
    location_type  curr = fixed[i]->location_;
    bool inside = true;
    for( unsigned int d = 0; d < dim; ++d ) if( fixed_x0[d] > curr[d] || curr[d] > fixed_x1[d] ) inside = false;
    if( inside ) fixed_inside.push_back( new feature_type( *(fixed[i]) ) );
  }
  vcl_cout << "Features in fixed window: " << fixed_inside.size() << vcl_endl;

  // extract moving features from the current region only
  // create new copies so that normalization of each dataset individually is possible
  feature_vector_type  moving_inside;
  for( typename feature_vector_type::size_type  i = 0; i < moving.size(); ++i ) {
    location_type  curr = moving[i]->location_;
    bool inside = true;
    for( unsigned int d = 0; d < dim; ++d ) if( moving_x0[d] > curr[d] || curr[d] > moving_x1[d] ) inside = false;
    if( inside ) moving_inside.push_back( new feature_type( *(moving[i]) ) );
  }
  vcl_cout << "Features in moving window: " << moving_inside.size() << vcl_endl;

	if( moving_inside.size() < 1000 || fixed_inside.size() < 1000 ) {
		vcl_cout << "Warning: Not enough features in the regions, fixed: " << fixed_inside.size() << " moving: " << moving_inside.size() << vcl_endl;
		constructed = false;
	}

  forward_ = new estimation_type( moving_inside, fixed, trans );
  backward_ = new estimation_type( fixed_inside, moving, trans_inverse );

  forward_->weight_by_strength( true );
  forward_->normalize_matches();
  backward_->weight_by_strength( true );
  backward_->normalize_matches();

  m_Threader = ThreaderType::New();
};


// Construct estimation object given moving and fixed feature sets,
// an initial transformation and moving region corners.
template < unsigned int dim, unsigned int dof >
rrl_estimation_symmetric_ICP_matching_all<dim, dof>::rrl_estimation_symmetric_ICP_matching_all( vcl_vector< feature_sptr_type >  const &  moving, 
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

  m_Threader = ThreaderType::New();
};


// Set the current moving and fixed regions. Only features withing the regions are used in matching.
template < unsigned int dim, unsigned int dof >
bool rrl_estimation_symmetric_ICP_matching_all<dim, dof>::set_regions( vcl_vector< feature_sptr_type >  const &  moving, 
                                                           vcl_vector< feature_sptr_type >  const &  fixed,
                                                           location_type  const &  moving_x0,
                                                           location_type  const &  moving_x1,
                                                           location_type  const &  fixed_x0,
                                                           location_type  const &  fixed_x1 )
{
  typedef vcl_vector< feature_sptr_type >  feature_vector_type;

  // only consider points that are in the current region
  //   
  // extract fixed features from the current region only
  // create new copies so that normalization of each dataset individually is possible
  feature_vector_type  fixed_inside;
  for( typename feature_vector_type::size_type  i = 0; i < fixed.size(); ++i ) {
    location_type  curr = fixed[i]->location_;
    bool inside = true;
    for( unsigned int d = 0; d < dim; ++d ) if( fixed_x0[d] > curr[d] || curr[d] > fixed_x1[d] ) inside = false;
    if( inside ) fixed_inside.push_back( new feature_type( *(fixed[i]) ) );
  }

  // extract moving features from the current region only
  // create new copies so that normalization of each dataset individually is possible
  feature_vector_type  moving_inside;
  for( typename feature_vector_type::size_type  i = 0; i < moving.size(); ++i ) {
    location_type  curr = moving[i]->location_;
    bool inside = true;
    for( unsigned int d = 0; d < dim; ++d ) if( moving_x0[d] > curr[d] || curr[d] > moving_x1[d] ) inside = false;
    if( inside ) moving_inside.push_back( new feature_type( *(moving[i]) ) );
  }

  this->set_features( moving_inside, fixed_inside );
  
  if( moving_inside.size() < 1000 || fixed_inside.size() < 1000 ) {
    vcl_cout << "Warning: Not enough features in the regions, fixed: " << fixed_inside.size() << " moving: " << moving_inside.size() << vcl_endl;
    return false;
  }

  return true;
}


// Set the current moving and fixed regions. Only features withing the regions are used in matching.
template < unsigned int dim, unsigned int dof >
void rrl_estimation_symmetric_ICP_matching_all<dim, dof>::set_features( vcl_vector< feature_sptr_type >  const &  moving_inside, 
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


// Set voronoi maps to forward and backward estimators to speed up matching.
template < unsigned int dim, unsigned int dof >
void rrl_estimation_symmetric_ICP_matching_all<dim, dof>::set_voronoi_maps( VoronoiMapImageType::Pointer const &  moving_voronoi_map,
                                                                            VoronoiMapImageType::Pointer const &  fixed_voronoi_map )
{
  forward_->set_fixed_voronoi_map( fixed_voronoi_map );
  backward_->set_fixed_voronoi_map( moving_voronoi_map );
}


// Set initial transform.
template < unsigned int dim, unsigned int dof >
void rrl_estimation_symmetric_ICP_matching_all<dim, dof>::set_transform( trans_sptr_type const &  trans )
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
void rrl_estimation_symmetric_ICP_matching_all<dim, dof>::initialize()
{
  this->forward_->initialize();
  this->backward_->initialize();
}


template < unsigned int dim, unsigned int dof >
void rrl_estimation_symmetric_ICP_matching_all<dim, dof>::find_closest_euclidean()
{
  forward_->find_closest_euclidean();
  backward_->find_closest_euclidean();
}


#define PROGRESS_OUTPUT( procedure, timer )  vcl_cout << "----> Computing " << procedure << " ... " << vcl_endl; timer.mark()
#define TIMER_OUTPUT( procedure, timer )  vcl_cout << "----> Computing " << procedure << " done in " << timer.real()/1000 << "." << timer.real()%1000 << " sec." << vcl_endl


// Run one EM iteration of estimation (parameters, weights, covariance, weights).
// Return true when converged.
template < unsigned int dim, unsigned int dof >
bool rrl_estimation_symmetric_ICP_matching_all<dim, dof>::one_iteration( void *  caller, 
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

  DEBUG( forward_->trans_->print( vcl_cout ) );
  DEBUG( backward_->trans_->print( vcl_cout ) );


  // convergence testing
  //
  double weighted_error_after = forward_->weighted_error();

  double curr_error_difference = weighted_error_after - weighted_error_before;
  
  // rate of change
  double diff = (weighted_error_after - weighted_error_before) / weighted_error_after ;

  DEBUG( vcl_cout << "Converged: after before rate_diff: " << weighted_error_after << " " << weighted_error_before << " " << diff << vcl_endl );

  bool converged = vcl_abs( diff ) < 1e-4;

  weighted_error_ = weighted_error_after;
  
  if( iteration > 0 ) {
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
  }
  DEBUG( vcl_cout << "Oscillation count: " << oscillation_count_ << vcl_endl );

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

template < unsigned int dim, unsigned int dof >
void rrl_estimation_symmetric_ICP_matching_all<dim, dof>::ThreadedGenerateDataForward( unsigned int ThreadId ) 
{


  for( typename vcl_vector< feature_sptr_type >::size_type  i = m_Indices[ThreadId]; i < m_Indices[ThreadId+1]; ++i ) {
  //for( typename vcl_vector< feature_sptr_type >::size_type  i = 0; i < end_of_forward; ++i ) {
    // moving point p
    //typename vcl_vector< feature_sptr_type >::size_type  ind = indices_vector[i];
    typename vcl_vector< feature_sptr_type >::size_type  ind = i;
    // vbl_smart_ptr is not thread safe -> don't use them
    feature_type*  p = forward_->matches_[ind]->from_.as_pointer();
  
    // Jacobian w.r.t. parameters
    vnl_matrix_fixed< coord_type, dim, dof >  Jth;
    forward_->trans_->jacobian_wrt_par_thread( p->location_, Jth );
    //vcl_cout << p->location_ << "  pc  " << trans_->center_moving_ << vcl_endl;

    // fixed point q
    feature_type*  q = forward_->matches_[ind]->to_[0].as_pointer();
    //q_loc.set_column( 0, q->location_ );

    const double &  weight = forward_->matches_[ind]->w_[0];
if( weight < 1e-4 ) continue;
    //AtA += Jth.transpose() * q->error_projector_ * Jth * weight;
    //Atb += Jth.transpose() * q->error_projector_ * q_loc * weight;

    vnl_matrix_fixed< coord_type, dof, dim >  Jth_q_proj;
    vnl_matrix_fixed< coord_type, dim, dim >  error_proj = q->error_projector_ * weight;

// THERE IS A BUG IN THE MACRO, it is cubersome anyways
//    COMPUTE_ATA( error_proj, q->location_, AtAs[ThreadId], Atbs[ThreadId] );

    // the first multiplication and transpose is done together
    // jacobian has many zeros, so some multiplications are saved
    Jth_q_proj( 0, 0 ) = Jth( 0, 0 ) * error_proj( 0, 0 );
    Jth_q_proj( 0, 1 ) = Jth( 0, 0 ) * error_proj( 0, 1 );
    Jth_q_proj( 0, 2 ) = Jth( 0, 0 ) * error_proj( 0, 2 );

    Jth_q_proj( 1, 0 ) = Jth( 0, 1 ) * error_proj( 0, 0 );
    Jth_q_proj( 1, 1 ) = Jth( 0, 1 ) * error_proj( 0, 1 );
    Jth_q_proj( 1, 2 ) = Jth( 0, 1 ) * error_proj( 0, 2 );

    Jth_q_proj( 2, 0 ) = Jth( 0, 2 ) * error_proj( 0, 0 );
    Jth_q_proj( 2, 1 ) = Jth( 0, 2 ) * error_proj( 0, 1 );
    Jth_q_proj( 2, 2 ) = Jth( 0, 2 ) * error_proj( 0, 2 );

    Jth_q_proj( 3, 0 ) = Jth( 1, 3 ) * error_proj( 1, 0 );
    Jth_q_proj( 3, 1 ) = Jth( 1, 3 ) * error_proj( 1, 1 );
    Jth_q_proj( 3, 2 ) = Jth( 1, 3 ) * error_proj( 1, 2 );
                        
    Jth_q_proj( 4, 0 ) = Jth( 1, 4 ) * error_proj( 1, 0 );
    Jth_q_proj( 4, 1 ) = Jth( 1, 4 ) * error_proj( 1, 1 );
    Jth_q_proj( 4, 2 ) = Jth( 1, 4 ) * error_proj( 1, 2 );
                        
    Jth_q_proj( 5, 0 ) = Jth( 1, 5 ) * error_proj( 1, 0 );
    Jth_q_proj( 5, 1 ) = Jth( 1, 5 ) * error_proj( 1, 1 );
    Jth_q_proj( 5, 2 ) = Jth( 1, 5 ) * error_proj( 1, 2 );

    Jth_q_proj( 6, 0 ) = Jth( 2, 6 ) * error_proj( 2, 0 );
    Jth_q_proj( 6, 1 ) = Jth( 2, 6 ) * error_proj( 2, 1 );
    Jth_q_proj( 6, 2 ) = Jth( 2, 6 ) * error_proj( 2, 2 );
                        
    Jth_q_proj( 7, 0 ) = Jth( 2, 7 ) * error_proj( 2, 0 );
    Jth_q_proj( 7, 1 ) = Jth( 2, 7 ) * error_proj( 2, 1 );
    Jth_q_proj( 7, 2 ) = Jth( 2, 7 ) * error_proj( 2, 2 );
                        
    Jth_q_proj( 8, 0 ) = Jth( 2, 8 ) * error_proj( 2, 0 );
    Jth_q_proj( 8, 1 ) = Jth( 2, 8 ) * error_proj( 2, 1 );
    Jth_q_proj( 8, 2 ) = Jth( 2, 8 ) * error_proj( 2, 2 );

    Jth_q_proj( 9, 0 ) = error_proj( 0, 0 );
    Jth_q_proj( 9, 1 ) = error_proj( 0, 1 );
    Jth_q_proj( 9, 2 ) = error_proj( 0, 2 );

    Jth_q_proj( 10, 0 ) = error_proj( 1, 0 );
    Jth_q_proj( 10, 1 ) = error_proj( 1, 1 );
    Jth_q_proj( 10, 2 ) = error_proj( 1, 2 );

    Jth_q_proj( 11, 0 ) = error_proj( 2, 0 );
    Jth_q_proj( 11, 1 ) = error_proj( 2, 1 );
    Jth_q_proj( 11, 2 ) = error_proj( 2, 2 );

    for (unsigned int i = 0; i < dof; ++i)
      for (unsigned int j = i; j < dof; ++j) {
        double accum = Jth_q_proj[i][0] * Jth[0][j];
        for (unsigned int k = 1; k < dim; ++k)
          accum += Jth_q_proj[i][k] * Jth[k][j];
        AtAs[ThreadId][i][j] += accum;
        //AtA[i][j] = AtA[j][i] += accum;
      }
    for (unsigned int i = 0; i < dof; ++i) {
      double accum = Jth_q_proj[i][0] * q->location_[0];
      for (unsigned int k = 1; k < dim; ++k)
        accum += Jth_q_proj[i][k] * q->location_[k];
      Atbs[ThreadId][i][0] += accum;
    }

  }



}


template < unsigned int dim, unsigned int dof >
void rrl_estimation_symmetric_ICP_matching_all<dim, dof>::ThreadedGenerateDataBackward( unsigned int ThreadId ) 
{



  for( typename vcl_vector< feature_sptr_type >::size_type  i = m_Indices[ThreadId]; i < m_Indices[ThreadId+1]; ++i ) {
  //for( typename vcl_vector< feature_sptr_type >::size_type  i = 0; i < end_of_forward; ++i ) {
    // moving point p
    //typename vcl_vector< feature_sptr_type >::size_type  ind = indices_vector[i];
    typename vcl_vector< feature_sptr_type >::size_type  ind = i;
    // vbl_smart_ptr is not thread safe -> don't use them
    feature_type*  p = backward_->matches_[ind]->from_.as_pointer();
  
    // fixed point q
    feature_type*  q = backward_->matches_[ind]->to_[0].as_pointer();

    // Jacobian w.r.t. parameters
    vnl_matrix_fixed< coord_type, dim, dof > Jth;
    forward_->trans_->jacobian_wrt_par_thread( q->location_, Jth );
    //vcl_cout << p->location_ << "  pc  " << trans_->center_moving_ << vcl_endl;

    //q_loc.set_column( 0, p->location_ );

    const double &  weight = backward_->matches_[ind]->w_[0];

if( weight < 1e-4 ) continue;

    //AtA += Jth.transpose() * p->error_projector_ * Jth * weight;
    //Atb += Jth.transpose() * p->error_projector_ * q_loc * weight;

    vnl_matrix_fixed< coord_type, dof, dim >  Jth_p_proj;
    vnl_matrix_fixed< coord_type, dim, dim >  error_proj = p->error_projector_ * weight;

// THERE IS A BUG IN THE MACRO, it is cubersome anyways
//    COMPUTE_ATA( error_proj, p->location_, AtAs[ThreadId], Atbs[ThreadId] );

//#if 0
    // the first multiplication and transpose is done together
    // jacobian has many zeros, so some multiplications are saved
    Jth_p_proj( 0, 0 ) = Jth( 0, 0 ) * error_proj( 0, 0 );
    Jth_p_proj( 0, 1 ) = Jth( 0, 0 ) * error_proj( 0, 1 );
    Jth_p_proj( 0, 2 ) = Jth( 0, 0 ) * error_proj( 0, 2 );

    Jth_p_proj( 1, 0 ) = Jth( 0, 1 ) * error_proj( 0, 0 );
    Jth_p_proj( 1, 1 ) = Jth( 0, 1 ) * error_proj( 0, 1 );
    Jth_p_proj( 1, 2 ) = Jth( 0, 1 ) * error_proj( 0, 2 );

    Jth_p_proj( 2, 0 ) = Jth( 0, 2 ) * error_proj( 0, 0 );
    Jth_p_proj( 2, 1 ) = Jth( 0, 2 ) * error_proj( 0, 1 );
    Jth_p_proj( 2, 2 ) = Jth( 0, 2 ) * error_proj( 0, 2 );

    Jth_p_proj( 3, 0 ) = Jth( 1, 3 ) * error_proj( 1, 0 );
    Jth_p_proj( 3, 1 ) = Jth( 1, 3 ) * error_proj( 1, 1 );
    Jth_p_proj( 3, 2 ) = Jth( 1, 3 ) * error_proj( 1, 2 );
                        
    Jth_p_proj( 4, 0 ) = Jth( 1, 4 ) * error_proj( 1, 0 );
    Jth_p_proj( 4, 1 ) = Jth( 1, 4 ) * error_proj( 1, 1 );
    Jth_p_proj( 4, 2 ) = Jth( 1, 4 ) * error_proj( 1, 2 );
                        
    Jth_p_proj( 5, 0 ) = Jth( 1, 5 ) * error_proj( 1, 0 );
    Jth_p_proj( 5, 1 ) = Jth( 1, 5 ) * error_proj( 1, 1 );
    Jth_p_proj( 5, 2 ) = Jth( 1, 5 ) * error_proj( 1, 2 );

    Jth_p_proj( 6, 0 ) = Jth( 2, 6 ) * error_proj( 2, 0 );
    Jth_p_proj( 6, 1 ) = Jth( 2, 6 ) * error_proj( 2, 1 );
    Jth_p_proj( 6, 2 ) = Jth( 2, 6 ) * error_proj( 2, 2 );
                        
    Jth_p_proj( 7, 0 ) = Jth( 2, 7 ) * error_proj( 2, 0 );
    Jth_p_proj( 7, 1 ) = Jth( 2, 7 ) * error_proj( 2, 1 );
    Jth_p_proj( 7, 2 ) = Jth( 2, 7 ) * error_proj( 2, 2 );
                        
    Jth_p_proj( 8, 0 ) = Jth( 2, 8 ) * error_proj( 2, 0 );
    Jth_p_proj( 8, 1 ) = Jth( 2, 8 ) * error_proj( 2, 1 );
    Jth_p_proj( 8, 2 ) = Jth( 2, 8 ) * error_proj( 2, 2 );

    Jth_p_proj( 9, 0 ) = error_proj( 0, 0 );
    Jth_p_proj( 9, 1 ) = error_proj( 0, 1 );
    Jth_p_proj( 9, 2 ) = error_proj( 0, 2 );

    Jth_p_proj( 10, 0 ) = error_proj( 1, 0 );
    Jth_p_proj( 10, 1 ) = error_proj( 1, 1 );
    Jth_p_proj( 10, 2 ) = error_proj( 1, 2 );

    Jth_p_proj( 11, 0 ) = error_proj( 2, 0 );
    Jth_p_proj( 11, 1 ) = error_proj( 2, 1 );
    Jth_p_proj( 11, 2 ) = error_proj( 2, 2 );

    for (unsigned int i = 0; i < dof; ++i)
      for (unsigned int j = i; j < dof; ++j) {
        double accum = Jth_p_proj[i][0] * Jth[0][j];
        for (unsigned int k = 1; k < dim; ++k)
          accum += Jth_p_proj[i][k] * Jth[k][j];
        AtAs[ThreadId][i][j] += accum;
        //AtA[i][j] = AtA[j][i] += accum;
      }
    for (unsigned int i = 0; i < dof; ++i) {
      double accum = Jth_p_proj[i][0] * p->location_[0];
      for (unsigned int k = 1; k < dim; ++k)
        accum += Jth_p_proj[i][k] * p->location_[k];
      Atbs[ThreadId][i][0] += accum;
    }

//#endif
  }





}





template < unsigned int dim, unsigned int dof >
ITK_THREAD_RETURN_TYPE rrl_estimation_symmetric_ICP_matching_all<dim, dof>::ThreaderCallbackForward( void * arg )
{
  typedef ThreaderType::ThreadInfoStruct  ThreadInfoType;
  ThreadInfoType * infoStruct = static_cast< ThreadInfoType * >( arg );
  const unsigned int threadId = infoStruct->ThreadID;
  //std::cout << "Thread = " << threadId << std::endl;

  ThreadStruct *str;
  str = (ThreadStruct *)(((ThreaderType::ThreadInfoStruct *)(arg))->UserData);

  str->estimation->ThreadedGenerateDataForward( threadId );

  return ITK_THREAD_RETURN_VALUE;
}




template < unsigned int dim, unsigned int dof >
ITK_THREAD_RETURN_TYPE rrl_estimation_symmetric_ICP_matching_all<dim, dof>::ThreaderCallbackBackward( void * arg )
{
  typedef ThreaderType::ThreadInfoStruct  ThreadInfoType;
  ThreadInfoType * infoStruct = static_cast< ThreadInfoType * >( arg );
  const unsigned int threadId = infoStruct->ThreadID;
  //std::cout << "Thread = " << threadId << std::endl;

  ThreadStruct *str;
  str = (ThreadStruct *)(((ThreaderType::ThreadInfoStruct *)(arg))->UserData);

  str->estimation->ThreadedGenerateDataBackward( threadId );

  return ITK_THREAD_RETURN_VALUE;
}


// Estimate initial covariance matrix.
template < unsigned int dim, unsigned int dof >
void rrl_estimation_symmetric_ICP_matching_all<dim, dof>::estimate_LS( bool estimate_parameters )
{
  double forward_scale = forward_->estimate_scale_and_assign_weight();
  if( forward_scale < 0.005 ) forward_scale = 0.005; // in normalized coordinates

  double backward_scale = backward_->estimate_scale_and_assign_weight();
  if( backward_scale < 0.005 ) backward_scale = 0.005; // in normalized coordinates

  // Normalize before estimation.
  //
  // Normalization parameters.
  location_type  center_moving( 0.0 ), center_fixed( 0.0 );
  coord_type  avg_rad_moving( 0.0 ), avg_rad_fixed( 0.0 );

  if( normalize_matches_ ) {
    // convert transformation and matches back to normalized coordinates
    forward_->cdcl_normalize_matches( center_moving, avg_rad_moving, center_fixed, avg_rad_fixed );
    forward_->trans_->normalize( avg_rad_moving, avg_rad_fixed, center_moving, center_fixed );

    // normalize backward transformation and matches with flipped radii and centers
    this->cdcl_normalize_matches_known( backward_->matches_, center_fixed, avg_rad_fixed, center_moving, avg_rad_moving );
    backward_->trans_->normalize( avg_rad_fixed, avg_rad_moving, center_fixed, center_moving );
  }


  if( forward_->weight_by_strength_ ) forward_->weight_by_strength();
  if( backward_->weight_by_strength_ ) backward_->weight_by_strength();

  forward_->weight_spatially();
  backward_->weight_spatially();

  // Estimate initial matrix by setting up a least squares problem
  // A x = q
  // A^T A
  vnl_matrix_fixed< coord_type, dof, dof >  AtA( 0.0 );
  vnl_matrix_fixed< coord_type, dof, 1 >    Atb( 0.0 );

  double forward_backward_ratio = double( forward_->matches_.size() ) / double( forward_->matches_.size() + backward_->matches_.size() );
  
  // FORWARD
  {
  typedef typename vcl_vector< typename vcl_vector< feature_sptr_type >::size_type >  indices_vector_type;
  indices_vector_type  indices_vector( forward_->matches_.size(), 0 );
  for( typename indices_vector_type::size_type n = 0; n < forward_->matches_.size(); ++n ) {
    indices_vector[n] = n;
  }
  vcl_random_shuffle( indices_vector.begin(), indices_vector.end() );
  
  //vnl_matrix_fixed< coord_type, dim, 1 >  q_loc;
  // go through all moving points
  //vcl_random_shuffle( forward_->matches_.begin(), forward_->matches_.end() ); // Warning: currently can't use shuffle because we would later erase some forward matches, when deleting backward matches addition
//  if( number_matches_ >= forward_->matches_.size() ) {
//    number_matches_ = forward_->matches_.size();
//    finest_level_ = true;
//  }
  typename vcl_vector< feature_sptr_type >::size_type  matches_used = number_matches_ * forward_backward_ratio;
  //matches_used = vcl_min( matches_used, forward_->matches_.size() );
  matches_used = forward_->matches_.size();
  //for( typename vcl_vector< feature_sptr_type >::size_type  i = end_of_forward; i < forward_->matches_.size(); ++i ) {


  m_Threader->SetNumberOfThreads( 1 ); // no multithreading
  //m_Threader->SetNumberOfThreads( 2 );


  // split the vector of matches into blocks for each thread (blocks stored as indices)
  unsigned int numberOfThreads = m_Threader->GetNumberOfThreads();
  unsigned int oneBlock = matches_used / numberOfThreads;
  m_Indices.clear();
  for( unsigned int n = 0; n < numberOfThreads; ++n ) {
    m_Indices.push_back( n * oneBlock );
  }
  m_Indices.push_back( matches_used );  // the last thread will have the remainder in case matches_used is not divisible by the numberOfThreads


  // BeforeThreadedGenerateData
  AtAs.resize( numberOfThreads );
  Atbs.resize( numberOfThreads );
  for( unsigned int n = 0; n < numberOfThreads; ++n ) {
    AtAs[n].fill( 0.0 );
    Atbs[n].fill( 0.0 );
  }

  // Set up the multithreaded processing
  ThreadStruct str;
  str.estimation = this;

  // ThreadedGenerateData
  m_Threader->SetSingleMethod( this->ThreaderCallbackForward, &str );
  m_Threader->SingleMethodExecute();

//this->ThreadedGenerateDataForward( 0 );
//this->ThreadedGenerateDataForward( 1 );
//this->ThreadedGenerateDataForward( 2 );
//this->ThreadedGenerateDataForward( 3 );












  //// the above COMPUTE_ATA filled only upper triangular -> copy to the lower triangle
  //for (unsigned int i = 0; i < dof; ++i)
  //  for (unsigned int j = i+1; j < dof; ++j)
  //    (*AtA)[j][i] = (*AtA)[i][j];

  }

  // BACKWARD -- reversing the roles of moving and fixed points
  {
  typedef typename vcl_vector< typename vcl_vector< feature_sptr_type >::size_type >  indices_vector_type;
  indices_vector_type  indices_vector( backward_->matches_.size(), 0 );
  for( typename indices_vector_type::size_type n = 0; n < backward_->matches_.size(); ++n ) {
    indices_vector[n] = n;
  }
  vcl_random_shuffle( indices_vector.begin(), indices_vector.end() );
  
  //vnl_matrix_fixed< coord_type, dim, 1 >  q_loc;
  // go through all moving points
  //vcl_random_shuffle( backward_->matches_.begin(), backward_->matches_.end() ); // Warning: currently can't use shuffle because we would later erase some forward matches, when deleting backward matches addition
//  if( number_matches_ >= backward_->matches_.size() ) {
//    number_matches_ = backward_->matches_.size();
//    finest_level_ = true;
//  }
  typename vcl_vector< feature_sptr_type >::size_type  matches_used = number_matches_ * ( 1.0 - forward_backward_ratio );
  //matches_used = vcl_min( matches_used, backward_->matches_.size() );
  matches_used = backward_->matches_.size();
  //for( typename vcl_vector< feature_sptr_type >::size_type  i = end_of_forward; i < backward_->matches_.size(); ++i ) {

  // split the vector of matches into blocks for each thread (blocks stored as indices)
  unsigned int numberOfThreads = m_Threader->GetNumberOfThreads();
  unsigned int oneBlock = matches_used / numberOfThreads;
  m_Indices.clear();
  for( unsigned int n = 0; n < numberOfThreads; ++n ) {
    m_Indices.push_back( n * oneBlock );
  }
  m_Indices.push_back( matches_used );  // the last thread will have the remainder in case matches_used is not divisible by the numberOfThreads


  // Set up the multithreaded processing
  ThreadStruct str;
  str.estimation = this;

  // ThreadedGenerateData
  m_Threader->SetSingleMethod( this->ThreaderCallbackBackward, &str );
  m_Threader->SingleMethodExecute();

//this->ThreadedGenerateDataBackward( 0 );
//this->ThreadedGenerateDataBackward( 1 );
//this->ThreadedGenerateDataBackward( 2 );
//this->ThreadedGenerateDataBackward( 3 );











  // AfterThreadedGenerateData
  for( unsigned int n = 0; n < numberOfThreads; ++n ) {
    AtA += AtAs[n];
    Atb += Atbs[n];
  }

  // the above filled only upper triangular -> copy to the lower triangle
  for (unsigned int i = 0; i < dof; ++i)
    for (unsigned int j = i+1; j < dof; ++j)
      AtA[j][i] = AtA[i][j];

  }

  vnl_svd< double > svd( AtA );
  vnl_matrix_fixed< double, dof, dof >  invAA = svd.inverse();

  // don't need to multiply by scale, since it is already included in the robust weights [yang:cvpr04]
  forward_->trans_->set_covariance( invAA );

  // assign inverse covariance
  // note: would have to be estimated
  //backward_->trans_->set_covariance( );

  if( estimate_parameters ) {
    vnl_matrix_fixed< double, dof, 1 >  new_params = invAA * Atb;
    //vcl_cout << "PARAMS: " << new_params << vcl_endl;

    forward_->trans_->set_parameterization( new_params.get_column( 0 ) );

    //backward_->trans_ = forward_->trans_->inverse();
  }

  if( normalize_matches_ ) {
    // convert transformation and matches back to unnormalized coordinates
    forward_->trans_->unnormalize( avg_rad_moving, avg_rad_fixed, center_moving, center_fixed );
    forward_->cdcl_unnormalize_matches( center_moving, avg_rad_moving, center_fixed, avg_rad_fixed );

    // unnormalize backward transformation and matches with flipped radii and centers
    backward_->trans_->unnormalize( avg_rad_fixed, avg_rad_moving, center_fixed, center_moving );
    backward_->cdcl_unnormalize_matches( center_fixed, avg_rad_fixed, center_moving, avg_rad_moving );
  }


  if( estimate_parameters ) {
    backward_->trans_ = forward_->trans_->inverse();
  }


  //vcl_cout << "Initial forward covariance: " << vcl_endl << forward_->trans_->get_covariance() << vcl_endl << "Scale: " << forward_scale << vcl_endl;
  //vcl_cout << "Initial backward covariance: " << vcl_endl << backward_->trans_->get_covariance() << vcl_endl << "Scale: " << backward_scale << vcl_endl;

}


template < unsigned int dim, unsigned int dof >
void rrl_estimation_symmetric_ICP_matching_all<dim, dof>::estimate_LS_backward( bool estimate_parameters )
{
  // swap forwad and backwad, estimate, and swap back
  typename estimation_type::sptr  forward = forward_;
  typename estimation_type::sptr  backward = backward_;

  forward_ = backward_;
  backward_ = forward;
  
  this->estimate_LS( estimate_parameters );

  forward_ = forward;
  backward_ = backward;

#if 0
  trans_sptr_type  backward_global_transform;

  double forward_scale = forward_->estimate_scale_and_assign_weight();
  if( forward_scale < 0.005 ) forward_scale = 0.005; // in normalized coordinates

  double backward_scale = backward_->estimate_scale_and_assign_weight();
  if( backward_scale < 0.005 ) backward_scale = 0.005; // in normalized coordinates


  // Normalize before estimation.
  //
  // Normalization parameters.
  location_type  center_moving( 0.0 ), center_fixed( 0.0 );
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
  if( forward_->weight_by_strength_ ) forward_->weight_by_strength();

  forward_->weight_spatially();
  backward_->weight_spatially();

  // Estimate initial matrix by setting up a least squares problem
  // A x = q
  // A^T A
  vnl_matrix_fixed< double, dof, dof >  AtA( 0.0 );
  vnl_matrix_fixed< double, dof, 1 >    Atb( 0.0 );
  
  // BACKWARD
  {
  vnl_matrix_fixed< double, dim, 1 >  q_loc;
  // go through all moving points
  //for( typename vcl_vector< feature_sptr_type >::size_type  i = end_of_backward; i < backward_->matches_.size(); ++i ) {
  for( typename vcl_vector< feature_sptr_type >::size_type  i = 0; i < backward_->matches_.size(); ++i ) {
    // moving point p
    feature_sptr_type  p = backward_->matches_[i]->from_;
  
    // Jacobian w.r.t. parameters
    vnl_matrix_fixed< double, dim, dof > const &  Jth = backward_->trans_->jacobian_wrt_par( p->location_ );
    //vcl_cout << p->location_ << "  pc  " << trans_->center_moving_ << vcl_endl;

    // fixed point q
    feature_sptr_type  q = backward_->matches_[i]->to_[0];
    q_loc.set_column( 0, q->location_ );

    double weight = backward_->matches_[i]->w_[0];

    AtA += Jth.transpose() * q->error_projector_ * Jth * weight;
    Atb += Jth.transpose() * q->error_projector_ * q_loc * weight;

    //vcl_cout << Jth << vcl_endl << AtA << vcl_endl << vcl_endl;

  }
  }

  // FORWARD
  {
  typedef typename vcl_vector< typename vcl_vector< feature_sptr_type >::size_type >  indices_vector_type;
  indices_vector_type  indices_vector( forward_->matches_.size(), 0 );
  for( typename indices_vector_type::size_type n = 0; n < forward_->matches_.size(); ++n ) {
    indices_vector[n] = n;
  }
  vcl_random_shuffle( indices_vector.begin(), indices_vector.end() );
  
  vnl_matrix_fixed< coord_type, dim, 1 >  q_loc;
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
  
    // fixed point q
    feature_sptr_type  q = forward_->matches_[ind]->to_[0];

    // Jacobian w.r.t. parameters
    vnl_matrix_fixed< coord_type, dim, dof > const &  Jth = backward_->trans_->jacobian_wrt_par( q->location_ );
    //vcl_cout << p->location_ << "  pc  " << trans_->center_moving_ << vcl_endl;


    q_loc.set_column( 0, p->location_ );

    double weight = forward_->matches_[ind]->w_[0];

    AtA += Jth.transpose() * p->error_projector_ * Jth * weight;
    Atb += Jth.transpose() * p->error_projector_ * q_loc * weight;

    //vcl_cout << Jth << vcl_endl << AtA << vcl_endl << vcl_endl;

  }
  }

  vnl_svd< coord_type > svd( AtA );
  vnl_matrix_fixed< double, dof, dof >  invAA = svd.inverse();

  // don't need to multiply by scale, since it is already included in the robust weights [yang:cvpr04]
  backward_->trans_->set_covariance( invAA );

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
#endif
}


// Normalize moving and fixed points which the matches were formed with
// Use supplied known centers and radii
template < unsigned int dim, unsigned int dof >
void rrl_estimation_symmetric_ICP_matching_all<dim, dof>::cdcl_normalize_matches_known( vcl_vector< match_sptr_type >        & matches,
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
//    p->covariance_ /= (avg_radius_moving*avg_radius_moving);

    // fixed point q
    feature_sptr_type  q = matches[i]->to_[0];
    // if it was inserted (i.e. not already there)
    if( fixed_matched.insert( q ).second ) {
      // convert locations to normalized coordinate system
      q->location_ -= center_fixed;
      q->location_ /= avg_radius_fixed;

      // convert covariances to normalized coordinate system
//        q->covariance_ /= (avg_radius_fixed*avg_radius_fixed);
    }
  }

}



#define RRL_ESTIMATION_SYMMETRIC_ICP_MATCHING_ALL_INSTANTIATE( dim, dof ) \
template                                                  \
class rrl_estimation_symmetric_ICP_matching_all< dim, dof >;
