#include <vcl_map.h>
#include <vcl_utility.h>
#include <vcl_algorithm.h>
#include <vcl_iomanip.h>
#include <vcl_iostream.h>
#include <vcl_iterator.h>
#include <vcl_cstdlib.h>

#include <vul/vul_timer.h>
#include <vnl/algo/vnl_svd.h>
#include <vnl/algo/vnl_cholesky.h>
#include <vnl/vnl_inverse.h>
#include <vnl/vnl_transpose.h>

#include <vnl/algo/vnl_symmetric_eigensystem.h>
#include <vnl/algo/vnl_generalized_eigensystem.h>

#include "cdcl_estimation.h"
#include "cdcl_utils.h"
#include "cdcl_trans_rigid3d.h"
#include "cdcl_lbfgs.h"

//#undef VTK_FOUND
// VTK is used to dump data for analysis in external programs.
#ifdef VTK_FOUND
#include "cdcl_utils_VTK.h"
#endif

//:
// \file
// \brief  Estimation handling the whole registration: 1) matching,
//         2) parameter estimate, 3) covariance estimate.
//         Feature location and parameter covariances are modeled.
//         The estimation is done with expectaction maximization and covariance
//         estimated as a free parameter. Derived from ICP since it's used
//         for initialization.
// \author Michal Sofka
// \date   May 2006


// Construct estimation object given moving and fixed feature sets
// and an initial transformation.
template < unsigned int dim, unsigned int dof >
cdcl_estimation<dim, dof>::cdcl_estimation( vcl_vector< feature_sptr_type > const &  moving, 
                                            vcl_vector< feature_sptr_type > const &  fixed,
                                            trans_sptr_type                 const &  trans )
  : cdcl_estimation_ICP<dim,dof>( moving, fixed, trans )
{
}


// Construct estimation object given moving and fixed feature sets
// and an initial transformation.
template < unsigned int dim, unsigned int dof >
cdcl_estimation<dim, dof>::cdcl_estimation( vcl_vector< vcl_vector< feature_sptr_type > > const &  moving, 
                                            vcl_vector< vcl_vector< feature_sptr_type > > const &  fixed,
                                            trans_sptr_type                               const &  trans,
                                            vcl_vector< double >                          const &  fixed_spacing )
  : cdcl_estimation_ICP<dim,dof>( moving, fixed, trans, fixed_spacing )
{ 
}


#define PROGRESS_OUTPUT( procedure, timer )  vcl_cout << "----> Computing " << procedure << " ... " << vcl_endl; timer.mark()
#define TIMER_OUTPUT( procedure, timer )  vcl_cout << "----> Computing " << procedure << " done in " << timer.real()/1000 << "." << timer.real()%1000 << " sec." << vcl_endl


// Initialize estimation (compute initial weights and covariance).
template < unsigned int dim, unsigned int dof >
void cdcl_estimation<dim, dof>::initialize()
{
  vul_timer timer;

  PROGRESS_OUTPUT( "initial weights", timer );
  this->find_closest_euclidean();
  TIMER_OUTPUT( "initial weights", timer );

  // estimate initial covariance
  PROGRESS_OUTPUT( "initial covariance", timer );
  this->estimate_LS();
  TIMER_OUTPUT( "initial covariance", timer );

  this->trans_->print( vcl_cout );
}


// Run one EM iteration of estimation (parameters, weights, covariance, weights).
// Return true when converged.
template < unsigned int dim, unsigned int dof >
bool cdcl_estimation<dim, dof>::one_iteration( void *  caller, 
                                               display_callback_type  display_points,
                                               unsigned int  iteration )
{
  vul_timer timer;

  unsigned int startEM = 5;

  vcl_cout << "******************** Iteration " << iteration << "********************" << vcl_endl;

  bool resolution_switch_allowed = iteration > startEM;
  PROGRESS_OUTPUT( "weights", timer );
  this->compute_weights( resolution_switch_allowed );
  TIMER_OUTPUT( "weights", timer );
  
  if( caller ) display_points( caller, this->moving_[this->res_], this->fixed_[this->res_], this->matches_, this->trans_, iteration );
  #ifdef VTK_FOUND
  cdcl_write_matches_VTK( this->moving_[this->res_], this->fixed_[this->res_], this->matches_, this->trans_, iteration );
  #endif

  trans_sptr_type  global_transform;
  if( this->trans_->is_type( cdcl_trans_rigid3d::type_id() ) ) { // if transform is rigid with small angle approximation
    // from now on, transform is incremental
    // save current transform for recomposing
    global_transform = this->trans_;
    this->trans_ = this->trans_->create_new(); // will create identity transformation

    this->trans_->set_covariance( global_transform->get_covariance() );

    // for rigid transform, we are estimating incremental transform
    // therefore transform points first
    vnl_matrix_fixed< double, dim, dim > const &  Jp = global_transform->jacobian_wrt_loc();
    // go through all moving points
    for( typename vcl_vector< match_sptr_type >::size_type  i = 0; i < this->matches_.size(); ++i ) {
      // moving point p
      feature_sptr_type  p = this->matches_[i]->from_;
      
      // reset the original from point
      // transform location and covariance
      const vnl_vector_fixed< double, dim >      location   = global_transform->map_loc( p->location_ );
      vnl_matrix_fixed< double, dim, dim > covariance = Jp * p->covariance_ * Jp.transpose();
      this->matches_[i]->from_ = new cdcl_feature< dim >( location, covariance );
    }
    this->trans_->set_parameterization( vnl_vector< double >( 6, 0.0 ) );
  }

  bool parameters_converged = false;
  bool covariance_converged = false;
  if( iteration > startEM ) {
    PROGRESS_OUTPUT( "parameters", timer );
    parameters_converged = this->estimate_parameters();
    TIMER_OUTPUT( "parameters", timer );

    PROGRESS_OUTPUT( "weights", timer );
    if( this->trans_->is_type( cdcl_trans_rigid3d::type_id() ) ) { // if transform is rigid with small angle approximation
      // new_transform is incremental rigid transform
      global_transform->recompose_increment( this->trans_ );
      this->compute_weights( resolution_switch_allowed, global_transform );


      // for rigid transform, we are estimating incremental transform
      // therefore transform points first
      vnl_matrix_fixed< double, dim, dim > const &  Jp = global_transform->jacobian_wrt_loc();
      // go through all moving points
      for( typename vcl_vector< match_sptr_type >::size_type  i = 0; i < this->matches_.size(); ++i ) {
        // moving point p
        feature_sptr_type  p = this->matches_[i]->from_;
        
        // reset the original from point
        // transform location and covariance
        const vnl_vector_fixed< double, dim >      location   = global_transform->map_loc( p->location_ );
        vnl_matrix_fixed< double, dim, dim > covariance = Jp * p->covariance_ * Jp.transpose();
        this->matches_[i]->from_ = new cdcl_feature< dim >( location, covariance );
      }

    }
    else this->compute_weights( resolution_switch_allowed );
    TIMER_OUTPUT( "weights", timer );
  }

  PROGRESS_OUTPUT( "covariance", timer );
  covariance_converged = this->estimate_covariance();
  TIMER_OUTPUT( "covariance", timer );

  if( this->trans_->is_type( cdcl_trans_rigid3d::type_id() ) ) { // if transform is rigid with small angle approximation
    global_transform->set_covariance( this->trans_->get_covariance() );
    this->trans_ = global_transform;
    // from now on, transform is not incremental
  }

  this->trans_->print( vcl_cout );

  //vcl_cout << "CONVERGENCE: " << parameters_converged << " " << covariance_converged << " " << trans_->get_covariance().frobenius_norm() << vcl_cout;
  return parameters_converged && covariance_converged && this->trans_->get_covariance().frobenius_norm() < 1e-3;
}


// Run estimation.
template < unsigned int dim, unsigned int dof >
bool cdcl_estimation<dim, dof>::run( void *  caller, 
                                     display_callback_type  display_points )
{
  this->initialize();

  bool converged = false;
  unsigned int max_iterations = 85;
  unsigned int iteration = 0;

  while( !converged && iteration < max_iterations ) {
    converged = this->one_iteration( caller, display_points, iteration );  
    ++iteration;
  }

  return converged;
}


// Estimate parameters.
// Return true when converged.
template < unsigned int dim, unsigned int dof >
bool cdcl_estimation<dim, dof>::estimate_parameters()
{
  // Normalize before estimation.
  //
  // Normalization parameters.
  vnl_vector_fixed< double, dim >  center_moving( 0.0 ), center_fixed( 0.0 );
  double avg_rad_moving( 0.0 ), avg_rad_fixed( 0.0 );

  if( this->normalize_matches_ ) {
    this->cdcl_normalize_matches( center_moving, avg_rad_moving, center_fixed, avg_rad_fixed );
    this->trans_->normalize( avg_rad_moving, avg_rad_fixed, center_moving, center_fixed );
  }

  // objective function
  obj_fun_par_type  f_par( this->matches_, this->trans_ );

  vnl_vector< double >  x = this->trans_->get_parameterization();  // can afford to copy the vector
  vnl_vector< double >  x_old = x;

  // estimator
  cdcl_lbfgs  minimizer( f_par );
  minimizer.set_max_function_evals( 5 );

  //minimizer.set_f_tolerance( 0.001 );
  //minimizer.set_x_tolerance( 1e-3 );
  minimizer.set_verbose( true );
  minimizer.set_check_derivatives( false );
  minimizer.set_trace( false );
  //minimizer.default_step_length = 0.1;
  bool converged = minimizer.minimize( x );
  bool parameters_converged = ( minimizer.get_start_error() - minimizer.get_end_error() ) < 1e-4;

  if( ( !converged && minimizer.get_failure_code() < 0 ) || x.size() == 0 ) {
    vcl_cerr << "Error during minimization, failure code: " << minimizer.get_failure_code() << " . Keeping the old value." << vcl_endl;
    this->trans_->set_parameterization( x_old );
    //exit( 1 );
  }
  else {
    this->trans_->set_parameterization( x );
  }

  if( this->normalize_matches_ ) {
    // convert transformation back to unnormalized coordinates
    this->trans_->unnormalize( avg_rad_moving, avg_rad_fixed, center_moving, center_fixed );
    this->cdcl_unnormalize_matches( center_moving, avg_rad_moving, center_fixed, avg_rad_fixed );
  }

  return parameters_converged;
}


// Estimate covariance.
// Return true when converged.
template < unsigned int dim, unsigned int dof >
bool cdcl_estimation<dim, dof>::estimate_covariance()
{
  // Normalize before estimation.
  //
  // Normalization parameters.
  vnl_vector_fixed< double, dim >  center_moving( 0.0 ), center_fixed( 0.0 );
  double avg_rad_moving( 0.0 ), avg_rad_fixed( 0.0 );

  if( this->normalize_matches_ ) {
    this->cdcl_normalize_matches( center_moving, avg_rad_moving, center_fixed, avg_rad_fixed );
    this->trans_->normalize( avg_rad_moving, avg_rad_fixed, center_moving, center_fixed );
  }

  // objective function
  obj_fun_cov_type  f_cov( this->matches_, this->trans_ );

  // number of independent parameters in the covariance matrix
  const unsigned int cov_dof = dof*(dof+1)/2;

  // assign upper triangular matrix, cholesky decomposition of the covariance, to x
  vnl_cholesky  chol( this->trans_->get_covariance() );
  vnl_matrix< double >  upper = chol.upper_triangle();
  vnl_vector< double >  x( cov_dof, 0.0 );
  unsigned int ind = 0;
  for( unsigned int d1 = 0; d1 < dof; ++d1 )
    for( unsigned int d2 = d1; d2 < dof; ++d2 ) {
      x( ind ) = upper( d1, d2 );
      //vcl_cout << d1 << " " << d2 << " " << upper( d1, d2 ) << vcl_endl;
      ++ind;
    }

  vnl_vector< double >  x_old = x;

  // estimator
  cdcl_lbfgs  minimizer( f_cov );
  minimizer.set_max_function_evals( 5 );
  minimizer.set_verbose( true );
  minimizer.set_check_derivatives( false );
  minimizer.set_trace( false );
  //minimizer.default_step_length = 0.1;
  bool converged = minimizer.minimize( x );
  bool covariance_converged = ( minimizer.get_start_error() - minimizer.get_end_error() ) < 1e-5;

  if( ( !converged && minimizer.get_failure_code() < 0 ) || x.size() == 0 ) {
      vcl_cerr << "Error during minimization, failure code: " << minimizer.get_failure_code() << " . Keeping the old value." << vcl_endl;
      x = x_old;
      //exit( 1 );
  }

  // assign parameterization into upper triangular matrix
  upper.fill( 0.0 );
  ind = 0;
  for( unsigned int d1 = 0; d1 < dof; ++d1 )
    for( unsigned int d2 = d1; d2 < dof; ++d2 ) {
      upper( d1, d2 ) = x( ind );
      ++ind;
    }

  // reconstruct the covariance from the cholesky decomposition of the covariance
  this->trans_->set_covariance( vnl_transpose( upper )*upper );

  //if( !minimized ) 
  //  vcl_abort();
  // covariance has already been reconstructed from the cholesky decomposition in the optimizer loop


  if( this->normalize_matches_ ) {
    // convert transformation back to unnormalized coordinates
    this->trans_->unnormalize( avg_rad_moving, avg_rad_fixed, center_moving, center_fixed );
    this->cdcl_unnormalize_matches( center_moving, avg_rad_moving, center_fixed, avg_rad_fixed );
  }

  return covariance_converged;
}


// Compute weights, use Mahalanobis distances.
template < unsigned int dim, unsigned int dof >
void cdcl_estimation<dim, dof>::compute_weights( bool resolution_switch_allowed, trans_sptr_type  trans )
{
  if( !trans ) trans = this->trans_;
  // Jacobian w.r.t. p (moving point)
  vnl_matrix_fixed< double, dim, dim > const &  Jp = trans->jacobian_wrt_loc();

  // multiplier on the largest eigenvalue to get the radius
  const double eval_multiplier = 5.0;

  // num matches randomly selected in each region
  const unsigned int max_matches_allowed = 100;
 
  // save all weights for normalization
  vcl_map< feature_sptr_type, double >  to_weights;
  vcl_vector< double > from_weights;

  // query point for kd-tree
  rsdl_point  query_pt( dim, 0 );

  bool stop = false;
  bool not_restarted = true;
  while( !stop ) {
    unsigned int min_matches = 99999;
    unsigned int max_matches = 0;
    vcl_vector< unsigned int >  matches_histogram;
    matches_histogram.resize( this->fixed_[this->res_].size()/10 + 1, 0 );
    vcl_vector< unsigned int > num_matches;
    double avg_largest_eval = 0.0;
    //vnl_matrix< double >  V( dim, dim );
    //vnl_vector< double >  Dp( dim ), Dq( dim );  
    this->matches_.clear();
    // go through all moving points
    for( typename vcl_vector< feature_sptr_type >::size_type  i = 0; i < this->moving_[this->res_].size(); ++i ) {
      // moving point p
      feature_sptr_type  p = this->moving_[this->res_][i];
      vnl_vector_fixed< double, dim >  mapped_moving = trans->map_loc( p->location_ );

      // Jacobian w.r.t. transformation parameters
      vnl_matrix_fixed< double, dim, dof > const &  Jth = trans->jacobian_wrt_par( p->location_ );

      // compute largest eigen value of the transfer error covariance
      vnl_matrix_fixed< double, dim, dim >  Cij_p = this->transfer_covar( Jth );
      avg_largest_eval += vcl_sqrt( cdcl_largest_evalue( Cij_p ) );

      // Cij_p covariance of the match without fixed feature covariance
      // this will serve as upper bound of our search range
      Cij_p += Jp * p->covariance_ * Jp.transpose();

      // select largest eigenvalue as radius
      // compute largest eigenvalue efficiently
      double radius = eval_multiplier*vcl_sqrt( cdcl_largest_evalue( Cij_p ) );

      // query kd-tree and find set of points that are potential matches
      query_pt.set_cartesian( mapped_moving );     
      vcl_vector< rsdl_point >  points;
      vcl_vector< int >  indices;
      this->kd_tree_[this->res_]->points_in_radius( query_pt, radius, points, indices );

      // sum of weights for normalization
      double sum_w = 0.0;
      if( indices.size() > 0 ) {
        // select randomly point from the set of retrieved vertices
        vcl_vector< int >  indices_pruned;
        indices_pruned.clear();
        // random shuffle the indices and then pick first max_matches_possible
        vcl_random_shuffle( indices.begin(), indices.end() );
        typename vcl_vector< int >::size_type  max_matches_possible = vcl_min( max_matches_allowed, (unsigned int)indices.size() );
        
        match_sptr_type  match = new cdcl_match<dim>;
        match->from_ = this->moving_[this->res_][i];

        // go through all potential points (within search range given by the covariance)
        for( typename vcl_vector< int >::size_type  j = 0; j < max_matches_possible; ++j ) {
          // fixed point q
          feature_sptr_type  q = this->fixed_[this->res_][indices[j]];

          // matching error
          vnl_vector_fixed< double, dim >  const &  err_vec = mapped_moving - q->location_;

          // covariance of a correspondence match
          vnl_matrix_fixed< double, dim, dim >  Cij = Cij_p + q->covariance_;
          vnl_matrix_fixed< double, dim, dim >  Cij_1 = vnl_inverse( Cij );

          vnl_matrix_fixed< double, dim, 1 >  err;
          err.set_column( 0, err_vec );

          vnl_matrix_fixed< double, 1, 1 >  residual = err.transpose() * Cij_1 * err;
          double w0 = weight( residual( 0, 0 ) );

          if( w0 > 0.0 ) {     
            //// similarity weight by the dot product of the normals
            ////
            //const vnl_matrix< double >  Cq = q->covariance_;
            //vnl_symmetric_eigensystem_compute( Cq, V, Dq );

            //// zeroth column is eigenvector corresponding to smallest eigen value
            //// this will be the point normal
            //vnl_vector< double >  q_normal = V.get_column( 0 );

            //const vnl_matrix< double >  Cp = Jp * p->covariance_ * vnl_tranpose( Jp );
            //vnl_symmetric_eigensystem_compute( Cp, V, Dp );

            //// zeroth column is eigenvector corresponding to smallest eigen value
            //// this will be the point normal
            //vnl_vector< double >  p_normal = V.get_column( 0 );

            //w0 = w0 * dot_product( p_normal, q_normal );//(Dq[0]/Dq[1]) / (Dp[0]/Dp[1]);

            //// similarity weight by the similarity between the covariances
            ////
            //const vnl_matrix< double >  Cq = q->covariance_;// + transfer_covar( Jth );
            //const vnl_matrix< double >  Cp = Jp * p->covariance_ * vnl_tranpose( Jp );//Cij_p;
            //vnl_generalized_eigensystem  geig( Cp, Cq );
            //double D0 = vcl_log( geig.D( 0, 0 ) ) * vcl_log( geig.D( 0, 0 ) );
            //double D1 = vcl_log( geig.D( 1, 1 ) ) * vcl_log( geig.D( 1, 1 ) );
            //w0 = w0 * (vcl_sqrt( D0 + D1 ) );


            match->to_.push_back( q );
            match->w_.push_back( w0 );

            typename vcl_map< feature_sptr_type, double >::iterator  it = to_weights.find( q );
            if( it == to_weights.end() ) to_weights[q] = w0;
            else to_weights[q] += w0;

            sum_w += w0;
          }
        }

        if( match->to_.size() > 0 ) {
          //vcl_cout << "rad: " << radius << " ind: " << indices.size() << " to: " << match->to_.size() << vcl_endl;
          this->matches_.push_back( match );
          from_weights.push_back( sum_w );
          sum_w = 0.0;

          if( match->to_.size() < min_matches ) min_matches = match->to_.size();
          if( match->to_.size() > max_matches ) max_matches = match->to_.size();
          ++(matches_histogram[match->to_.size()/10]);
          num_matches.push_back( match->to_.size() );
        }
      }
     
    }

    avg_largest_eval /= this->matches_.size();
    
    //vcl_vector< unsigned int >::iterator  loc = num_matches.begin() + (num_matches.size()-1)/2 + 1;
    //vcl_nth_element( num_matches.begin(), loc, num_matches.end() );
    //unsigned int med_matches = *loc;


    //vcl_cout << "Switch resol, curr=" << this->res_ 
    //         << ",min=" << min_matches 
    //         << ",max=" << max_matches << vcl_endl;
    //vcl_cout << "Switch eval, spc=" << this->fixed_spacing_[this->res_]
    //         << ",avg lg=" << avg_largest_eval << vcl_endl;

    //vcl_cout << "Switch histogram: ";
    //vcl_copy( matches_histogram.begin(), matches_histogram.end(), vcl_ostream_iterator<unsigned int>(vcl_cout, " ") );
    //vcl_cout << vcl_endl;

    stop = true;
//    if( avg_largest_eval < fixed_spacing_[res_]/1.0 && not_restarted ) {
    if( avg_largest_eval < 2.0*this->fixed_spacing_[this->res_] && resolution_switch_allowed & not_restarted ) {
      if( this->res_ > 0 ) {
        --this->res_;
        stop = false;
        not_restarted = false;
      }
    }

    if( this->res_+1 <= this->max_res_ && avg_largest_eval > 2.0*this->fixed_spacing_[this->res_+1] && resolution_switch_allowed & not_restarted ) {
      ++this->res_;
      stop = false;
      not_restarted = false;
    }



  }
    
  // normalize weights
  // weight between each two points is normalized by
  // the sum of weights of all correspondences with the from point
  // and by the sum of weights of all correpondences with the to point
  vcl_vector< vcl_vector< double > >  norm_weights;
  for( typename vcl_vector< match_sptr_type >::size_type  i = 0; i < this->matches_.size(); ++i ) {
    // sum of weights of all correspondences with the from point
    double wi = from_weights[i];

    for( typename vcl_vector< match_sptr_type >::size_type  j = 0; j < this->matches_[i]->to_.size(); ++j ) {
      // sum of weights of all correpondences with the to point
      double wj = to_weights[this->matches_[i]->to_[j]];

      this->matches_[i]->w_[j] = this->matches_[i]->w_[j] * this->matches_[i]->w_[j] / ( wi * wj );
      //matches_[i]->w_[j] = matches_[i]->w_[j] / wi;
    }
  }

  //this->prune_matches( 20 );

  // covariance is probably too small for any matches too be found, use closest points then
  if( this->matches_.size() == 0 ) this->find_closest_euclidean();
}


namespace {
template <unsigned int dim>
class bigger_on_weight
{
public:
  bool operator() ( const vcl_pair< double, typename cdcl_feature< dim >::sptr >&  left, 
                    const vcl_pair< double, typename cdcl_feature< dim >::sptr >&  right )
  {
    return ( left.first > right.first );
  }
};
}


// Prune matches, keep only leave_only strongest.
template < unsigned int dim, unsigned int dof >
void cdcl_estimation<dim, dof>::prune_matches( unsigned int leave_only )
{
  for( typename vcl_vector< match_sptr_type >::size_type  i = 0; i < this->matches_.size(); ++i ) {
    unsigned int keep = vcl_min( leave_only, (unsigned int)this->matches_[i]->w_.size() );
    // pairs of weights and matched features
    vcl_vector< vcl_pair< double, feature_sptr_type > >  pairs;
    for( typename vcl_vector< match_sptr_type >::size_type  j = 0; j < this->matches_[i]->w_.size(); ++j ) {
      vcl_pair< double, feature_sptr_type >  pair( this->matches_[i]->w_[j], this->matches_[i]->to_[j] );
      pairs.push_back( pair );
    }
    typedef typename vcl_vector< vcl_pair< double, feature_sptr_type > >::iterator  pair_iterator_type;
    pair_iterator_type  nth = pairs.begin() + keep;
    vcl_nth_element( pairs.begin(), nth, pairs.end(), bigger_on_weight<dim>() );

    //// renormalize to sum 1
    //double sum_w = 0.0;
    double sum_w = 1.0;
    //for( vcl_vector< vcl_pair< double, feature_sptr_type > >::iterator  it = pairs.begin(); it != nth; ++it ) {
    //  sum_w += it->first;
    //}

    // assign new matches
    this->matches_[i]->w_.clear();
    this->matches_[i]->to_.clear();
    for( pair_iterator_type  it = pairs.begin(); it != nth; ++it ) {
      this->matches_[i]->w_.push_back( it->first / sum_w );
      this->matches_[i]->to_.push_back( it->second );
    }
  }
}



#define CDCL_ESTIMATION_INSTANTIATE( dim, dof ) \
template                                        \
class cdcl_estimation< dim, dof >;
