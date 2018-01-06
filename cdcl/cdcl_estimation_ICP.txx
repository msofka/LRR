#include <vcl_utility.h>
#include <vcl_algorithm.h>
#include <vcl_iomanip.h>
#include <vcl_iostream.h>
#include <vcl_iterator.h>
#include <vcl_cstdlib.h>
#include <vcl_memory.h>
#include <vcl_set.h>

#include <vul/vul_timer.h>

#include <vnl/vnl_inverse.h>
#include <vnl/vnl_transpose.h>
#include <vnl/algo/vnl_svd.h>
#include <vnl/algo/vnl_symmetric_eigensystem.h>
#include <vnl/algo/vnl_generalized_eigensystem.h>

#include <rrel/rrel_muset_obj.h>

#include "cdcl_estimation_ICP.h"
#include "cdcl_utils.h"
#include "cdcl_trans_rigid3d.h"

#undef VTK_FOUND
// VTK is used to dump data for analysis in external programs.
#ifdef VTK_FOUND
#include "cdcl_utils_VTK.h"
#endif

//:
// \file
// \brief  Estimation handling the whole registration: 1) matching,
//         2) parameter estimate, 3) covariance estimate.
//         ICP with Euclidean distances for matching and normal distances
//         for estimation is used.
// \author Michal Sofka
// \date   May 2006


// Construct estimation object given moving and fixed feature sets
// and an initial transformation.
template < unsigned int dim, unsigned int dof >
cdcl_estimation_ICP<dim, dof>::cdcl_estimation_ICP( vcl_vector< feature_sptr_type > const &  moving, 
                                                    vcl_vector< feature_sptr_type > const &  fixed,
                                                    trans_sptr_type                 const &  trans )
  : trans_( trans ),
    normalize_matches_( false ),
    weight_by_strength_( false ),
    scale_estimated_( false )
{
  if( dim == 3 ) { // multiresolution for 3d registration
    // remove constness
    vcl_vector< typename cdcl_feature< dim >::sptr >  fixed_fea = fixed;
    vcl_vector< typename cdcl_feature< dim >::sptr >  moving_fea = moving;
    vnl_vector_fixed< double, dim >  center_fixed, center_moving;
    double  avg_rad_fixed, avg_rad_moving;

    cdcl_normalize_data( moving_fea, center_moving, avg_rad_moving ); // data need to be normalized for computing covariances inside subsample (or pass in centers)

    if( trans->is_type( cdcl_trans_rigid3d::type_id() ) ) {
      // RIGID
      // Same radius for fixed and moving
      avg_rad_fixed = avg_rad_moving;
      cdcl_normalize_data_known_radius( fixed_fea, center_fixed, avg_rad_fixed );
    }
    else {
      // AFFINE
      cdcl_normalize_data( fixed_fea, center_fixed, avg_rad_fixed );
    }

    // normalize initial transform
    trans_->normalize( avg_rad_moving, avg_rad_fixed, center_moving, center_fixed );

    vcl_vector< double > dummy;
    subsample_data_fine_covariances( moving_fea, moving_, dummy );
    subsample_data_fine_covariances( fixed_fea, fixed_, fixed_spacing_ );

    //cdcl_unnormalize_data( fixed_fea, center_fixed, avg_rad_moving );
    //cdcl_unnormalize_data( moving_fea, center_moving, avg_rad_moving );
  }
  else {
    moving_.push_back( moving );
    fixed_.push_back( fixed );
    fixed_spacing_.push_back( 1.0 );  // this is for switching resolutions and gets set in subsample
  }
  
  this->build_kd_trees();
}


// Construct estimation object given moving and fixed feature sets
// and an initial transformation.
template < unsigned int dim, unsigned int dof >
cdcl_estimation_ICP<dim, dof>::cdcl_estimation_ICP( vcl_vector< vcl_vector< feature_sptr_type > > const &  moving, 
                                                    vcl_vector< vcl_vector< feature_sptr_type > > const &  fixed,
                                                    trans_sptr_type                               const &  trans,
                                                    vcl_vector< double >                          const &  fixed_spacing )
  : moving_( moving ),
    fixed_( fixed ),
    fixed_spacing_( fixed_spacing ),
    trans_( trans ),
    normalize_matches_( false ),
    weight_by_strength_( false ),
    scale_estimated_( false )
{ 
  this->build_kd_trees();
}


// Build kd trees, one for each resolution level of the fixed image.
template < unsigned int dim, unsigned int dof >
void cdcl_estimation_ICP<dim, dof>::build_kd_trees()
{
  //save_data( moving_, "moving" );
  //save_data( fixed_, "fixed" );

  // set the coarsest resolution
  res_ = max_res_ = vcl_min( moving_.size(), fixed_.size() ) - 1;


  // build a set of kd_tree's, one for each resolution level of the fixed image
  for( typename vcl_vector< vcl_vector< feature_sptr_type > >::size_type n = 0; n < fixed_.size(); ++n ) {
    // build kd_tree from fixed points
    vcl_vector< rsdl_point > points;
    rsdl_point  pt( dim, 0 );
    for( typename vcl_vector< feature_sptr_type >::size_type j = 0; j < fixed_[n].size(); ++j ) {
      pt.set_cartesian( fixed_[n][j]->location_ );
      points.push_back( pt );
    }

    vcl_cout << "kdtree size: " << points.size() << vcl_endl;
    kd_tree_.push_back( new rsdl_kd_tree( points ) );
  }
}

// Compute median error using current matches.
template < unsigned int dim, unsigned int dof >
double
cdcl_estimation_ICP<dim, dof>::median_error()
{
  // for decomposition to get the eigenvectors from the covariances
  vnl_matrix< double >  V( dim, dim );
  vnl_vector< double >  D( dim );

  vcl_vector< double >  errors;
  // go through all moving points
  for( typename vcl_vector< match_sptr_type >::size_type  i = 0; i < matches_.size(); ++i ) {
    // moving point p
    feature_sptr_type  p = matches_[i]->from_;   
    vnl_vector_fixed< double, dim >  mapped_moving = trans_->map_loc( p->location_ );

    // go through all matches (fixed points)
    for( typename vcl_vector< feature_sptr_type >::size_type  j = 0; j < matches_[i]->to_.size(); ++j ) {
      // fixed point q
      feature_sptr_type  q = matches_[i]->to_[j];

      // matching error
      vnl_vector_fixed< double, dim >  const &  err_vec = mapped_moving - q->location_;
      vnl_matrix_fixed< double, dim, 1 >  err;
      err.set_column( 0, err_vec );

      // normal distance
      //
      const vnl_matrix< double >  A = q->covariance_;
      vnl_symmetric_eigensystem_compute( A, V, D );

      // zeroth column is eigenvector corresponding to smallest eigen value
      // this will be the point normal
      vnl_vector< double >  evect_sm = V.get_column( 0 );

      // normal distance error projector
      vnl_matrix_fixed< double, dim, dim >  err_proj = outer_product( evect_sm, evect_sm );

      vnl_matrix_fixed< double, 1, 1 >  residual = err.transpose() * err_proj * err;

      errors.push_back( vcl_sqrt( residual( 0, 0 ) ) );

    }
  }

  vcl_vector< double >::iterator  loc = errors.begin() + (errors.size()-1)/2 + 1;
  vcl_nth_element( errors.begin(), loc, errors.end() );
  double med_error = *loc;

  return med_error;
}


// Compute RMS error using current matches.
template < unsigned int dim, unsigned int dof >
double
cdcl_estimation_ICP<dim, dof>::RMS_error()
{
  // for decomposition to get the eigenvectors from the covariances
  vnl_matrix< double >  V( dim, dim );
  vnl_vector< double >  D( dim );
  vnl_matrix_fixed< double, dim, dim >  Jp = trans_->jacobian_wrt_loc();

  double error = 0.0;
  unsigned int count = 0;
  // go through all moving points
  for( typename vcl_vector< match_sptr_type >::size_type  i = 0; i < matches_.size(); ++i ) {
    // moving point p
    feature_sptr_type  p = matches_[i]->from_;   
    vnl_vector_fixed< double, dim >  mapped_moving = trans_->map_loc( p->location_ );

    // go through all matches (fixed points)
    for( typename vcl_vector< feature_sptr_type >::size_type  j = 0; j < matches_[i]->to_.size(); ++j ) {
      // fixed point q
      feature_sptr_type  q = matches_[i]->to_[j];

      // matching error
      vnl_vector_fixed< double, dim >  const &  err_vec = mapped_moving - q->location_;
      vnl_matrix_fixed< double, dim, 1 >  err;
      err.set_column( 0, err_vec );      

      // normal distance
      //
      const vnl_matrix< double >  A = q->covariance_;
      vnl_symmetric_eigensystem_compute( A, V, D );

      // zeroth column is eigenvector corresponding to smallest eigen value
      // this will be the point normal
      vnl_vector< double >  evect_sm = V.get_column( 0 );

      // normal distance error projector
      vnl_matrix_fixed< double, dim, dim >  err_proj = outer_product( evect_sm, evect_sm );

      //err_proj = q->covariance_ + Jp * p->covariance_ * Jp.transpose();
      //err_proj = vnl_inverse( err_proj );

      //err_proj.fill( 0.0 );
      //err_proj( 0, 0 ) = 1.0;
      //err_proj( 1, 1 ) = 1.0;
      //err_proj( 2, 2 ) = 1.0;

      vnl_matrix_fixed< double, 1, 1 >  residual = err.transpose() * err_proj * err;

      error += residual( 0, 0 );
    }
    count += matches_[i]->to_.size();
  }

  error = vcl_sqrt( error / double( count ) );

  return error;
}


// Compute weighted alignment error using current matches.
template < unsigned int dim, unsigned int dof >
double
cdcl_estimation_ICP<dim, dof>::weighted_error()
{
  // for decomposition to get the eigenvectors from the covariances
  vnl_matrix< double >  V( dim, dim );
  vnl_vector< double >  D( dim );

  double error = 0.0;
  double sum_weight = 0.0;
  // go through all moving points
  for( typename vcl_vector< match_sptr_type >::size_type  i = 0; i < matches_.size(); ++i ) {
    // moving point p
    feature_sptr_type  p = matches_[i]->from_;   
    vnl_vector_fixed< double, dim >  mapped_moving = trans_->map_loc( p->location_ );

    // go through all matches (fixed points)
    for( typename vcl_vector< feature_sptr_type >::size_type  j = 0; j < matches_[i]->to_.size(); ++j ) {
      // fixed point q
      feature_sptr_type  q = matches_[i]->to_[j];

      // matching error
      vnl_vector_fixed< double, dim >  const &  err_vec = mapped_moving - q->location_;
      vnl_matrix_fixed< double, dim, 1 >  err;
      err.set_column( 0, err_vec );      

      // normal distance
      //
      const vnl_matrix< double >  A = q->covariance_;
      vnl_symmetric_eigensystem_compute( A, V, D );

      // zeroth column is eigenvector corresponding to smallest eigen value
      // this will be the point normal
      vnl_vector< double >  evect_sm = V.get_column( 0 );

      // normal distance error projector
      vnl_matrix_fixed< double, dim, dim >  err_proj = outer_product( evect_sm, evect_sm );

      double weight = matches_[i]->w_[j];
      vnl_matrix_fixed< double, 1, 1 >  residual = err.transpose() * err_proj * err;

      error += vcl_sqrt( residual( 0, 0 ) ) * weight;
      sum_weight += weight;
    }
  }

  error /= sum_weight;

  return error;
}


// Compute RMS error using current matches.
template < unsigned int dim, unsigned int dof >
void
cdcl_estimation_ICP<dim, dof>::compute_transfer_error_covariance()
{
  vnl_matrix_fixed< double, dim, dim >  covarJ( 0.0 );
  vnl_matrix_fixed< double, dof, dof >  covar = trans_->get_covariance();
  // go through all moving points
  for( typename vcl_vector< match_sptr_type >::size_type  i = 0; i < matches_.size(); ++i ) {
    // moving point p
    feature_sptr_type  p = matches_[i]->from_;   
    vnl_matrix_fixed< double, dim, dof > const &  Jth = trans_->jacobian_wrt_par( p->location_ );

    // transfer error covariance
    covarJ += Jth * covar * Jth.transpose();
  }

  covarJ *= ( 1.0 / double( matches_.size() ) );
  trans_->set_covarianceJ( covarJ );
}


#define PROGRESS_OUTPUT( procedure, timer )  vcl_cout << "----> Computing " << procedure << " ... " << vcl_endl; timer.mark()
#define TIMER_OUTPUT( procedure, timer )  vcl_cout << "----> Computing " << procedure << " done in " << timer.real()/1000 << "." << timer.real()%1000 << " sec." << vcl_endl


// Run estimation.
template < unsigned int dim, unsigned int dof >
bool cdcl_estimation_ICP<dim, dof>::run( void *  caller, 
                                         display_callback_type  display_points )
{
  this->initialize();

  bool converged = false;
  unsigned int max_iterations = 30;
  unsigned int iteration = 0;

  while( !converged && iteration < max_iterations ) {
    converged = this->one_iteration( caller, display_points, iteration );  
    ++iteration;
  }

  return converged;
}


// Run one EM iteration of estimation (parameters, weights, covariance, weights).
// Return true when converged.
template < unsigned int dim, unsigned int dof >
bool cdcl_estimation_ICP<dim, dof>::one_iteration( void *  caller, 
                                                   display_callback_type  display_points,
                                                   unsigned int  iteration )
{
  vul_timer timer;

  vcl_cout << "******************** Iteration " << iteration << "********************" << vcl_endl;

  PROGRESS_OUTPUT( "weights", timer );
  this->find_closest_euclidean();
  TIMER_OUTPUT( "weights", timer );
  
  if( caller ) display_points( caller, moving_[res_], fixed_[res_], matches_, trans_, iteration );
  #ifdef VTK_FOUND
  cdcl_write_matches_VTK( moving_[res_], fixed_[res_], matches_, trans_, iteration );
  #endif

  double RMS_before = this->RMS_error();

  trans_sptr_type  global_transform;
  if( trans_->is_type( cdcl_trans_rigid3d::type_id() ) ) { // if transform is rigid with small angle approximation
    // from now on, transform is incremental
    // save current transform for recomposing
    global_transform = trans_;
    trans_ = trans_->create_new(); // will create identity transformation

    trans_->set_covariance( global_transform->get_covariance() );

    // for rigid transform, we are estimating incremental transform
    // therefore transform points first
    vnl_matrix_fixed< double, dim, dim > const &  Jp = global_transform->jacobian_wrt_loc();
    // go through all moving points
    for( typename vcl_vector< match_sptr_type >::size_type  i = 0; i < matches_.size(); ++i ) {
      // moving point p
      feature_sptr_type  p = matches_[i]->from_;
      
      // reset the original from point
      // transform location and covariance
      const vnl_vector_fixed< double, dim >       location   = global_transform->map_loc( p->location_ );
      vnl_matrix_fixed< double, dim, dim > covariance = Jp * p->covariance_ * Jp.transpose();
      matches_[i]->from_ = new cdcl_feature< dim >( location, covariance );
    }
    trans_->set_parameterization( vnl_vector< double >( 6, 0.0 ) );
  }

  PROGRESS_OUTPUT( "parameters", timer );
  this->estimate_LS( true );
  TIMER_OUTPUT( "parameters", timer );

  if( trans_->is_type( cdcl_trans_rigid3d::type_id() ) ) { // if transform is rigid with small angle approximation
    // new_transform is incremental rigid transform
    global_transform->recompose_increment( trans_ );
    global_transform->set_covariance( trans_->get_covariance() );
    trans_ = global_transform;
    // from now on, transform is not incremental
  }

  trans_->print( vcl_cout );
  double RMS_after = this->RMS_error();

  ++iteration;

  //vnl_matrix_fixed< double, dim, dim >  const &  A = trans_->get_A();
  //double scale = vcl_sqrt( A(0,0)*A(0,0) + A(0,1)*A(0,1) );
  //bool converged = scale < 0.1 || scale > 10.0;
  bool converged = vcl_abs( RMS_before - RMS_after ) < 1e-4;

  return converged;
}


// Initialize estimation (compute initial weights and covariance).
template < unsigned int dim, unsigned int dof >
void cdcl_estimation_ICP<dim, dof>::initialize()
{
  trans_->print( vcl_cout );
}


template <class T>
double
rrel_util_median_abs_dev_scale( const T& begin,  const T& end, int dof = 1 )
{
  long count = long(end - begin); // VC6 & 7 has broken iterator_traits
  assert( count > 0);

  if ( count <= dof )
    return 0;

  for ( T i=begin; i!=end; ++i ) {
    *i = vcl_abs( *i );
  }
  T loc = begin + ((count-dof)/2 + dof);
  vcl_nth_element( begin, loc, end );
  return 1.4826 * (1 + 5.0/(count-dof)) * *loc;
}


// Compute weights, use Euclidean distances.
template < unsigned int dim, unsigned int dof >
void cdcl_estimation_ICP<dim, dof>::find_closest_euclidean()
{
  // query point for kd-tree
  rsdl_point  query_pt( dim, 0 );
  // use n nearest neighbors
  int n = 1;
  vcl_vector< rsdl_point >  points( n );
  vcl_vector< int >  indices( n );

  // find matches based on Euclidean distance, and store them
  //
  matches_.clear();
  // go through all moving points
  for( typename vcl_vector< feature_sptr_type >::size_type  i = 0; i < moving_[res_].size(); ++i ) {
    // moving point p
    feature_sptr_type  p = moving_[res_][i];
    vnl_vector_fixed< double, dim >  mapped_moving = trans_->map_loc( p->location_ );
   
    query_pt.set_cartesian( mapped_moving );
    // query kd-tree and find set of nearest points as potential matches
    kd_tree_[res_]->n_nearest( query_pt, n, points, indices );

    if( indices.size() > 0 ) {
      match_sptr_type  match = new cdcl_match<dim>;
      match->from_ = p;
      // go through all nearest neighbors
      for( typename vcl_vector< int >::size_type  j = 0; j < indices.size(); ++j ) {
        // fixed point q
        feature_sptr_type  q = fixed_[res_][indices[j]];
    
        match->to_.push_back( q );
        // temporarily set weight to 1.0
        match->w_.push_back( 1.0 );
      }
      
      matches_.push_back( match );
    }
  }

}

#define NORMAL_DISTANCES_FOR_ICP 1

// Estimate scale from current matches and assign robust weight.
template < unsigned int dim, unsigned int dof >
double cdcl_estimation_ICP<dim, dof>::estimate_scale_and_assign_weight()
{
  // for decomposition to get the eigenvectors from the covariances
  vnl_matrix< double >  V( dim, dim );
  vnl_vector< double >  D( dim );
  vcl_vector< double >  distances;
  vnl_matrix_fixed< double, dim, dim >  Jp = trans_->jacobian_wrt_loc();
  // compute residuals and store them to compute scale later
  //
  //matches_.clear();
  // go through all moving points
  double sum_weighted_error = 0.0;
  double sum_weights = 0.0;
  for( typename vcl_vector< match_sptr_type >::size_type  i = 0; i < matches_.size(); ++i ) {
    feature_sptr_type  p = matches_[i]->from_;

    vnl_vector_fixed< double, dim >  mapped_moving = trans_->map_loc( p->location_ );

    // go through all matches
    for( typename vcl_vector< match_sptr_type >::size_type  j = 0; j < matches_[i]->w_.size(); ++j ) {
      feature_sptr_type  q = matches_[i]->to_[j];
      // matching error
      vnl_vector_fixed< double, dim >  const &  err_vec = mapped_moving - q->location_;
      vnl_matrix_fixed< double, dim, 1 >  err;
      err.set_column( 0, err_vec );

#if NORMAL_DISTANCES_FOR_ICP
      // normal distance
      //
      const vnl_matrix< double >  A = q->covariance_;
      vnl_symmetric_eigensystem_compute( A, V, D );

      // zeroth column is eigenvector corresponding to smallest eigen value
      // this will be the point normal
      vnl_vector< double >  evect_sm = V.get_column( 0 );

      // normal distance error projector
      vnl_matrix_fixed< double, dim, dim >  err_proj = outer_product( evect_sm, evect_sm );

      vnl_matrix_fixed< double, 1, 1 >  residual = err.transpose() * err_proj * err;
      //vcl_cout << j << ": " << err_vec.two_norm() << "  " << residual( 0, 0 ) << "  " << matches_[i]->w_[j] << vcl_endl;
#else
      // Euclidean distance
      double residual = vnl_tranpose( err ) * err;
#endif

      if( scale_estimated_ ) {
        sum_weighted_error += matches_[i]->w_[j] * residual( 0, 0 );
        sum_weights += matches_[i]->w_[j];
      }
       
      // temporarily store squared residual
      matches_[i]->w_[j] = residual( 0, 0 );
      distances.push_back( vcl_sqrt( residual( 0, 0 ) ) );
    }

  }


  double scale = 0.0;
  if( !scale_estimated_ ) {
    // assuming that distance is Gaussian distributed
    // compute maximum absolute deviation to estimate scale
    // scale is std. dev. of error residuals
    //
    scale = rrel_util_median_abs_dev_scale( distances.begin(), distances.end() );
vcl_cout << "SCALE FROM MAD: " << scale;
    
    // muse and unwgted_scale_est are used in the first iteration.
    // do NOT use sk refinement
    const bool max_lookup_table_size = 0;   //  Prevents costly and unnecessary preliminary generation of a table
    const bool use_sk_refinement_in_muse = false;
    vcl_auto_ptr< rrel_objective >  obj( new rrel_muset_obj( max_lookup_table_size,
                         use_sk_refinement_in_muse ) );
    scale = obj->scale( distances.begin(), distances.end() );
vcl_cout << "  SCALE FROM MUSE: " << scale << std::endl;
    scale_estimated_ = true;
  }
  else {
    //const double epsilon = 1e-16;
    scale = vcl_sqrt( sum_weighted_error / sum_weights );
  }

  double constant = 2.0;
  //double constant = 1.0;
  scale = constant * scale;
  vcl_cout << "SCALE: " << scale << " residuals size: " << distances.size() << vcl_endl;

  // now that scale has been estimated, compute robust weight
  // keep only matches that have non-zero robust weight
  vcl_vector< match_sptr_type >  matches = matches_;
  matches_.clear();
  //matches_ = matches;
  // go through all moving points
  double sum_grad_w = 0.0;
  for( typename vcl_vector< match_sptr_type >::size_type  i = 0; i < matches.size(); ++i ) {
    feature_sptr_type  p = matches[i]->from_;
    
    match_sptr_type  match = new cdcl_match<dim>;
    match->from_ = p;

    // go through all matches
    for( typename vcl_vector< match_sptr_type >::size_type  j = 0; j < matches[i]->to_.size(); ++j ) {
      feature_sptr_type  q = matches[i]->to_[j];

      // retrieve temporarily stored squared residual from matches[i]->w_[j] and compute robust weight
      double scale2 = scale*scale;
      double w = weight_BT( matches[i]->w_[j] / scale2 ) / scale2;  // scale*scale because residual has been squared

      if( this->weight_by_strength_ ) {
        double grad_w = p->strength_ * q->strength_;
        w *= grad_w;
        sum_grad_w += grad_w;
      }

// Use all matches for ICP because:
// Matches are normalized, then this->estimate_scale_and_assign_weight(); is called where only matches with nonzero weights are kept.
// But normalization normalizes them before, then match sets are changed and unnormalization DOES NOT unnormalize all of them (because they got dropped -- they had zero weight).
//      if( w > 0.0 ) {
        match->to_.push_back( q );
        match->w_.push_back( w );
//      }      
    }

    if( match->to_.size() > 0 ) {
      matches_.push_back( match );
    }

  }

  // normalize the gradient weight
  if( this->weight_by_strength_ ) {
    for( typename vcl_vector< match_sptr_type >::size_type  i = 0; i < matches_.size(); ++i ) {

      // go through all matches
      for( typename vcl_vector< match_sptr_type >::size_type  j = 0; j < matches_[i]->to_.size(); ++j ) {
        matches_[i]->w_[j] *= ( 1.0 / sum_grad_w );
      }

    }
  }


  return scale;
}


// Estimate initial covariance matrix.
template < unsigned int dim, unsigned int dof >
void cdcl_estimation_ICP<dim, dof>::estimate_LS( bool estimate_parameters )
{
  // Normalize before estimation.
  //
  // Normalization parameters.
  vnl_vector_fixed< double, dim >  center_moving( 0.0 ), center_fixed( 0.0 );
  double avg_rad_moving( 0.0 ), avg_rad_fixed( 0.0 );

  if( normalize_matches_ ) {
    this->cdcl_normalize_matches( center_moving, avg_rad_moving, center_fixed, avg_rad_fixed );
    this->trans_->normalize( avg_rad_moving, avg_rad_fixed, center_moving, center_fixed );
  }

  double scale = this->estimate_scale_and_assign_weight();
  if( scale < 0.005 ) scale = 0.005; // in normalized coordinates

  // Estimate initial matrix by setting up a least squares problem
  // A x = q
  // A^T A
  vnl_matrix_fixed< double, dof, dof >  AtA( 0.0 );
  vnl_matrix_fixed< double, dof, 1 >    Atb( 0.0 );
  
  vnl_matrix< double >  V( dim, dim );
  vnl_vector< double >  D( dim );
  vnl_matrix_fixed< double, dim, dim >  Jp = trans_->jacobian_wrt_loc();
  vnl_matrix_fixed< double, dim, 1 >  q_loc;
  double sum_w = 0.0;
  // go through all moving points
  for( typename vcl_vector< feature_sptr_type >::size_type  i = 0; i < matches_.size(); ++i ) {
    // moving point p
    feature_sptr_type  p = matches_[i]->from_;
    vnl_vector_fixed< double, dim >  mapped_moving = trans_->map_loc( p->location_ );

    // Jacobian w.r.t. parameters
    vnl_matrix_fixed< double, dim, dof > const &  Jth = trans_->jacobian_wrt_par( p->location_ );
    //vcl_cout << p->location_ << "  pc  " << trans_->center_moving_ << vcl_endl;

    for( typename vcl_vector< feature_sptr_type >::size_type  j = 0; j < matches_[i]->to_.size(); ++j ) {
      // fixed point q
      feature_sptr_type  q = matches_[i]->to_[j];
      if( trans_->is_type( cdcl_trans_rigid3d::type_id() ) ) { // if transform is rigid with small angle approximation
        // rigid is: p' = p + R p + t
        q_loc.set_column( 0, q->location_-p->location_ );
      }
      else {
        q_loc.set_column( 0, q->location_ );
      }

      double weight = matches_[i]->w_[j];

#if NORMAL_DISTANCES_FOR_ICP
      // normal distance
      
      const vnl_matrix< double >  A = q->covariance_;
      vnl_symmetric_eigensystem_compute( A, V, D );

      // zeroth column is eigenvector corresponding to smallest eigen value
      // this will be the point normal
      vnl_vector< double >  evect_sm = V.get_column( 0 );

      // normal distance error projector
      vnl_matrix_fixed< double, dim, dim >  err_proj = outer_product( evect_sm, evect_sm );
      
      AtA += Jth.transpose() * err_proj * Jth * weight;
      Atb += Jth.transpose() * err_proj * q_loc * weight;
#else
      // Euclidean distance
      AtA += weight * Jth.transpose() * Jth;
      Atb += weight * Jth.transpose() * q_loc;
#endif

      sum_w += weight;

      //vcl_cout << Jth << vcl_endl << AtA << vcl_endl << vcl_endl;

    }
  }
  //AtA *= ( 1.0 / sum_w );
  //Atb *= ( 1.0 / sum_w );

  vnl_svd< double > svd( AtA );
  vnl_matrix_fixed< double, dof, dof >  invAA = svd.inverse();

  // don't need to multiply by scale, since it is already included in the robust weights [yang:cvpr04]
  trans_->set_covariance( /*scale*scale **/ /** 2.0*/ /**/ invAA );
//removing scale*scale here and putting it to weight computation does not produce the same result - why?, check weight based on strength normalization
  if( estimate_parameters ) {
    vnl_matrix_fixed< double, dof, 1 >  new_params = invAA * Atb;
    //vcl_cout << "PARAMS: " << new_params << vcl_endl;
    trans_->set_parameterization( new_params.get_column( 0 ) );
  }

  if( normalize_matches_ ) {
    // convert transformation back to unnormalized coordinates
    this->trans_->unnormalize( avg_rad_moving, avg_rad_fixed, center_moving, center_fixed );
    this->cdcl_unnormalize_matches( center_moving, avg_rad_moving, center_fixed, avg_rad_fixed );
  }

  vcl_cout << "Initial covariance: " << vcl_endl << trans_->get_covariance() << vcl_endl << "Scale: " << scale << vcl_endl;

}


// Normalize moving and fixed points which the matches were formed with
template < unsigned int dim, unsigned int dof >
void cdcl_estimation_ICP<dim, dof>::cdcl_normalize_matches(
                             vnl_vector_fixed< typename feature_type::coord_type, dim >                   & center_moving, 
                             typename feature_type::coord_type                                            & avg_radius_moving, 
                             vnl_vector_fixed< typename feature_type::coord_type, dim >                   & center_fixed, 
                             typename feature_type::coord_type                                            & avg_radius_fixed )
{
  center_moving.fill( 0.0 );
  avg_radius_moving = 0.0;
  center_fixed.fill( 0.0 );
  avg_radius_fixed = 0.0;

  // set of fixed points that got matched to moving points
  // utilize unique container concept, i.e. no two points are identical
  // we need to make sure that points don't get normalized twice
  vcl_set< feature_sptr_type >  fixed_matched;

  // compute data center of mass
  //
  // go through all moving points
  for( typename vcl_vector< match_sptr_type >::size_type  i = 0; i < this->matches_.size(); ++i ) {
    // moving point p
    feature_sptr_type  p = this->matches_[i]->from_;
    center_moving += p->location_;

    // go through all matches (fixed points)
    for( typename vcl_vector< match_sptr_type >::size_type  j = 0; j < this->matches_[i]->to_.size(); ++j ) {
      // fixed point q
      feature_sptr_type  q = this->matches_[i]->to_[j];
      // if it can be inserted (i.e. not already there), add to center
      if( fixed_matched.insert( q ).second ) center_fixed += q->location_;
    }
  }
  center_moving /= this->matches_.size();
  center_fixed /= fixed_matched.size();


  // center the data and compute average radius
  //
  vnl_vector_fixed< typename feature_type::coord_type, dim >  centered;
  for( typename vcl_vector< match_sptr_type >::size_type  i = 0; i < this->matches_.size(); ++i ) {
    // moving point p
    feature_sptr_type  p = this->matches_[i]->from_;
    p->location_ -= center_moving;
    avg_radius_moving += p->location_.magnitude();
  }
  avg_radius_moving /= this->matches_.size();

  // go through all matches (fixed points)
  for( typename vcl_set< feature_sptr_type >::iterator  it = fixed_matched.begin(); it != fixed_matched.end(); ++it ) {
    // fixed point q
    feature_sptr_type  q = *it;
    q->location_ -= center_fixed;
    avg_radius_fixed += q->location_.magnitude();
  }
  avg_radius_fixed /= fixed_matched.size();;


  // assign normalized data
  //
  for( typename vcl_vector< match_sptr_type >::size_type  i = 0; i < this->matches_.size(); ++i ) {
    // moving point p
    feature_sptr_type  p = this->matches_[i]->from_;

    // convert locations to normalized coordinate system
    p->location_ /= avg_radius_moving;

    // convert covariances to normalized coordinate system
    p->covariance_ /= (avg_radius_moving*avg_radius_moving);
  }

  // go through all matches (fixed points)
  for( typename vcl_set< feature_sptr_type >::iterator  it = fixed_matched.begin(); it != fixed_matched.end(); ++it ) {
    // fixed point q
    feature_sptr_type  q = *it;
    // convert locations to normalized coordinate system
    q->location_ /= avg_radius_fixed;

    // convert covariances to normalized coordinate system
    q->covariance_ /= (avg_radius_fixed*avg_radius_fixed);
  }

}


// Unormalize moving and fixed points which the matches were formed with
template < unsigned int dim, unsigned int dof >
void cdcl_estimation_ICP<dim, dof>::cdcl_unnormalize_matches(
                               vnl_vector_fixed< typename feature_type::coord_type, dim >                   const & center_moving, 
                               typename feature_type::coord_type                                            const & avg_radius_moving, 
                               vnl_vector_fixed< typename feature_type::coord_type, dim >                   const & center_fixed, 
                               typename feature_type::coord_type                                            const & avg_radius_fixed )
{
  // set of fixed points that got matched to moving points
  // utilize unique container concept, i.e. no two points are identical
  // we need to make sure that points don't get unnormalized twice
  vcl_set< feature_sptr_type >  fixed_matched;

  // remove radius normalization and centering of the data
  //
  for( typename vcl_vector< match_sptr_type >::size_type  i = 0; i < this->matches_.size(); ++i ) {
    // moving point p
    feature_sptr_type  p = this->matches_[i]->from_;

    // convert locations to unnormalized coordinate system
    p->location_ *= avg_radius_moving;
    p->location_ += center_moving;

    // convert covariances to unnormalized coordinate system
    p->covariance_ *= (avg_radius_moving*avg_radius_moving);

    // go through all matches (fixed points)
    for( typename vcl_vector< match_sptr_type >::size_type  j = 0; j < this->matches_[i]->to_.size(); ++j ) {
      // fixed point q
      feature_sptr_type  q = this->matches_[i]->to_[j];

      // if it was inserted (i.e. not already there), unnormalize
      if( fixed_matched.insert( q ).second ) {
        // convert locations to unnormalized coordinate system
        q->location_ *= avg_radius_fixed;
        q->location_ += center_fixed;

        // convert covariances to unnormalized coordinate system
        q->covariance_ *= (avg_radius_fixed*avg_radius_fixed);
      }
    }
  }

}


#define CDCL_ESTIMATION_ICP_INSTANTIATE( dim, dof ) \
template                                            \
class cdcl_estimation_ICP< dim, dof >;
