#include <vcl_utility.h>
#include <vcl_algorithm.h>
#include <vcl_iomanip.h>
#include <vcl_iostream.h>
#include <vcl_iterator.h>
#include <vcl_cstdlib.h>
#include <vcl_set.h>

#include <vul/vul_timer.h>
#include <vnl/algo/vnl_svd.h>
#include <vnl/vnl_inverse.h>
#include <vnl/vnl_transpose.h>
#include <vnl/vnl_random.h>
#include <vnl/vnl_cross.h>

#include <vnl/algo/vnl_symmetric_eigensystem.h>
#include <vnl/algo/vnl_generalized_eigensystem.h>

#include <rrel/rrel_muset_obj.h>
#include <rsdl/rsdl_dist.h>
#include <vcl_memory.h>

#include "rrl_estimation_ICP_matching_all.h"
#include <cdcl/cdcl_utils.h>
#include <cdcl/cdcl_keypoint.h>

#include "itkBoundingBox.h"
#include "itkDanielssonDistanceMapImageFilter.h"

#undef VTK_FOUND
// VTK is used to dump data for analysis in external programs.
#ifdef VTK_FOUND
#include <cdcl/cdcl_utils_VTK.h>
#endif


//#define DEBUG( statement ) statement
#define DEBUG( statement ) 0


//:
// \file
// \brief  Estimation handling the whole registration: 1) matching,
//         2) parameter estimate, 3) covariance estimate.
//         ICP with Euclidean distances for matching and normal distances
//         for estimation is used. Features cdcl_feature_ICP are used.
//         This is the same as cdcl_estimation_ICP_matching_all, but with
//         different feature type and taking out all unncecessary computation.
// \author Michal Sofka
// \date   Feb 2008


#define USE_KD_TREES 0


// Construct estimation object given moving and fixed feature sets
// and an initial transformation.
template < unsigned int dim, unsigned int dof >
rrl_estimation_ICP_matching_all<dim, dof>::rrl_estimation_ICP_matching_all( vcl_vector< vcl_vector< feature_sptr_type > > const &  moving, 
                                                                          vcl_vector< vcl_vector< feature_sptr_type > > const &  fixed,
                                                                          trans_sptr_type                               const &  trans,
                                                                          vcl_vector< double >                          const &  fixed_spacing )
  : moving_( moving ),
    fixed_( fixed ),
    fixed_spacing_( fixed_spacing ),
    trans_( trans ),
    normalize_matches_( false ),
    weight_by_strength_( false ),
    scale_estimated_( false ),
    fixed_voronoi_map_( 0 ),
    number_points_( 1000 )
{
  // set the coarsest resolution
  res_ = max_res_ = vcl_min( moving_.size(), fixed_[feature_type::CORNER].size() ) - 1;

  #if USE_KD_TREES
  this->build_kd_trees();
  #endif

  if( number_points_ > moving_[res_].size() ) {
    number_points_ = moving_[res_].size();
  }
}


// Construct estimation object given moving and fixed feature sets
// and an initial transformation.
template < unsigned int dim, unsigned int dof >
rrl_estimation_ICP_matching_all<dim, dof>::rrl_estimation_ICP_matching_all( vcl_vector< feature_sptr_type > const &  moving, 
                                                    vcl_vector< feature_sptr_type > const &  fixed,
                                                    trans_sptr_type                 const &  trans )
  : trans_( trans ),
    normalize_matches_( false ),
    weight_by_strength_( false ),
    scale_estimated_( false ),
    fixed_voronoi_map_( 0 ),
    number_points_( 1000 )
{
  moving_.push_back( moving );
  fixed_.push_back( fixed );
  fixed_spacing_.push_back( 1.0 );  // this is for switching resolutions and gets set in subsample

  // set the coarsest resolution
  res_ = max_res_ = vcl_min( moving_.size(), fixed_[feature_type::CORNER].size() ) - 1;

  #if USE_KD_TREES
  this->build_kd_trees();
  #endif

  if( number_points_ > moving_[res_].size() ) {
    number_points_ = moving_[res_].size();
  }
}


template < unsigned int dim, unsigned int dof >
bool rrl_estimation_ICP_matching_all<dim, dof>::set_number_points( unsigned int number_points )
{
  vcl_cout << "number_points: " << number_points << vcl_endl;
  if( number_points >= moving_[res_].size() ) {
    number_points_ = moving_[res_].size();
    vcl_cout << "number_points_: " << number_points_ << vcl_endl;
    return true;
  }
  else {
    number_points_ = number_points;
    vcl_cout << "number_points_: " << number_points_ << vcl_endl;
    return false;
  }
}


// Build kd trees, one for each resolution level of the fixed image.
template < unsigned int dim, unsigned int dof >
void rrl_estimation_ICP_matching_all<dim, dof>::build_kd_trees()
{
  //save_data( moving_, "moving" );
  //save_data( fixed_, "fixed" );

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


// Set initial transform.
template < unsigned int dim, unsigned int dof >
void
rrl_estimation_ICP_matching_all<dim, dof>::
set_transform( trans_sptr_type const &  trans )
{
  trans_ = trans;

  // new transform has been set -> will need to reestimate scale
  scale_estimated_ = false;

  number_points_ = 1000;
  if( number_points_ > moving_[res_].size() ) {
    number_points_ = moving_[res_].size();
  }
}


// Compute median error using current matches.
template < unsigned int dim, unsigned int dof >
double
rrl_estimation_ICP_matching_all<dim, dof>::median_error()
{
  vcl_vector< double >  errors;
  // go through all moving points
  for( typename vcl_vector< match_sptr_type >::size_type  i = 0; i < matches_.size(); ++i ) {
    // moving point p
    feature_sptr_type  p = matches_[i]->from_;   
    vnl_vector_fixed< double, dim >  mapped_moving = trans_->map_loc( p->location_ );

    // fixed point q
    feature_sptr_type  q = matches_[i]->to_[0];

    // matching error
    vnl_vector_fixed< double, dim >  const &  err_vec = mapped_moving - q->location_;
    vnl_matrix_fixed< double, dim, dim >  err_proj = q->error_projector_;

    double residual2 = dot_product( err_vec.as_vector().post_multiply( err_proj ), err_vec );
    if( residual2 < 0.0 ) residual2 = 0.0; // it can happen that the residual is very very small negative number

    errors.push_back( vcl_sqrt( residual2 ) );

  }

  vcl_vector< double >::iterator  loc = errors.begin() + (errors.size()-1)/2 + 1;
  vcl_nth_element( errors.begin(), loc, errors.end() );
  double med_error = *loc;

  return med_error;
}


// Compute RMS error using current matches.
template < unsigned int dim, unsigned int dof >
double
rrl_estimation_ICP_matching_all<dim, dof>::RMS_error()
{
  vnl_matrix_fixed< double, dim, dim >  Jp = trans_->jacobian_wrt_loc();

  double error = 0.0;
  unsigned int count = 0;
  // go through all moving points
  for( typename vcl_vector< match_sptr_type >::size_type  i = 0; i < matches_.size(); ++i ) {
    // moving point p
    feature_sptr_type  p = matches_[i]->from_;   
    vnl_vector_fixed< double, dim >  mapped_moving = trans_->map_loc( p->location_ );

    // fixed point q
    feature_sptr_type  q = matches_[i]->to_[0];

    // matching error
    vnl_vector_fixed< double, dim > const &  err_vec = mapped_moving - q->location_;

    vnl_matrix_fixed< double, dim, dim > const &  err_proj = q->error_projector_;

    // e^T  error_projector  e
    double residual2 = dot_product( err_vec.as_vector().post_multiply( err_proj ), err_vec );
    if( residual2 < 0.0 ) residual2 = 0.0; // it can happen that the residual is very very small negative number

    error += residual2;
  }

  if( count != 0 ) error = vcl_sqrt( error / double( matches_.size() ) );

  return error;
}


// Compute weighted alignment error using current matches.
template < unsigned int dim, unsigned int dof >
double
rrl_estimation_ICP_matching_all<dim, dof>::weighted_error()
{
  double error = 0.0;
  double sum_weight = 0.0;
  // go through all moving points
  for( typename vcl_vector< match_sptr_type >::size_type  i = 0; i < matches_.size(); ++i ) {
    // moving point p
    feature_sptr_type  p = matches_[i]->from_;   
    vnl_vector_fixed< double, dim >  mapped_moving = trans_->map_loc( p->location_ );

    // fixed point q
    feature_sptr_type  q = matches_[i]->to_[0];

    // matching error
    vnl_vector_fixed< double, dim >  const &  err_vec = mapped_moving - q->location_;
    vnl_matrix_fixed< double, dim, dim >  err_proj = q->error_projector_;

    double residual2 = dot_product( err_vec.as_vector().post_multiply( err_proj ), err_vec );
    if( residual2 < 0.0 ) residual2 = 0.0; // it can happen that the residual is very very small negative number

    double weight = matches_[i]->w_[0];

    error += vcl_sqrt( residual2 ) * weight;
    sum_weight += weight;
  }

  if( sum_weight != 0 ) error /= sum_weight;

  return error;
}


template < unsigned int dim, unsigned int dof >
double
rrl_estimation_ICP_matching_all<dim, dof>::sheet_angles()
{
  // for decomposition to get the eigenvectors from the covariances
  vnl_matrix< double >  V( dim, dim );
  vnl_vector< double >  D( dim );

  double avg_angle_difference = 0.0;
  double sum_weight = 0.0;
  // go through all moving points
  for( typename vcl_vector< match_sptr_type >::size_type  i = 0; i < matches_.size(); ++i ) {
    // moving point p
    feature_sptr_type  p = matches_[i]->from_;   

    if( p->shape_ != feature_type::SHEET ) continue;

    // fixed point q
    feature_sptr_type  q = matches_[i]->to_[0];

    if( q->shape_ != feature_type::SHEET ) continue;

    vnl_vector_fixed< double, dim >  mapped_moving = trans_->map_loc( p->location_ );

    // map normal direction
    vnl_vector_fixed< double, dim >  p_loc_nor = p->location_ + p->directions_[0];
    vnl_vector_fixed< double, dim >  p_loc_nor_mapped = this->trans_->map_loc( p_loc_nor );

    vnl_vector_fixed< double, dim >  p_dir_mapped = ( p_loc_nor_mapped - mapped_moving ).normalize();

    double angle_difference = vcl_acos( dot_product( p_dir_mapped, q->directions_[0] ) );

    if( angle_difference > vnl_math::pi_over_2 )  angle_difference = vnl_math::pi - angle_difference;

    double weight = matches_[i]->w_[0];

    avg_angle_difference += weight * angle_difference;
    sum_weight += weight;

  }

  avg_angle_difference /= sum_weight;

  return avg_angle_difference;
}


template < unsigned int dim, unsigned int dof >
double
rrl_estimation_ICP_matching_all<dim, dof>::tube_angles()
{
  // for decomposition to get the eigenvectors from the covariances
  vnl_matrix< double >  V( dim, dim );
  vnl_vector< double >  D( dim );

  double avg_angle_difference = 0.0;
  double sum_weight = 0.0;
  // go through all moving points
  for( typename vcl_vector< match_sptr_type >::size_type  i = 0; i < matches_.size(); ++i ) {
    // moving point p
    feature_sptr_type  p = matches_[i]->from_;   

    if( p->shape_ != feature_type::TUBE ) continue;

    // fixed point q
    feature_sptr_type  q = matches_[i]->to_[0];

    if( q->shape_ != feature_type::TUBE ) continue;

    vnl_vector_fixed< double, dim >  mapped_moving = trans_->map_loc( p->location_ );

    // compute tangent directions
    vnl_vector_fixed< coord_type, dim >  p_tangent = vnl_cross_3d( p->directions_[0], p->directions_[1] );
    vnl_vector_fixed< coord_type, dim >  q_tangent = vnl_cross_3d( q->directions_[0], q->directions_[1] );

    // map tangent direction
    vnl_vector_fixed< double, dim >  p_loc_tan = p->location_ + p_tangent;
    vnl_vector_fixed< double, dim >  p_loc_tan_mapped = this->trans_->map_loc( p_loc_tan );

    vnl_vector_fixed< double, dim >  p_dir_mapped = ( p_loc_tan_mapped - mapped_moving ).normalize();

    //double angle_difference = vcl_acos( dot_product( p_tangent, q_tangent ) );
    double angle_difference = vcl_acos( dot_product( p_dir_mapped, q_tangent ) );

    if( angle_difference > vnl_math::pi_over_2 )  angle_difference = vnl_math::pi - angle_difference;

    double weight = matches_[i]->w_[0];

    avg_angle_difference += weight * angle_difference;
    sum_weight += weight;

  }

  avg_angle_difference /= sum_weight;

  return avg_angle_difference;
}


// Compute RMS error using current matches.
template < unsigned int dim, unsigned int dof >
void
rrl_estimation_ICP_matching_all<dim, dof>::compute_transfer_error_covariance()
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
bool rrl_estimation_ICP_matching_all<dim, dof>::run( void *  caller, 
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
bool rrl_estimation_ICP_matching_all<dim, dof>::one_iteration( void *  caller, 
                                                               display_callback_type  display_points,
                                                               unsigned int  iteration )
{
  vul_timer timer;

  vcl_cout << "******************** Iteration " << iteration << "********************" << vcl_endl;

  PROGRESS_OUTPUT( "weights", timer );
  this->find_closest_euclidean();
  TIMER_OUTPUT( "weights", timer );
  
  // WARNING: only sheets are passed into the callback (this does not matter for some callbacks which only use matches_ to display points)
  if( caller ) display_points( caller, moving_[res_], fixed_[res_], matches_, trans_, iteration );
  #ifdef VTK_FOUND
  cdcl_write_matches_VTK( moving_[res_], fixed_[res_], matches_, trans_, iteration );
  #endif

  double RMS_before = this->RMS_error();

  PROGRESS_OUTPUT( "parameters", timer );
  this->estimate_LS( true );
  TIMER_OUTPUT( "parameters", timer );

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
void rrl_estimation_ICP_matching_all<dim, dof>::initialize()
{
  DEBUG( trans_->print( vcl_cout ) );
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
void rrl_estimation_ICP_matching_all<dim, dof>::find_closest_euclidean()
{
  VoronoiMapImageType::SizeType  voronoiSize = fixed_voronoi_map_->GetLargestPossibleRegion().GetSize();

  if( !fixed_voronoi_map_ ) vcl_cout << "Error: Fixed Voronoi Map must be set!" << vcl_endl;

  // query point for kd-tree
  rsdl_point  query_pt( dim, 0 );
  // use n nearest neighbors
  int n = 1;
  vcl_vector< rsdl_point >  points( n );
  vcl_vector< int >  indices( n );

  // we will pick number_points_ number of moving points randomly
  vcl_random_shuffle( moving_[res_].begin(), moving_[res_].end() );

  // find matches based on Euclidean distance, and store them
  //
  matches_.clear();
  vcl_cout << "Using: " << number_points_ << " moving points." << vcl_endl;
  // go through all moving points
//  for( typename vcl_vector< feature_sptr_type >::size_type  i = 0; i < moving_[res_].size(); ++i ) {
  for( typename vcl_vector< feature_sptr_type >::size_type  i = 0; i < number_points_; ++i ) {
    // moving point p
    feature_sptr_type  p = moving_[res_][i];
    vnl_vector_fixed< double, dim >  mapped_moving = trans_->map_loc( p->location_ );
   
    #if USE_KD_TREE
    
    query_pt.set_cartesian( mapped_moving );
    // query kd-tree and find set of nearest points as potential matches
    kd_tree_[res_]->n_nearest( query_pt, n, points, indices );
    feature_sptr_type  q = fixed_[res_][indices[0]];
    
    #else

      itk::Point< double, 3 >  mapped_moving_itk;
      mapped_moving_itk[0] = mapped_moving[0];
      mapped_moving_itk[1] = mapped_moving[1];
      mapped_moving_itk[2] = mapped_moving[2];
      VoronoiMapImageType::IndexType  mapped_moving_index;
      fixed_voronoi_map_->TransformPhysicalPointToIndex( mapped_moving_itk, mapped_moving_index );

      for( unsigned int d = 0; d < 3; ++d ) {
        if( mapped_moving_index[d] < 0 ) mapped_moving_index[d] = 0;
        if( mapped_moving_index[d] >= VoronoiMapImageType::IndexValueType( voronoiSize[d] ) ) mapped_moving_index[d] = VoronoiMapImageType::IndexValueType( voronoiSize[d] - 1 );
      }
      unsigned int fixed_ind = fixed_voronoi_map_->GetPixel( mapped_moving_index );

      feature_sptr_type  q = fixed_[res_][fixed_ind];
    
    #endif

    // map normal direction
    vnl_vector_fixed< double, dim >  p_loc_nor = p->location_ + p->directions_[0];
    vnl_vector_fixed< double, dim >  p_loc_nor_mapped = this->trans_->map_loc( p_loc_nor );

    vnl_vector_fixed< double, dim >  p_dir_mapped = ( p_loc_nor_mapped - mapped_moving ).normalize();


      //vcl_cout << fixed_ind << "  " << q << "  " << p << vcl_endl;
  //    if( p->shape_ == feature_type::SHEET && q->shape_ == feature_type::SHEET ) if( vcl_acos( dot_product( p_dir_mapped, q->directions_[0] ) ) > 45.0 * vnl_math::pi / 180.0 ) { vcl_cout << "skipping" << vcl_endl; continue; }
    //if( p->shape_ == feature_type::CORNER || q->shape_ == feature_type::CORNER ) continue;
      match_sptr_type  match = new match_type;
      //feature_sptr_type  p_clone = dynamic_cast< feature_type* >( p->clone().as_pointer() );
      match->from_ = p;

      //match->to_.push_back( dynamic_cast< feature_type* >( q->clone().as_pointer() ) );
      match->to_.push_back( q );
      // temporarily set weight to 1.0
      match->w_.push_back( 1.0 );
      
      matches_.push_back( match );


  }

}

#define NORMAL_DISTANCES_FOR_ICP 1

// Estimate scale from current matches and assign robust weight.
template < unsigned int dim, unsigned int dof >
double rrl_estimation_ICP_matching_all<dim, dof>::estimate_scale_and_assign_weight()
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

  //vcl_cout << "ESTIMATING SCALE: " << vcl_endl;
  //trans_->print( vcl_cout );
  for( typename vcl_vector< match_sptr_type >::size_type  i = 0; i < matches_.size(); ++i ) {
    feature_sptr_type  p = matches_[i]->from_;

    vnl_vector_fixed< double, dim >  mapped_moving = trans_->map_loc( p->location_ );

    feature_sptr_type  q = matches_[i]->to_[0];

    // matching error
    vnl_vector_fixed< double, dim >  const &  err_vec = mapped_moving - q->location_;

    vnl_matrix_fixed< double, dim, dim >  err_proj = q->error_projector_;
    double residual2 = dot_product( err_vec.as_vector().post_multiply( err_proj ), err_vec );
    if( residual2 < 0.0 ) residual2 = 0.0; // it can happen that the residual is very very small negative number

    //vcl_cout << mapped_moving << "  " << q->location_ << "  " << residual2 << vcl_endl;
    if( scale_estimated_ ) {
      sum_weighted_error += matches_[i]->w_[0] * residual2;
      sum_weights += matches_[i]->w_[0];
    }
     
    // temporarily store squared residual
    matches_[i]->w_[0] = residual2;
    distances.push_back( vcl_sqrt( residual2 ) );
  }

  double scale = 0.0;
  if( !scale_estimated_ ) {
    // assuming that distance is Gaussian distributed
    // compute maximum absolute deviation to estimate scale
    // scale is std. dev. of error residuals
    //
    //scale = rrel_util_median_abs_dev_scale( distances.begin(), distances.end() );
    //vcl_cout << "SCALE FROM MAD: " << scale;
    
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

  //double constant = 2.0;
  double constant = 1.0;
  double sum_grad_w = 0.0;

  scale = constant * scale;
  vcl_cout << "SCALE: " << scale << " residuals size: " << distances.size() << vcl_endl;

  // now that scale has been estimated, compute robust weight
  // go through all moving points
  for( typename vcl_vector< match_sptr_type >::size_type  i = 0; i < matches_.size(); ++i ) {
    feature_sptr_type  p = matches_[i]->from_;
    
    feature_sptr_type  q = matches_[i]->to_[0];

    // retrieve temporarily stored squared residual from matches[i]->w_[j] and compute robust weight
    double scale2 = scale*scale;
    double w = weight_BT( matches_[i]->w_[0] / scale2 ) / scale2;  // 1st scale*scale because residual has been squared, 2nd so that we don't have to do it later
    //double w = weight_Cauchy( matches_[i]->w_[0] / scale2 ) / scale2;  // 1st scale*scale because residual has been squared, 2nd so that we don't have to do it later

    matches_[i]->w_[0] = w;

    //if( weight_by_strength_ && w != 0 ) { // the later condition can be removed once we don't store points with zero weights
    //  double grad_w = p->strength_ * q->strength_;
    //  matches_[i]->w_[0] *= grad_w;
    //  sum_grad_w += grad_w;
    //}


  }



// FIX THIS FOR SPEEDUP -- RIGHT NOW KEEPING ALL MATCHES (EVEN THE ONES THAT ARE ZERO) -- read note below
// HOW -- just copy over all non-zero (but after all thre scales and weights have been computed; pay attention to the note below
//
//  //double constant = 2.0;
//  double constant = 1.0;
//  scale = constant * scale;
//  vcl_cout << "SCALE: " << scale << " residuals size: " << distances.size() << vcl_endl;
//
//  // now that scale has been estimated, compute robust weight
//  // keep only matches that have non-zero robust weight
//  vcl_vector< match_sptr_type >  matches = matches_;
//  matches_.clear();
//  // go through all moving points
//  double sum_grad_w = 0.0;
//  for( typename vcl_vector< match_sptr_type >::size_type  i = 0; i < matches.size(); ++i ) {
//    feature_sptr_type  p = matches[i]->from_;
//    
//    match_sptr_type  match = new match_type;
//    match->from_ = p;
//
//    if( q->shape_ != shape ) continue;
//
//    // go through all matches
//    for( typename vcl_vector< match_sptr_type >::size_type  j = 0; j < matches[i]->to_.size(); ++j ) {
//      feature_sptr_type  q = matches[i]->to_[j];
//
//      // retrieve temporarily stored squared residual from matches[i]->w_[j] and compute robust weight
//      double scale2 = scale*scale;
//      double w = weight_BT( matches[i]->w_[j] / scale2 ) / scale2;  // 1st scale*scale because residual has been squared, 2nd so that we don't have to do it later
//
//// Use all matches for ICP because:
//// Matches are normalized, then this->estimate_scale_and_assign_weight(); is called where only matches with nonzero weights are kept.
//// But normalization normalizes them before, then match sets are changed and unnormalization DOES NOT unnormalize all of them (because they got dropped -- they had zero weight).
////      if( w > 0.0 ) {
//        match->to_.push_back( q );
//        match->w_.push_back( w );
////      }      
//    }
//
//  }


  return scale;
}


template < unsigned int dim, unsigned int dof >
void rrl_estimation_ICP_matching_all<dim, dof>::weight_by_strength()
{

  //double sum_grad_w = 0.0;
  vcl_vector< double >  gradient_weights;
  gradient_weights.reserve( matches_.size() );
  for( typename vcl_vector< match_sptr_type >::size_type  i = 0; i < matches_.size(); ++i ) {
    feature_sptr_type  p = matches_[i]->from_;
    
    feature_sptr_type  q = matches_[i]->to_[0];

    if( matches_[i]->w_[0] != 0 ) { // this condition can be removed once we don't store points with zero weights
      //double grad_w = vcl_sqrt( p->strength_ * q->strength_ );
      //if( grad_w > 300.0 ) grad_w = 300.0;
      double grad_w = p->strength_ * q->strength_;
      if( grad_w > 90000.0 ) grad_w = 90000.0;
      matches_[i]->w_[0] *= grad_w;
      gradient_weights.push_back( grad_w );
      //sum_grad_w += grad_w;
    }

  }

  double grad_scale = rrel_util_median_abs_dev_scale( gradient_weights.begin(), gradient_weights.end() );
  vcl_cout << "Gradient scale: " << grad_scale << vcl_endl;
  for( typename vcl_vector< match_sptr_type >::size_type  i = 0; i < matches_.size(); ++i ) {

    matches_[i]->w_[0] *= ( 1.0 / grad_scale );  // ( 1.0 / sum_grad_w );
  }

}


template < unsigned int dim, unsigned int dof >
void rrl_estimation_ICP_matching_all<dim, dof>::weight_spatially()
{
  DEBUG( vcl_cout << "Weighting spatially..." << vcl_endl );
  vcl_vector< double >  spatial_weights;
  for( typename vcl_vector< match_sptr_type >::size_type  i = 0; i < matches_.size(); ++i ) {

    feature_sptr_type  q = matches_[i]->to_[0];

    // locations must be normalized, so center of the region is already at \vect{ 0 }
    double dist_to_center = q->location_.magnitude();
    double spatial_weight = 1.0 / ( 1.0 + 5.0*dist_to_center*dist_to_center );
    matches_[i]->w_[0] *= spatial_weight;
    spatial_weights.push_back( spatial_weight );

  }
  double spatial_scale = rrel_util_median_abs_dev_scale( spatial_weights.begin(), spatial_weights.end() );
  
  for( typename vcl_vector< match_sptr_type >::size_type  i = 0; i < matches_.size(); ++i ) {
  
    matches_[i]->w_[0] *= ( 1.0 / spatial_scale );
  }
  
}


// Estimate initial covariance matrix.
template < unsigned int dim, unsigned int dof >
void rrl_estimation_ICP_matching_all<dim, dof>::estimate_LS( bool estimate_parameters )
{
  // Normalize before estimation.
  //
  // Normalization parameters.
  vnl_vector_fixed< double, dim >  center_moving( 0.0 ), center_fixed( 0.0 );
  double avg_rad_moving( 0.0 ), avg_rad_fixed( 0.0 );

  if( normalize_matches_ ) {
    cdcl_normalize_matches( center_moving, avg_rad_moving, center_fixed, avg_rad_fixed );
    trans_->normalize( avg_rad_moving, avg_rad_fixed, center_moving, center_fixed );
  }

  double scale = this->estimate_scale_and_assign_weight();
  if( scale < 0.005 ) scale = 0.005; // in normalized coordinates



  if( weight_by_strength_ ) this->weight_by_strength();

  // Estimate initial matrix by setting up a least squares problem
  // A x = q
  // A^T A
  vnl_matrix_fixed< double, dof, dof >  AtA( 0.0 );
  vnl_matrix_fixed< double, dof, 1 >    Atb( 0.0 );
  
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

    // fixed point q
    feature_sptr_type  q = matches_[i]->to_[0];
    q_loc.set_column( 0, q->location_ );

    double weight = matches_[i]->w_[0];

    AtA += Jth.transpose() * q->error_projector_ * Jth * weight;
    Atb += Jth.transpose() * q->error_projector_ * q_loc * weight;

    sum_w += weight;

    //vcl_cout << Jth << vcl_endl << AtA << vcl_endl << vcl_endl;

  }
  //AtA *= ( 1.0 / sum_w );
  //Atb *= ( 1.0 / sum_w );

  vnl_svd< double > svd( AtA );
  vnl_matrix_fixed< double, dof, dof >  invAA = svd.inverse();

  // don't need to multiply by scale, since it is already included in the robust weights [yang:cvpr04]
  trans_->set_covariance( /*scale*scale **/ /** 2.0*/ /**/ /*scale_corners * scale_corners **/ invAA );

  if( estimate_parameters ) {
    vnl_matrix_fixed< double, dof, 1 >  new_params = invAA * Atb;
    //vcl_cout << "PARAMS: " << new_params << vcl_endl;
    trans_->set_parameterization( new_params.get_column( 0 ) );
  }

  if( normalize_matches_ ) {
    // convert transformation back to unnormalized coordinates
    trans_->unnormalize( avg_rad_moving, avg_rad_fixed, center_moving, center_fixed );
    this->cdcl_unnormalize_matches( center_moving, avg_rad_moving, center_fixed, avg_rad_fixed );
  }

  vcl_cout << "Initial covariance: " << vcl_endl << trans_->get_covariance() << vcl_endl << "Scale: " << scale << vcl_endl;

}


// Normalize moving and fixed points which the matches were formed with
template < unsigned int dim, unsigned int dof >
void rrl_estimation_ICP_matching_all<dim, dof>::cdcl_normalize_matches(
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

    // fixed point q
    feature_sptr_type  q = this->matches_[i]->to_[0];
    // if it can be inserted (i.e. not already there), add to center
    if( fixed_matched.insert( q ).second ) center_fixed += q->location_;
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
//    p->covariance_ /= (avg_radius_moving*avg_radius_moving);
  }

  // go through all matches (fixed points)
  for( typename vcl_set< feature_sptr_type >::iterator  it = fixed_matched.begin(); it != fixed_matched.end(); ++it ) {
    // fixed point q
    feature_sptr_type  q = *it;
    // convert locations to normalized coordinate system
    q->location_ /= avg_radius_fixed;

    // convert covariances to normalized coordinate system
//    q->covariance_ /= (avg_radius_fixed*avg_radius_fixed);
  }

}


// Unormalize moving and fixed points which the matches were formed with
template < unsigned int dim, unsigned int dof >
void rrl_estimation_ICP_matching_all<dim, dof>::cdcl_unnormalize_matches(
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
//    p->covariance_ *= (avg_radius_moving*avg_radius_moving);

    // fixed point q
    feature_sptr_type  q = this->matches_[i]->to_[0];

    // if it was inserted (i.e. not already there), unnormalize
    if( fixed_matched.insert( q ).second ) {
      // convert locations to unnormalized coordinate system
      q->location_ *= avg_radius_fixed;
      q->location_ += center_fixed;

    // convert covariances to unnormalized coordinate system
//        q->covariance_ *= (avg_radius_fixed*avg_radius_fixed);
    }
  }

}


#define RRL_ESTIMATION_ICP_MATCHING_ALL_INSTANTIATE( dim, dof ) \
template                                            \
class rrl_estimation_ICP_matching_all< dim, dof >;
