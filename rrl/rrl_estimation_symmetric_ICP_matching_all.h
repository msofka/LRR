#ifndef rrl_estimation_symmetric_ICP_matching_all_h_
#define rrl_estimation_symmetric_ICP_matching_all_h_

#include <cdcl/cdcl_estimation_symmetric_ICP_matching_all.h>
#include <rrl/rrl_estimation_ICP_matching_all.h>

#include "itkMultiThreader.h"


//:
// \file
// \brief  Symmetric estimation.
//         The idea is to have forward and backward estimation object.
//         Constraints from inverse transform are included in the estimation.
//         This is the same as cdcl_estimation_symmetric_ICP but with 
//         rrl_estimation_symmetric_ICP_matching_all as the estimation type.
// \author Michal Sofka
// \date   Feb 2008



template < unsigned int dim, unsigned int dof >
class rrl_estimation_symmetric_ICP_matching_all : public cdcl_estimation_symmetric_ICP_matching_all< dim, dof >
{
public:
  typedef vbl_smart_ptr< rrl_estimation_symmetric_ICP_matching_all >  sptr;
  typedef rrl_estimation_ICP_matching_all< dim, dof >         estimation_type;

  typedef rrl_estimation_symmetric_ICP_matching_all< dim, dof >  superclass_type;

  typedef typename estimation_type::feature_type           feature_type;
  typedef typename estimation_type::feature_sptr_type      feature_sptr_type;
  typedef typename estimation_type::match_type             match_type;
  typedef typename estimation_type::match_sptr_type        match_sptr_type;
  typedef typename estimation_type::trans_sptr_type        trans_sptr_type;
  typedef typename feature_type::location_type             location_type;
  typedef typename feature_type::coord_type                coord_type;
  //typedef typename estimation_type::display_callback_type  display_callback_type;

  typedef itk::MultiThreader  ThreaderType;


  typedef  unsigned int  VoronoiMapPixelType;
  typedef itk::Image< VoronoiMapPixelType, 3 >   VoronoiMapImageType;
 
  typedef void   (*display_callback_type) ( void *  caller, 
                                            vcl_vector< feature_sptr_type > const &  moving_,
                                            vcl_vector< feature_sptr_type > const &  fixed_,
                                            vcl_vector< match_sptr_type >   const &  matches_,
                                            trans_sptr_type                 const &  transform,
                                            unsigned int                             iteration );

  rrl_estimation_symmetric_ICP_matching_all() {}

  // Construct estimation object given moving and fixed feature sets,
  // and an initial transformation. All movind and fixed points are used.
  rrl_estimation_symmetric_ICP_matching_all( vcl_vector< feature_sptr_type >  const &  moving, 
                             vcl_vector< feature_sptr_type >  const &  fixed,
                             trans_sptr_type                  const &  trans,
                             trans_sptr_type                  const &  trans_inverse );

  // Construct estimation object given moving and fixed feature sets,
  // an initial transformation, and moving region corners.
	// Constructed will be false if there are not enough points in either of the regions.
  rrl_estimation_symmetric_ICP_matching_all( vcl_vector< feature_sptr_type >  const &  moving, 
                             vcl_vector< feature_sptr_type >  const &  fixed,
                             trans_sptr_type                  const &  trans,
                             trans_sptr_type                  const &  trans_inverse,
                             location_type  const &  moving_x0,
                             location_type  const &  moving_x1,
                             location_type  const &  fixed_x0,
                             location_type  const &  fixed_x1,
														 bool constructed );


  // Construct estimation object given moving and fixed feature sets,
  // and features in fixe/moving image that will be used to match against moving/fixed image.
  rrl_estimation_symmetric_ICP_matching_all( vcl_vector< feature_sptr_type >  const &  moving, 
                                              vcl_vector< feature_sptr_type >  const &  fixed,
                                              trans_sptr_type                  const &  trans,
                                              trans_sptr_type                  const &  trans_inverse,
                                              vcl_vector< feature_sptr_type >  const &  moving_inside, 
                                              vcl_vector< feature_sptr_type >  const &  fixed_inside );

  // Set the current moving and fixed regions. Only features withing the regions are used in matching.
  // Note that we are passing in the features again as the fixed set is separated based on shape in cdcl_estimation_ICP_matching_all
  // and the original set of features is not available.
  bool set_regions( vcl_vector< feature_sptr_type >  const &  moving, 
                    vcl_vector< feature_sptr_type >  const &  fixed,
                    location_type  const &  moving_x0,
                    location_type  const &  moving_x1,
                    location_type  const &  fixed_x0,
                    location_type  const &  fixed_x1 );

  // Set features from fixe/moving image which will be used to match against moving/fixed image.
  void set_features( vcl_vector< feature_sptr_type >  const &  moving_inside, 
                     vcl_vector< feature_sptr_type >  const &  fixed_inside );

  // Set voronoi maps to forward and backward estimators to speed up matching.
  void set_voronoi_maps( VoronoiMapImageType::Pointer const &  moving_voronoi_map,
    VoronoiMapImageType::Pointer const &  fixed_voronoi_map );

  // Estimate parameters. Return true when converged.
  virtual bool estimate_parameters() { return false;/*to be implemented*/};

  // Estimate covariance. Return true when converged.
  virtual bool estimate_covariance() { return false;/*to be implemented*/};

  // Run estimation.
  // HACK: This is the same as the parent, but the display_callback_type is different.
  virtual bool run( void *  caller = 0, 
    typename cdcl_estimation_abs< dim, dof, cdcl_feature<dim> >::display_callback_type  display_points = 0 ) { return false; };

  // Initialize estimation (compute initial weights and covariance).
  virtual void initialize();

  // Compute weights, use Euclidean distances.
  void find_closest_euclidean();

  //void estimate_scale_and_assign_weight() { forward_->estimate_scale_and_assign_weight(); backward_->estimate_scale_and_assign_weight(); }

  // Estimate covariance matrix using least squares (normal distances between points).
  // Optionally estimate parameters as well.
  void estimate_LS( bool estimate_parameters );

  // Run one EM iteration of estimation (parameters, weights, covariance, weights).
  // Return true when converged.
  virtual bool one_iteration( void *  caller, 
                              display_callback_type  display_points,
                              unsigned int  iteration );

  // Set initial transform.
  void set_transform( trans_sptr_type const &  trans );

  // Get current transform.
  trans_sptr_type get_transform() { return forward_->get_transform(); };

  // Get current backward transform.
  // Usually this is used when covariance of the backward transform has been estimated (by estimate_LS_backward( false ); )
  trans_sptr_type get_transform_backward() { return backward_->get_transform(); };

  // Specify whether to weight matches by strength value of features.
  void weight_by_strength( bool weight_by_strength ) { forward_->weight_by_strength( weight_by_strength );
                                                       backward_->weight_by_strength( weight_by_strength );
                                                       /*weight_by_strength_ = weight_by_strength; */}

  // Set normalize_matches_ to true if the matched points should be normalized before estimation.
  virtual void normalize_matches() { forward_->normalize_matches();
                                     backward_->normalize_matches();
                                     normalize_matches_ = true; };

  // Compute RMS error using current matches.
  double RMS_error() { return forward_->RMS_error(); };

  // Return true if estimation oscillates.
  bool oscillated() { return oscillation_count_ > 10; };

  // Compute weighted alingment error.
  double weighted_error() { return forward_->weighted_error(); };

  // Compute transfer error covariance using current set of matches.
  void compute_transfer_error_covariance() { forward_->compute_transfer_error_covariance();  backward_->compute_transfer_error_covariance(); };

  // Estimate backward transform using symmetric matching (used for covariance estimation).
  void estimate_LS_backward( bool estimate_parameters );

  // Agreement of the angles using SHEET features and their normals.
  double sheet_angles() { return forward_->sheet_angles(); };

  // Agreement of the angles using TUBE features and their normals.
  double tube_angles() { return forward_->tube_angles(); };


  // Compute weighted alingment error.
  double weighted_error_backward() { return backward_->weighted_error(); };

  // Agreement of the angles using SHEET features and their normals.
  double sheet_angles_backward() { return backward_->sheet_angles(); };

  // Agreement of the angles using TUBE features and their normals.
  double tube_angles_backward() { return backward_->tube_angles(); };

protected:

  // Normalize moving and fixed points which the matches were formed with
  // Use supplied known centers and radii
  void cdcl_normalize_matches_known( vcl_vector< match_sptr_type >        & matches,
                                     vnl_vector_fixed< typename feature_type::coord_type, dim >                   const & center_moving, 
                                     typename feature_type::coord_type                                            const & avg_radius_moving, 
                                     vnl_vector_fixed< typename feature_type::coord_type, dim >                   const & center_fixed, 
                                     typename feature_type::coord_type                                            const & avg_radius_fixed );

  // Callback for multithreader to compute AtA matrix during estimation.
  // Static function used as a "callback" by the MultiThreader.  The threading
  // library will call this routine for each thread.
  static ITK_THREAD_RETURN_TYPE  ThreaderCallbackForward( void * arg );
  static ITK_THREAD_RETURN_TYPE  ThreaderCallbackBackward( void * arg );

  void ThreadedGenerateDataForward( unsigned int ThreadId );
  void ThreadedGenerateDataBackward( unsigned int ThreadId );

  // Forward estimation object, estimation of forward matches and covariances.
  vbl_smart_ptr< estimation_type >  forward_;

  // Backward estimation object, estimation of forward matches and covariances.
  vbl_smart_ptr< estimation_type >  backward_;

  // True if the matched points should be normalized before estimation.
  bool normalize_matches_;

  // Current number of oscillations, used for convergence testing.
  unsigned int oscillation_count_;

  // Previous error difference (used for oscillation counting).
  double error_difference_;

  // Current weighted error.
  double weighted_error_;

  // Current number of matches used in estimation (multiresolution).
  unsigned int number_matches_;

  // Estimating at the finest level of the multiresolution (using maximum number of available matches).
  bool finest_level_;

  // Indices of the matches_ array separating it into chunks for multithreading.
  std::vector< typename vcl_vector< match_sptr_type >::size_type >  m_Indices;

  // Multithreading of the estimation.
  ThreaderType::Pointer  m_Threader;

  vcl_vector< vnl_matrix_fixed< coord_type, dof, dof > >  AtAs;
  vcl_vector< vnl_matrix_fixed< coord_type, dof, 1 > >    Atbs;

  // Internal structure used for passing image data into the threading library.
  struct ThreadStruct
  {
    sptr  estimation;
  };

};




#endif
