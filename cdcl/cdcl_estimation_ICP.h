#ifndef cdcl_estimation_ICP_h_
#define cdcl_estimation_ICP_h_

#include <vbl/vbl_ref_count.h>
#include <vbl/vbl_smart_ptr.h>
 
#include <rsdl/rsdl_kd_tree.h>

#include "cdcl_estimation_abs.h"
#include "cdcl_feature.h"
#include "cdcl_match.h"
#include "cdcl_trans.h"

//:
// \file
// \brief  Estimation handling the whole registration: 1) matching,
//         2) parameter estimate, 3) covariance estimate.
//         ICP with Euclidean distances for matching and normal distances
//         for estimation is used.
// \author Michal Sofka
// \date   May 2006


template < unsigned int dim, unsigned int dof >
class cdcl_estimation_ICP : public cdcl_estimation_abs< dim, dof >
{
public:

  typedef vbl_smart_ptr< cdcl_estimation_ICP >  sptr;

  typedef cdcl_estimation_abs< dim, dof >                  superclass_type;
  typedef typename superclass_type::feature_type           feature_type;
  typedef typename superclass_type::feature_sptr_type      feature_sptr_type;
  typedef typename superclass_type::match_sptr_type        match_sptr_type;
  typedef typename superclass_type::trans_type             trans_type;
  typedef typename superclass_type::trans_sptr_type        trans_sptr_type;
  typedef typename superclass_type::display_callback_type  display_callback_type;
 
  cdcl_estimation_ICP(): normalize_matches_( false ), weight_by_strength_( false ) {}

  // Construct estimation object given moving and fixed feature sets
  // and an initial transformation.
  cdcl_estimation_ICP( vcl_vector< feature_sptr_type > const &  moving, 
                       vcl_vector< feature_sptr_type > const &  fixed,
                       trans_sptr_type                 const &  trans );

  // Construct estimation object given moving and fixed multiresolution feature sets
  // and an initial transformation.
  // These need to have been normalized and have proper covariances.
  // This way, estimation object does not need to know anything about normalization.
  cdcl_estimation_ICP( vcl_vector< vcl_vector< feature_sptr_type > > const &  moving, 
                       vcl_vector< vcl_vector< feature_sptr_type > > const &  fixed,
                       trans_sptr_type                               const &  trans,
                       vcl_vector< double >                          const &  fixed_spacing );

  // Run estimation.
  virtual bool run( void *  caller = 0, 
                    display_callback_type  display_points = 0 );

  // Initialize estimation (compute initial weights and covariance).
  virtual void initialize();

  // Run one EM iteration of estimation (parameters, weights, covariance, weights).
  // Return true when converged.
  bool virtual one_iteration( void *  caller, 
                              display_callback_type  display_points,
                              unsigned int  iteration );

  // Set initial transform, will need to reestimate scale.
  void set_transform( trans_sptr_type const &  trans ) { trans_ = trans; scale_estimated_ = false; };

  // Get current transform.
  trans_sptr_type get_transform() { return trans_; };

  // Set normalize_matches_ to true if the matched points should be normalized before estimation.
  virtual void normalize_matches() { normalize_matches_ = true; };

  // Set moving data, useful when using different moving data, but same fixed data.
  // In that case, fixed kd-tree does not need to be created again.
  void set_moving_data( vcl_vector< vcl_vector< feature_sptr_type > > const &  moving ) { moving_ = moving; }

  // Compute RMS error using current matches.
  double RMS_error();

  // Compute weighted alingment error.
  double weighted_error();

  // Compute transfer error covariance using current set of matches.
  void compute_transfer_error_covariance();

  // Specify whether to weight matches by strength value of features.
  void weight_by_strength( bool weight_by_strength ) { weight_by_strength_ = weight_by_strength; }

protected:
  // Build kd trees, one for each resolution level of the fixed image.
  void build_kd_trees();

  // Estimate scale from current matches and assign robust weight.
  double estimate_scale_and_assign_weight();

  // Compute weights, use Euclidean distances.
  void find_closest_euclidean();

  // Estimate covariance matrix using least squares (normal distances between points).
  // Optionally estimate parameters as well.
  void estimate_LS( bool estimate_parameters = false );

  // Compute median error using current matches.
  double median_error();

  // Normalize moving and fixed points which the matches were formed with.
  void cdcl_normalize_matches( vnl_vector_fixed< typename feature_type::coord_type, dim >  & center_moving, 
                               typename feature_type::coord_type                           & avg_radius_moving, 
                               vnl_vector_fixed< typename feature_type::coord_type, dim >  & center_fixed, 
                               typename feature_type::coord_type                           & avg_radius_fixed );

  // Unormalize moving and fixed points which the matches were formed with.
  void cdcl_unnormalize_matches( vnl_vector_fixed< typename feature_type::coord_type, dim >                   const & center_moving, 
                                 typename feature_type::coord_type                                            const & avg_radius_moving, 
                                 vnl_vector_fixed< typename feature_type::coord_type, dim >                   const & center_fixed, 
                                 typename feature_type::coord_type                                            const & avg_radius_fixed );

  // Moving feature set, vector of fine-to-coarse resolutions.
  vcl_vector< vcl_vector< feature_sptr_type > >  moving_;

  // Fixed feature set, vector of fine-to-coarse resolutions.
  vcl_vector< vcl_vector< feature_sptr_type > >  fixed_;

  // Current set of matches.
  vcl_vector< match_sptr_type >  matches_;

  // Current transformation estimate.
  trans_sptr_type  trans_;

  // Kd-tree formed from fixed feature points, vector of fine-to-coarse resolutions.
  vcl_vector< vbl_smart_ptr< rsdl_kd_tree > >  kd_tree_;

  // Current resolution level.
  unsigned int res_;

  // Maximum resolution level (coarsest).
  unsigned int max_res_;

  // Spacing between points in the fixed dataset, used for switching resolutions.
  vcl_vector< double >  fixed_spacing_;

  // True if the matched points should be normalized before estimation.
  bool normalize_matches_;

  // Weight the matches by the feature strength.
  bool weight_by_strength_;

  // We have estimated scale already (in the first iteration).
  bool scale_estimated_;

};




#endif
