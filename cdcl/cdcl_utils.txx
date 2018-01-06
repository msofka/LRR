
#ifndef cdcl_utils_txx_
#define cdcl_utils_txx_

#include <vcl_iostream.h>
#include <vcl_fstream.h>
#include <vcl_iomanip.h>
#include <vcl_sstream.h>
#include <vcl_algorithm.h>
#include <vcl_set.h>
#include <vcl_vector.h>
#include <vcl_iterator.h>

#include <rsdl/rsdl_dist.h>
#include <rsdl/rsdl_kd_tree.h>

#include <vnl/vnl_math.h>
#include <vnl/vnl_random.h>
#include <vnl/vnl_transpose.h>
#include <vnl/vnl_vector_fixed.h>
#include <vnl/algo/vnl_symmetric_eigensystem.h>

#include "cdcl_utils.h"

//:
// \file
// \brief  Utility functions for cdcl class.
// \author Michal Sofka
// \date   May 2006


// Parse raw matrix data into a feature set.
template < unsigned int dim >
int cdcl_parse_raw_data( vnl_matrix< typename cdcl_feature< dim >::coord_type >  const &  raw_matrix, 
                         vcl_vector< vbl_smart_ptr< cdcl_feature< dim > > > &  feature_set )
{
  typedef typename cdcl_feature<dim>::coord_type  coord_type;

  feature_set.clear();
  // parse the raw matrix into feature set
  if( raw_matrix.cols() == 5 || raw_matrix.cols() == 9 ) {
    for( unsigned  r = 0; r < raw_matrix.rows(); ++r ) {
      // parse location
      vnl_vector_fixed< coord_type, dim >  location;
      for( int d = 0; d < dim; ++d ) {
        location[d] = raw_matrix( r, d );
      }

      // parse covariance
      vnl_matrix_fixed< coord_type, dim, dim >  covariance;
      unsigned int ind = dim;
      for( unsigned int d1 = 0; d1 < dim; ++d1 )
        for( unsigned int d2 = d1; d2 < dim; ++d2 ) {
          covariance( d2, d1 ) = covariance( d1, d2 ) = raw_matrix( r, ind );
          ++ind;
        }

      //vcl_cout << "location: " << vcl_endl << location << vcl_endl
      //         << "covariance: " << vcl_endl << covariance << vcl_endl;

      // create feature vector
      feature_set.push_back( new cdcl_feature< dim >( location, covariance ) );
    }
    // compute location covariances, if these lines are included, covariances read from the file are discarded
    // be careful about this since the covariances in a file might have been computed from images, while here only feature points are used
    vnl_vector_fixed< coord_type, dim >  dummy_min_pt, dummy_max_pt;
    double mean_dist = data_spacing_and_bounding_box( feature_set, dummy_min_pt, dummy_max_pt );
    vcl_vector< typename vcl_vector< vbl_smart_ptr< cdcl_feature< dim > > >::size_type >  dummy;
    location_covariance( feature_set, mean_dist, dummy );
  }
  else if( raw_matrix.cols() == 2 || raw_matrix.cols() == 3 ) { // no covariance available from file
    for( unsigned  r = 0; r < raw_matrix.rows(); ++r ) {
      // parse location
      vnl_vector_fixed< coord_type, dim >  location;
      for( unsigned int d = 0; d < dim; ++d ) {
        location[d] = raw_matrix( r, d );
      }

      // initialize covariance to 0
      vnl_matrix_fixed< coord_type, dim, dim >  covariance( 0.0 );

      // create feature vector
      feature_set.push_back( new cdcl_feature< dim >( location, covariance ) );
    }
    // compute location covariances
    vnl_vector_fixed< coord_type, dim >  dummy_min_pt, dummy_max_pt;
    double mean_dist = data_spacing_and_bounding_box( feature_set, dummy_min_pt, dummy_max_pt );

    vcl_vector< typename vcl_vector< vbl_smart_ptr< cdcl_feature< dim > > >::size_type >  dummy;
    location_covariance( feature_set, mean_dist, dummy );
  }
  else {
    vcl_cerr << "Error: Don't know how to parse data file with " << raw_matrix.cols() << " columns." << vcl_endl;
    return 1;
  }

  return 0;
}


// Compute center and average radius of a feature set.
template < unsigned int dim >
void cdcl_radius_center( vcl_vector< vbl_smart_ptr< cdcl_feature< dim > > >  const &  feature_set,
                         vnl_vector_fixed< typename cdcl_feature<dim>::coord_type, dim >                           &  center,
                         typename cdcl_feature<dim>::coord_type                                                    &  avg_radius )
{
  // compute data center of mass
  typedef typename vcl_vector< typename cdcl_feature< dim >::sptr >::size_type  feature_set_size_type;
  center.fill( 0.0 );
  for( feature_set_size_type i = 0; i < feature_set.size(); ++i ) {
    center += feature_set[i]->location_;
  }
  center /= feature_set.size();

  // center the data and compute average radius
  avg_radius = 0.0;
  for( feature_set_size_type  i = 0; i < feature_set.size(); ++i ) {
    feature_set[i]->location_ -= center;
    avg_radius += feature_set[i]->location_.magnitude();
  }
  avg_radius /= feature_set.size();
}


// Normalize data in a feature set to have location with zero mean and avg. radius of one.
template < unsigned int dim >
void cdcl_normalize_data( vcl_vector< vbl_smart_ptr< cdcl_feature< dim > > >  &  feature_set,
                          vnl_vector_fixed< typename cdcl_feature<dim>::coord_type, dim >                     &  center,
                          typename cdcl_feature<dim>::coord_type                                              &  avg_radius )
{
  // compute data center of mass
  typedef typename vcl_vector< typename cdcl_feature< dim >::sptr >::size_type  feature_set_size_type;
  center.fill( 0.0 );
  for( feature_set_size_type i = 0; i < feature_set.size(); ++i ) {
    center += feature_set[i]->location_;
  }
  center /= feature_set.size();

  // center the data and compute average radius
  avg_radius = 0.0;
  for( feature_set_size_type  i = 0; i < feature_set.size(); ++i ) {
    feature_set[i]->location_ -= center;
    avg_radius += feature_set[i]->location_.magnitude();
  }
  avg_radius /= feature_set.size();

  // assign normalized data
  for( feature_set_size_type  i = 0; i < feature_set.size(); ++i ) {
    // convert locations to normalized coordinate system
    feature_set[i]->location_ /= avg_radius;

    // convert covariances to normalized coordinate system
    feature_set[i]->covariance_ /= (avg_radius*avg_radius);
  }
}


// Normalize data in a feature set to have location with zero mean and avg. radius of one.
// Use given radius.
template < unsigned int dim >
void cdcl_normalize_data_known_radius( vcl_vector< vbl_smart_ptr< cdcl_feature< dim > > >        &  feature_set,
                                       vnl_vector_fixed< typename cdcl_feature<dim>::coord_type, dim >                           &  center,
                                       typename cdcl_feature<dim>::coord_type                                              const &  avg_radius )
{
  // compute data center of mass
  typedef typename vcl_vector< typename cdcl_feature< dim >::sptr >::size_type  feature_set_size_type;
  center.fill( 0.0 );
  for( feature_set_size_type i = 0; i < feature_set.size(); ++i ) {
    center += feature_set[i]->location_;
  }
  center /= feature_set.size();

  // assign normalized data
  for( feature_set_size_type  i = 0; i < feature_set.size(); ++i ) {
    // center the data
    feature_set[i]->location_ -= center;
    // convert locations to normalized coordinate system
    feature_set[i]->location_ /= avg_radius;

    // convert covariances to normalized coordinate system
    feature_set[i]->covariance_ /= (avg_radius*avg_radius);
  }
}


// Normalize data in a feature set to have location with zero mean and avg. radius of one.
// Use given center and radius.
template < unsigned int dim >
void cdcl_normalize_data_known( vcl_vector< vbl_smart_ptr< cdcl_feature< dim > > >        &  feature_set,
                                vnl_vector_fixed< typename cdcl_feature<dim>::coord_type, dim >                     const &  center,
                                typename cdcl_feature<dim>::coord_type                                              const &  avg_radius )
{
  // assign normalized data
  typedef typename vcl_vector< typename cdcl_feature< dim >::sptr >::size_type  feature_set_size_type;
  for( feature_set_size_type  i = 0; i < feature_set.size(); ++i ) {
    // center the data
    feature_set[i]->location_ -= center;
    // convert locations to normalized coordinate system
    feature_set[i]->location_ /= avg_radius;

    // convert covariances to normalized coordinate system
    feature_set[i]->covariance_ /= (avg_radius*avg_radius);
  }
}


// Unnormalize data in a feature set.
template < unsigned int dim >
void cdcl_unnormalize_data( vcl_vector< vbl_smart_ptr< cdcl_feature< dim > > > const &  feature_set,
                            vnl_vector_fixed< typename cdcl_feature<dim>::coord_type, dim >                      const &  center,
                            typename cdcl_feature<dim>::coord_type                                               const &  avg_radius )
{
  typedef typename vcl_vector< typename cdcl_feature< dim >::sptr >::size_type  feature_set_size_type;

  // assign unnormalized data
  for( feature_set_size_type  i = 0; i < feature_set.size(); ++i ) {
    // convert locations to normalized coordinate system
    feature_set[i]->location_ = feature_set[i]->location_ * avg_radius + center;

    // convert covariances to unnormalized coordinate system
    feature_set[i]->covariance_ *= (avg_radius*avg_radius);
  }
}


// Normalize moving and fixed points which the matches were formed with
template < unsigned int dim >
void cdcl_normalize_matches( vcl_vector< vbl_smart_ptr< cdcl_match< dim > > >  & matches, 
                             vnl_vector_fixed< typename cdcl_feature<dim>::coord_type, dim >                   & center_moving, 
                             typename cdcl_feature<dim>::coord_type                                            & avg_radius_moving, 
                             vnl_vector_fixed< typename cdcl_feature<dim>::coord_type, dim >                   & center_fixed, 
                             typename cdcl_feature<dim>::coord_type                                            & avg_radius_fixed )
{
  typedef typename cdcl_feature< dim >::sptr     feature_sptr_type;
  typedef typename cdcl_match< dim >::sptr       match_sptr_type;

  center_moving.fill( 0.0 );
  avg_radius_moving = 0.0;
  center_fixed.fill( 0.0 );
  avg_radius_fixed = 0.0;

  // set of fixed points that got matched to moving points
  // utilize unique container concept, i.e. no two points are identical
  // we need to make sure that points don't get normalized twice
  vcl_set< vbl_smart_ptr< cdcl_feature< dim > > >  fixed_matched;

  // compute data center of mass
  //
  // go through all moving points
  for( typename vcl_vector< match_sptr_type >::size_type  i = 0; i < matches.size(); ++i ) {
    // moving point p
    feature_sptr_type  p = matches[i]->from_;
    center_moving += p->location_;

    // go through all matches (fixed points)
    for( typename vcl_vector< match_sptr_type >::size_type  j = 0; j < matches[i]->to_.size(); ++j ) {
      // fixed point q
      feature_sptr_type  q = matches[i]->to_[j];
      // if it can be inserted (i.e. not already there), add to center
      if( fixed_matched.insert( q ).second ) center_fixed += q->location_;
    }
  }
  center_moving /= matches.size();
  center_fixed /= fixed_matched.size();


  // center the data and compute average radius
  //
  vnl_vector_fixed< typename cdcl_feature<dim>::coord_type, dim >  centered;
  for( typename vcl_vector< match_sptr_type >::size_type  i = 0; i < matches.size(); ++i ) {
    // moving point p
    feature_sptr_type  p = matches[i]->from_;
    p->location_ -= center_moving;
    avg_radius_moving += p->location_.magnitude();
  }
  avg_radius_moving /= matches.size();

  // go through all matches (fixed points)
  for( typename vcl_set< vbl_smart_ptr< cdcl_feature< dim > > >::iterator  it = fixed_matched.begin(); it != fixed_matched.end(); ++it ) {
    // fixed point q
    feature_sptr_type  q = *it;
    q->location_ -= center_fixed;
    avg_radius_fixed += q->location_.magnitude();
  }
  avg_radius_fixed /= fixed_matched.size();;


  // assign normalized data
  //
  for( typename vcl_vector< match_sptr_type >::size_type  i = 0; i < matches.size(); ++i ) {
    // moving point p
    feature_sptr_type  p = matches[i]->from_;

    // convert locations to normalized coordinate system
    p->location_ /= avg_radius_moving;

    // convert covariances to normalized coordinate system
    p->covariance_ /= (avg_radius_moving*avg_radius_moving);
  }

  // go through all matches (fixed points)
  for( typename vcl_set< vbl_smart_ptr< cdcl_feature< dim > > >::iterator  it = fixed_matched.begin(); it != fixed_matched.end(); ++it ) {
    // fixed point q
    feature_sptr_type  q = *it;
    // convert locations to normalized coordinate system
    q->location_ /= avg_radius_fixed;

    // convert covariances to normalized coordinate system
    q->covariance_ /= (avg_radius_fixed*avg_radius_fixed);
  }

}

// Normalize moving and fixed points which the matches were formed with
// Use supplied known centers and radii
template < unsigned int dim >
void cdcl_normalize_matches_known( vcl_vector< vbl_smart_ptr< cdcl_match< dim > > >        & matches, 
                                   vnl_vector_fixed< typename cdcl_feature<dim>::coord_type, dim >                   const & center_moving, 
                                   typename cdcl_feature<dim>::coord_type                                            const & avg_radius_moving, 
                                   vnl_vector_fixed< typename cdcl_feature<dim>::coord_type, dim >                   const & center_fixed, 
                                   typename cdcl_feature<dim>::coord_type                                            const & avg_radius_fixed )
{
  typedef typename cdcl_feature< dim >::sptr     feature_sptr_type;
  typedef typename cdcl_match< dim >::sptr       match_sptr_type;

  // set of fixed points that got matched to moving points
  // utilize unique container concept, i.e. no two points are identical
  // we need to make sure that points don't get normalized twice
  vcl_set< vbl_smart_ptr< cdcl_feature< dim > > >  fixed_matched;

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


// Unormalize moving and fixed points which the matches were formed with
template < unsigned int dim >
void cdcl_unnormalize_matches( vcl_vector< vbl_smart_ptr< cdcl_match< dim > > >        & matches, 
                               vnl_vector_fixed< typename cdcl_feature<dim>::coord_type, dim >                   const & center_moving, 
                               typename cdcl_feature<dim>::coord_type                                            const & avg_radius_moving, 
                               vnl_vector_fixed< typename cdcl_feature<dim>::coord_type, dim >                   const & center_fixed, 
                               typename cdcl_feature<dim>::coord_type                                            const & avg_radius_fixed )
{
  typedef typename cdcl_feature< dim >::sptr     feature_sptr_type;
  typedef typename cdcl_match< dim >::sptr       match_sptr_type;

  // set of fixed points that got matched to moving points
  // utilize unique container concept, i.e. no two points are identical
  // we need to make sure that points don't get unnormalized twice
  vcl_set< vbl_smart_ptr< cdcl_feature< dim > > >  fixed_matched;

  // remove radius normalization and centering of the data
  //
  for( typename vcl_vector< match_sptr_type >::size_type  i = 0; i < matches.size(); ++i ) {
    // moving point p
    feature_sptr_type  p = matches[i]->from_;

    // convert locations to unnormalized coordinate system
    p->location_ *= avg_radius_moving;
    p->location_ += center_moving;

    // convert covariances to unnormalized coordinate system
    p->covariance_ *= (avg_radius_moving*avg_radius_moving);

    // go through all matches (fixed points)
    for( typename vcl_vector< match_sptr_type >::size_type  j = 0; j < matches[i]->to_.size(); ++j ) {
      // fixed point q
      feature_sptr_type  q = matches[i]->to_[j];

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


namespace {
// Compute eigenvalues, sorted smallest to largest.
inline void cdcl_evalues( vnl_matrix_fixed< cdcl_feature<2>::coord_type, 2, 2 > const &  C, vnl_vector_fixed< cdcl_feature<2>::coord_type, 2 > & evalues )
{
  cdcl_feature<2>::coord_type temp = (( C(0,0) - C(1,1) )*( C(0,0) - C(1,1) ) + 4*C(0,1)*C(0,1) );
  cdcl_feature<2>::coord_type sqrt_temp = vcl_sqrt( temp );

  cdcl_feature<2>::coord_type Lpp = (C(0,0) + C(1,1) - sqrt_temp)/2.0;
  cdcl_feature<2>::coord_type Lqq = (C(0,0) + C(1,1) + sqrt_temp)/2.0;

  // sort the eigenvalues: lambda1 < lambda2
  if( Lpp < Lqq ) {
    evalues[0] = Lpp;
    evalues[1] = Lqq;
  }
  else {
    evalues[0] = Lqq;
    evalues[1] = Lpp;
  }
}

// Compute eigenvalues, sorted smallest to largest.
inline void cdcl_evalues( vnl_matrix_fixed< cdcl_feature<3>::coord_type, 3, 3 > const &  C, vnl_vector_fixed< cdcl_feature<3>::coord_type, 3 > & evalues )
{
  double l1, l2, l3;
  vnl_symmetric_eigensystem_compute_eigenvals( C( 0, 0 ), C( 0, 1 ), C( 0, 2 ),
                                                          C( 1, 1 ), C( 1, 2 ),
                                                                     C( 2, 2 ),
                                               l1, l2, l3 );
  evalues[0] = l1; // smallest
  evalues[2] = l2;
  evalues[3] = l3; // largest
}

}


// Compute location covariance given a set of points.
// Pass in spacing between points to compute radius in which the covariance is computed around each point.
template < unsigned int dim >
void location_covariance( vcl_vector< vbl_smart_ptr< cdcl_feature< dim > > >                                    &  feature_points,
                          const typename cdcl_feature<dim>::coord_type                                          &  spacing,
                          vcl_vector< typename vcl_vector< vbl_smart_ptr< cdcl_feature< dim > > >::size_type >  &  indices )
{
  //vcl_cout << "Computing location covariance for " << feature_points.size() << " points with spacing " << spacing << vcl_endl;
  typedef typename cdcl_feature< dim >::sptr  feature_sptr_type;
  typedef typename vcl_vector< feature_sptr_type >::size_type  feature_set_size_type;
  typedef typename cdcl_feature<dim>::coord_type  coord_type;

  // build kd_tree from points
  vcl_vector< rsdl_point > points;
  rsdl_point  pt( dim, 0 );
  for( feature_set_size_type i = 0; i < feature_points.size(); ++i ) {
    //pt.set_cartesian( feature_points[i]->location_ );
    for( unsigned int d = 0; d < dim; ++d ) pt.cartesian( d ) = feature_points[i]->location_[d];
    points.push_back( pt );
  }

  vbl_smart_ptr< rsdl_kd_tree >  kd_tree = new rsdl_kd_tree( points );
  double sigma, radius;
  if( dim == 3 ) {  // for bunny
    sigma  = 1.5*spacing;
    radius = 1.5*sigma;
  }
  else {  // for H example
    sigma  = 1.5*spacing;
    radius = 1.5*sigma;
  }
 
  vcl_vector< feature_set_size_type >  noise_indices;
  for( feature_set_size_type i = 0; i < feature_points.size(); ++i ) {
    // query kd-tree and find set of points that are in a neighborhood
    rsdl_point  query_pt( dim, 0 );
    //query_pt.set_cartesian( feature_points[i]->location_ );
    for( unsigned int d = 0; d < dim; ++d ) query_pt.cartesian( d ) = feature_points[i]->location_[d];


    vcl_vector< rsdl_point >  points;
    vcl_vector< int >  indices;
    kd_tree->points_in_radius( query_pt, radius, points, indices );

    //assert( indices.size() > 0 );
    if( indices.size() <= 1 ) { // we have most probably found a lonely noise point (if size is 1, then only 1 point found (equal to query point), no neighbors)
      noise_indices.push_back( i );
      continue;
    }

    // compute mean location
    vnl_vector_fixed< coord_type, dim >  mean_location( 0.0 );
    for( typename vcl_vector< int >::size_type  j = 0; j < indices.size(); ++j ) {
      // current point
      feature_sptr_type  curr = feature_points[indices[j]];
      mean_location += curr->location_;
    }
    mean_location /= indices.size();

    // compute covariance as sum of outer products of centered coordinates
    vnl_matrix_fixed< typename cdcl_feature<dim>::coord_type, dim, dim >  covar( 0.0 );
    double gauss_weights_sum = 0.0;
    for( typename vcl_vector< int >::size_type  j = 0; j < indices.size(); ++j ) {
      // current point
      feature_sptr_type  curr = feature_points[indices[j]];
      vnl_vector_fixed< coord_type, dim >  centered = curr->location_ - mean_location;
      double weight = cdcl_gauss_weight( centered, sigma );
      covar += coord_type( weight ) * outer_product( centered, centered );
      gauss_weights_sum += weight;
    }
    covar *= ( 1.0 / gauss_weights_sum );

    //if( cdcl_largest_evalue( covar ) > 900.0 * spacing ) {  // std. dev. must be smaller than 30
    //  noise_indices.push_back( i );
    //  continue;
    //}

    feature_points[i]->covariance_ = covar;
  }

  if( noise_indices.size() > 0 ) {
#define RANDOM_COV_FOR_NOISE 0
#if RANDOM_COV_FOR_NOISE
    vcl_cerr << "There were " << noise_indices.size() << " potential noise points found, their covariance will be random ..." << vcl_endl;
    // compute average largest and smallest eigenvalues
    double  avg_sqrt_evalue = 0.0;
    vnl_vector_fixed< double, dim >  curr_evalue( 0.0 );
    for( feature_set_size_type i = 0; i < feature_points.size(); ++i ) {
      cdcl_evalues( feature_points[i]->covariance_, curr_evalue );
      for( unsigned int d = 0; d < dim; ++d )
        if( curr_evalue[d] < 0.0 ) vcl_cout << curr_evalue[d] << vcl_endl << feature_points[i]->covariance_ << vcl_endl;
        else avg_sqrt_evalue += vcl_sqrt( curr_evalue[d] );
    }
    avg_sqrt_evalue /= ( ( feature_points.size() - noise_indices.size() ) * dim );  // subtract noise indices size which currently have zero covariance

    vnl_random  random;
    for( typename vcl_vector< feature_set_size_type >::size_type  i = 0; i < noise_indices.size(); ++i ) {
      // random eigenvalue based on average eigenvalues, insert into a diagonal matrix
      vnl_matrix_fixed< double, dim, dim >  avg_values_diag( 0.0 );
      for( unsigned int d = 0; d < dim; ++d ) {
        // get 2*sigma around current sqrt of eigenvalues for covariance
        // average will still stay the same if taken between 0 and 2.0*sqrt(eigenvalue)
        double new_sigma = random.drand32( 1e-5, 2.0*avg_sqrt_evalue );
        avg_values_diag( d, d ) = new_sigma*new_sigma;
      }
      vcl_cout << vcl_endl;

      // angle defining orientation of the local coordinate system
      double ang = random.drand32( -vnl_math::pi, vnl_math::pi );

      // matrix of eigenvectors
      vnl_matrix_fixed< double, dim, dim >  evectors;
      if( dim == 2 ) {
        evectors( 0, 0 ) =  vcl_cos( ang );
        evectors( 0, 1 ) = -vcl_sin( ang );
        evectors( 1, 0 ) =  vcl_sin( ang );
        evectors( 1, 1 ) =  vcl_cos( ang );
      }
      else {
        vcl_cout << "Error: Missing implementation: add eigenvectors for 3 dimensions." << vcl_endl;
        vcl_abort();
      }
      // construct by recomposing eigen decomposition
      feature_points[noise_indices[i]]->covariance_ = evectors * avg_values_diag * evectors.transpose();
    }
#else
    vcl_cerr << "There were " << noise_indices.size() << " potential noise points found, removing them from the point set ... " << vcl_endl;
    typename vcl_vector< feature_sptr_type >  features_wout_noise;
    // copy only non-noise points
    noise_indices.push_back( feature_points.size() ); // add the last index
    typename vcl_vector< feature_set_size_type >::size_type  start_from = 0;
    for( typename vcl_vector< feature_set_size_type >::size_type  i = 0; i < noise_indices.size(); ++i ) {
      for( feature_set_size_type j = start_from; j < noise_indices[i]; ++j ) {
        features_wout_noise.push_back( feature_points[j] );
      }
      start_from = noise_indices[i]+1;
    }
    // assign non-noise points to the original feature vector
//commented this out because indices to points_for_cov in subsample_fine_covariances are w.r.t. original dataset (wout removing any points here)
//    feature_points = features_wout_noise;

    noise_indices.erase( noise_indices.begin() + noise_indices.size() - 1 ); // remove the last index
    indices = noise_indices;
#endif

  }

}


// Compute data spacing as a mean distance between points.
// Compute bounding box of all data points along the way.
template < unsigned int dim >
double data_spacing_and_bounding_box( vcl_vector< vbl_smart_ptr< cdcl_feature< dim > > >  const &  feature_points,
                                      vnl_vector_fixed< typename cdcl_feature<dim>::coord_type, dim >                         &  min_pt,
                                      vnl_vector_fixed< typename cdcl_feature<dim>::coord_type, dim >                         &  max_pt )
{
  typedef typename cdcl_feature< dim >::sptr  feature_sptr_type;
  typedef typename vcl_vector< feature_sptr_type >::size_type  feature_set_size_type;

  // compute mean distance between points, this will denote resolution of the data
  // do this using a kd tree to search for nearest neighbors
  // find bounding box of the data set
  //
  // build kd_tree from feature points
  vcl_vector< rsdl_point > points;
  rsdl_point  pt( dim, 0 );
  for( feature_set_size_type  j = 0; j < feature_points.size(); ++j ) {
    //pt.set_cartesian( feature_points[j]->location_ );
    for( unsigned int d = 0; d < dim; ++d ) pt.cartesian( d ) = feature_points[j]->location_[d];
    points.push_back( pt );
  } 
  rsdl_kd_tree  kd_tree( points );

  // compute mean distance between points and bounding box
  double mean_dist = 0.0;
  min_pt = feature_points[0]->location_;
  max_pt = feature_points[0]->location_;
  for( feature_set_size_type  j = 0; j < feature_points.size(); ++j ) {
    vnl_vector_fixed< typename cdcl_feature<dim>::coord_type, dim > curr_pt = feature_points[j]->location_;
    //pt.set_cartesian( curr_pt );
    for( unsigned int d = 0; d < dim; ++d ) pt.cartesian( d ) = curr_pt[d];

    // query kd-tree and find set of nearest points
    vcl_vector< rsdl_point >  points;
    vcl_vector< int >  indices;
    // use 2 nearest neighbors, first is the actual point (since the query point is from the tree)
    // we are interested in second (point closest to the query tree)
    int n = 2;
    kd_tree.n_nearest( pt, n, points, indices );

    double dist = rsdl_dist( pt, points[1] );//vcl_sqrt( vnl_vector_ssd( feature_points[j]->location_, feature_points[indices[0]] ) );
    mean_dist += dist;

    // search for the bounding box
    for( unsigned int d = 0; d < dim; ++d ) {
      if( curr_pt[d] > max_pt[d] ) max_pt[d] = curr_pt[d];
      if( curr_pt[d] < min_pt[d] ) min_pt[d] = curr_pt[d];
    }
  }
  mean_dist /= feature_points.size();

  return mean_dist;
}
                     

// Subsample data for a multiresolution computation.
template < unsigned int dim >
void subsample_data( vcl_vector< vbl_smart_ptr< cdcl_feature< dim > > >               const &  feature_points,
                     vcl_vector< vcl_vector< vbl_smart_ptr< cdcl_feature< dim > > > >       &  multires_points,
                     vcl_vector< typename cdcl_feature<dim>::coord_type >                   &  mean_spacing,
                     const unsigned int                                                        features_at_coarsest )
{
  typedef typename cdcl_feature< dim >::sptr  feature_sptr_type;
  typedef typename vcl_vector< feature_sptr_type >::size_type  feature_set_size_type;
  typedef typename cdcl_feature<dim>::coord_type  coord_type;

  vcl_vector< feature_sptr_type >  curr_points = feature_points;
  multires_points.clear();
  multires_points.push_back( curr_points );

  // coordinates of the bounding box of all points
  vnl_vector_fixed< coord_type, dim >  min_pt, max_pt;
  double mean_dist = data_spacing_and_bounding_box( feature_points, min_pt, max_pt );

  mean_spacing.clear();
  mean_spacing.push_back( mean_dist );

  vcl_cout << "Level " << multires_points.size()-1 << ": " << curr_points.size() << " points, spacing: " << mean_dist << vcl_endl;

  while( curr_points.size() > features_at_coarsest ) {

    // build kd_tree from feature points
    vcl_vector< rsdl_point > points;
    rsdl_point  pt( dim, 0 );
    for( feature_set_size_type  j = 0; j < curr_points.size(); ++j ) {
      //pt.set_cartesian( curr_points[j]->location_ );
      for( unsigned int d = 0; d < dim; ++d ) pt.cartesian( d ) = curr_points[j]->location_[d];
      points.push_back( pt );
    } 
    rsdl_kd_tree  kd_tree( points );

    // box size is equal to the resolution doubled
    double box_size = 2.0*mean_dist;

    vcl_vector< feature_sptr_type >  one_level;
    one_level.clear();

    // number of bins at each dimension
    vnl_vector_fixed< unsigned int, dim >  bin_count;
    // total number of bins
    unsigned int num_bins = 1;
    for( unsigned int d = 0; d < dim; ++d ) {
      bin_count[d] = static_cast<unsigned int>( vcl_ceil( ( max_pt[d] - min_pt[d] ) * ( 1.0/box_size ) ) );
      num_bins *= bin_count[d];
    }

    vnl_random  rand_gen( 9667566 );

    // slide bounding box of current resolution and randomly pick point in it
    //
    // current indices of the bounding box in the multiresolution grid (sliding box in a loop below)
    vnl_vector_fixed< unsigned int, dim >  c( 0u );
    for( unsigned int n = 1; n <= num_bins; ++n ) {

      // determine coordinates of the bounding box  
      rsdl_point  min_bb( dim ), max_bb( dim );
      for( unsigned int d = 0; d < dim; ++d ) {
        min_bb.cartesian(d) = min_pt[d] + ( c[d]   ) * box_size;
        max_bb.cartesian(d) = min_pt[d] + ( c[d]+1 ) * box_size;
      }
      //vcl_cout << min_bb << " " << max_bb << vcl_endl;

      // retrieve points in the current bounding box
      vcl_vector< rsdl_point >  points;
      vcl_vector< int >  indices;
      rsdl_bounding_box  box( min_bb, max_bb );
      kd_tree.points_in_bounding_box( box, points, indices );

      if( indices.size() > 0 ) {
        // select randomly point from the bounding box
        unsigned int rand_ind = rand_gen.lrand32( 0, indices.size()-1 );
        one_level.push_back( curr_points[indices[rand_ind]] );
      }

      // increment indices of the bounding box in the multiresolution grid
      unsigned int offset = n;
      for( unsigned int d = 0; d < dim; ++d ) {
        c[d] = offset % bin_count[d];
        offset = offset / bin_count[d];

        // for 3D
        //c[0] = n % offset1;
        //c[1] = ( n / offset1 ) % offset2;
        //c[2] = ( ( n / offset1 ) / offset2 ) % offset3;
      }

      //vcl_cout << c[0] << " " << c[1] << " " << c[2] << vcl_endl;
    }

    // recompute location covariance for the current level
    vcl_vector< unsigned int >  dummy;
    //location_covariance( one_level, box_size, dummy );







    multires_points.push_back( one_level );
    curr_points = one_level;

    // compute mean distance and boudning box for the current level
    mean_dist = data_spacing_and_bounding_box( curr_points, min_pt, max_pt );
    mean_spacing.push_back( mean_dist );

    vcl_cout << "Level " << multires_points.size()-1 << ": " << one_level.size() << " points, spacing: " << mean_dist << vcl_endl;
  }

}


// Subsample data for a multiresolution computation.
// Use covariances computed from the finest level (but different radius of the neighborhood).
template < unsigned int dim >
void subsample_data_fine_covariances( vcl_vector< vbl_smart_ptr< cdcl_feature< dim > > >               const &  feature_points,
                                      vcl_vector< vcl_vector< vbl_smart_ptr< cdcl_feature< dim > > > >       &  multires_points,
                                      vcl_vector< typename cdcl_feature<dim>::coord_type >                   &  mean_spacing )
{
  typedef typename cdcl_feature< dim >::sptr  feature_sptr_type;
  typedef typename vcl_vector< feature_sptr_type >::size_type  feature_set_size_type;
  typedef typename cdcl_feature<dim>::coord_type  coord_type;

  vcl_vector< feature_sptr_type >  curr_points = feature_points;
  multires_points.clear();
  multires_points.push_back( curr_points );

  // coordinates of the bounding box of all points
  vnl_vector_fixed< coord_type, dim >  min_pt, max_pt;
  double mean_dist = data_spacing_and_bounding_box( feature_points, min_pt, max_pt );

  mean_spacing.clear();
  mean_spacing.push_back( mean_dist );

  vcl_cout << "Level " << multires_points.size()-1 << ": " << curr_points.size() << " points, spacing: " << mean_dist << vcl_endl;

  const unsigned int features_at_coarsest = 300;
  // make a copy of the data which will be used over and over for computing covariances at each level
  vcl_vector< feature_sptr_type >  points_for_cov;
  if( curr_points.size() > features_at_coarsest ) {
    for( unsigned int i = 0; i < curr_points.size(); ++i ) {
      feature_sptr_type  feature_copy = new cdcl_feature<dim>( curr_points[i]->location_, curr_points[i]->covariance_ );
      points_for_cov.push_back( feature_copy );
    }
  }
  // vector of indices to the finest level
  // the first level is finest, so the vector contains all indices
  // as we subsample, this vector will contain less and less indices to the finest level
  vcl_vector< feature_set_size_type >  curr_indices;
  for( typename vcl_vector< feature_sptr_type >::size_type  t = 0; t < curr_points.size(); ++t ) curr_indices.push_back( t );
  while( curr_points.size() > features_at_coarsest ) {

    // box size is equal to the resolution doubled
    double box_size = 2.0*mean_dist;

    // source for covariance computation is always the finest level, but a different neighborhood size is used
    vcl_vector< feature_set_size_type >  noise_indices;
    location_covariance( points_for_cov, box_size, noise_indices );


    #if !RANDOM_COV_FOR_NOISE
    // remove noise indices from curr_indices and corresponding points in curr_points
    if( noise_indices.size() > 0 ) {
      for( typename vcl_vector< feature_set_size_type >::size_type  i = 0; i < noise_indices.size(); ++i ) {
        typename vcl_vector< feature_set_size_type >::iterator  found = vcl_find( curr_indices.begin(), curr_indices.end(), noise_indices[i] );
        if( found != curr_indices.end() ) {
          // this is assuming that curr_indices and curr_points have corresponding points at each array position
          unsigned int index = ( found - curr_indices.begin() ) / sizeof( feature_set_size_type );
          curr_indices.erase( found );
          curr_points.erase( curr_points.begin() + index );
        }
      }

      //typename vcl_vector< feature_sptr_type >  curr_points_wout_noise;
      //// copy only non-noise points
      //noise_indices.push_back( curr_points.size() ); // add the last index
      //typename vcl_vector< feature_set_size_type >::size_type  start_from = 0;
      //for( typename vcl_vector< feature_set_size_type >::size_type  i = 0; i < noise_indices.size(); ++i ) {
      //  for( feature_set_size_type j = start_from; j < noise_indices[i]; ++j ) {
      //    curr_points_wout_noise.push_back( curr_points[j] );
      //  }
      //  start_from = noise_indices[i]+1;
      //}
      //noise_indices.erase( noise_indices.begin() + noise_indices.size() - 1 ); // remove the last index

      //// assign non-noise points to the original feature vector
      //curr_points = curr_points_wout_noise;

    }
    #endif








    // build kd_tree from feature points
    vcl_vector< rsdl_point > points;
    rsdl_point  pt( dim, 0 );
    for( feature_set_size_type  j = 0; j < curr_points.size(); ++j ) {
      //pt.set_cartesian( curr_points[j]->location_ );
      for( unsigned int d = 0; d < dim; ++d ) pt.cartesian( d ) = curr_points[j]->location_[d];
      points.push_back( pt );
    } 
    rsdl_kd_tree  kd_tree( points );

    vcl_vector< feature_sptr_type >  one_level;
    one_level.clear();
    vcl_vector< feature_set_size_type >  one_level_indices_to_finest;
    one_level_indices_to_finest.clear();


    // number of bins at each dimension
    vnl_vector_fixed< unsigned int, dim >  bin_count;
    // total number of bins
    unsigned int num_bins = 1;
    for( unsigned int d = 0; d < dim; ++d ) {
      bin_count[d] = static_cast<unsigned int>( vcl_ceil( ( max_pt[d] - min_pt[d] ) * ( 1.0/box_size ) ) );
      num_bins *= bin_count[d];
    }

    vnl_random  rand_gen( 9667566 );

    // slide bounding box of current resolution and randomly pick point in it
    //
    // current indices of the bounding box in the multiresolution grid (sliding box in a loop below)
    vnl_vector_fixed< unsigned int, dim >  c( 0u );
    for( unsigned int n = 1; n <= num_bins; ++n ) {

      // determine coordinates of the bounding box  
      rsdl_point  min_bb( dim ), max_bb( dim );
      for( unsigned int d = 0; d < dim; ++d ) {
        min_bb.cartesian(d) = min_pt[d] + ( c[d]   ) * box_size;
        max_bb.cartesian(d) = min_pt[d] + ( c[d]+1 ) * box_size;
      }
      //vcl_cout << min_bb << " " << max_bb << vcl_endl;

      // retrieve points in the current bounding box
      vcl_vector< rsdl_point >  points;
      vcl_vector< int >  indices;
      rsdl_bounding_box  box( min_bb, max_bb );
      kd_tree.points_in_bounding_box( box, points, indices );

      if( indices.size() > 0 ) {

    //    // select randomly point from the bounding box
    //    unsigned int rand_ind = rand_gen.lrand32( 0, indices.size()-1 );
    //    unsigned int index_to_points = indices[rand_ind];
    //    // get the covariance from the finest level where we recomputed it (with different neighborhood size)
    //    unsigned int index_to_finest = curr_indices[index_to_points];

    //    typename vcl_vector< feature_set_size_type >::iterator  found = vcl_find( curr_indices.begin(), curr_indices.end(), noise_indices[i] );



    //    for( vcl_vector< int >::size_type  l = 0; l < indices.size(); ++l ) {
    //      for( typename vcl_vector< feature_set_size_type >::size_type  i = 0; i < noise_indices.size(); ++i ) {
    //        typename vcl_vector< feature_set_size_type >::iterator  found = vcl_find( curr_indices.begin(), curr_indices.end(), noise_indices[i] );
    //        if( found != curr_indices.end() ) {
    ////          unsigned int index = ( found - curr_indices.begin() ) / sizeof( feature_set_size_type );
    //          curr_indices.erase( found );
    ////          curr_points.erase( curr_points.begin() + index );
    //        }
    //    }


        // select randomly point from the bounding box
        unsigned int rand_ind = rand_gen.lrand32( 0, indices.size()-1 );
        unsigned int index_to_points = indices[rand_ind];
        feature_sptr_type  picked_point = curr_points[index_to_points];

        // get the covariance from the finest level where we recomputed it (with different neighborhood size)
        unsigned int index_to_finest = curr_indices[index_to_points];
        //if( picked_point->covariance_ != points_for_cov[index_to_finest]->covariance_ ) vcl_cout << "^^^^^^^^^^^^^^COVARIANCES DIFFERENT^^^^^^^^^^^^^^^" << vcl_endl;
        //vcl_cout << "OLD_COV: " << vcl_endl << picked_point->covariance_ << vcl_endl
        //         << "NEW_COV: " << vcl_endl << points_for_cov[index_to_finest]->covariance_ << vcl_endl;
        picked_point->covariance_ = points_for_cov[index_to_finest]->covariance_;
        one_level.push_back( picked_point );

        // save indices to the finest level (for next subsampled level)
        one_level_indices_to_finest.push_back( index_to_finest );
      }

      // increment indices of the bounding box in the multiresolution grid
      unsigned int offset = n;
      for( unsigned int d = 0; d < dim; ++d ) {
        c[d] = offset % bin_count[d];
        offset = offset / bin_count[d];

        // for 3D
        //c[0] = n % offset1;
        //c[1] = ( n / offset1 ) % offset2;
        //c[2] = ( ( n / offset1 ) / offset2 ) % offset3;
      }

      //vcl_cout << c[0] << " " << c[1] << " " << c[2] << vcl_endl;
    }

    multires_points.push_back( one_level );
    curr_points = one_level;
    curr_indices = one_level_indices_to_finest;

    // compute mean distance and bounding box for the current level
    mean_dist = data_spacing_and_bounding_box( curr_points, min_pt, max_pt );
    mean_spacing.push_back( mean_dist );

    vcl_cout << "Level " << multires_points.size()-1 << ": " << one_level.size() << " points, spacing: " << mean_dist << vcl_endl;
  }

}


// Save multiresolution data.
template < unsigned int dim >
void save_data( vcl_vector< vcl_vector< vbl_smart_ptr< cdcl_feature< dim > > > > const &  multires_points,
                vcl_string  prefix = "level" )
{
  typedef typename vcl_vector< vcl_vector< vbl_smart_ptr< cdcl_feature< dim > > > >::size_type  multires_points_size_type;
  typedef typename vcl_vector< vbl_smart_ptr< cdcl_feature< dim > > >::size_type                points_size_type;

  for( multires_points_size_type  n = 0; n < multires_points.size(); ++n ) {
    // prepare file name
    vcl_ostringstream  points_name;
    points_name << prefix << "_" << vcl_setw( 2 ) << vcl_setfill( '0' ) << n << ".dat";
    vcl_cout << "Saving " << points_name.str() << vcl_endl;

    vcl_ofstream  ostr( points_name.str().c_str() );
    for( points_size_type  i = 0; i < multires_points[n].size(); ++i ) {
      ostr << multires_points[n][i] << vcl_endl;
    }
    ostr.close();
  }

}

// Load multiresolution data return true if successful, false otherwise.
template < unsigned int dim >
bool load_data( vcl_vector< vcl_vector< vbl_smart_ptr< cdcl_feature< dim > > > > &  multires_points,
                vcl_string  prefix = "level" )
{
  typedef typename vcl_vector< vcl_vector< vbl_smart_ptr< cdcl_feature< dim > > > >::size_type  multires_points_size_type;
  typedef typename vcl_vector< vbl_smart_ptr< cdcl_feature< dim > > >::size_type                points_size_type;

  vcl_cout << "WARNING: Only one level loading implemented..." << vcl_endl;
  //for( multires_points_size_type  n = 0; n < multires_points.size(); ++n ) {
  //  // prepare file name
  //  vcl_ostringstream  points_name;
  //  points_name << prefix << vcl_setw( 2 ) << vcl_setfill( '0' ) << n << ".dat";
  //  vcl_cout << "Loading " << points_name.str() << vcl_endl;

  //  vcl_ifstream  istr( points_name.str().c_str() );
  //  vcl_vector< vbl_smart_ptr< cdcl_feature< dim > > >  features;
  //  while( !istr.eof() ) {
  //    vbl_smart_ptr< cdcl_feature< dim > >  feature;
  //    istr >> feature;
  //    features.push_back( feature );
  //  }
  //  istr.close();
  //}


  multires_points_size_type  n = 0;

  // prepare file name
  vcl_ostringstream  points_name;
  points_name << prefix << "_" << vcl_setw( 2 ) << vcl_setfill( '0' ) << n << ".dat";
  vcl_cout << "Loading " << points_name.str() << vcl_endl;

  vnl_matrix< typename cdcl_feature<dim>::coord_type >  raw_pts;

  // load ascii data
  vcl_ifstream  data_file( points_name.str().c_str() );
  bool have_read = raw_pts.read_ascii( data_file );
  data_file.close();

  if( !have_read ) {
    vcl_cerr << "Error: Can't read data from " << points_name.str() << vcl_endl;
    return false;
  }

  vcl_vector< vbl_smart_ptr< cdcl_feature< dim > > >  one_level;
  cdcl_parse_raw_data( raw_pts, one_level );
  multires_points.push_back( one_level );

  return true;
}


#define CDCL_UTILS_INSTANTIATE( dim )                                                                                              template                                                                                                                           int cdcl_parse_raw_data<dim>( vnl_matrix< cdcl_feature<dim>::coord_type >      const &  raw_matrix,                                                          vcl_vector< vbl_smart_ptr< cdcl_feature< dim > > >              &  feature_set );                         template                                                                                                                           void cdcl_radius_center<dim>( vcl_vector< vbl_smart_ptr< cdcl_feature< dim > > >  const &  feature_set,                                                               vnl_vector_fixed< cdcl_feature<dim>::coord_type, dim >    &  center,                                                                    cdcl_feature<dim>::coord_type                             &  avg_radius );                                template                                                                                                                           void cdcl_normalize_data<dim>( vcl_vector< vbl_smart_ptr< cdcl_feature< dim > > >     &  feature_set,                                                                      vnl_vector_fixed< cdcl_feature<dim>::coord_type, dim > &  center,                                                                           cdcl_feature<dim>::coord_type                          &  avg_radius );                                      template                                                                                                                           void cdcl_normalize_data_known_radius<dim>( vcl_vector< vbl_smart_ptr< cdcl_feature< dim > > >        &  feature_set,                                                               vnl_vector_fixed< cdcl_feature<dim>::coord_type, dim >    &  center,                                                                    cdcl_feature<dim>::coord_type                       const &  avg_radius );                  template                                                                                                                           void cdcl_unnormalize_data<dim>( vcl_vector< vbl_smart_ptr< cdcl_feature< dim > > >      const &  feature_set,                                                                vnl_vector_fixed< cdcl_feature<dim>::coord_type, dim >  const &  center,                                                                     cdcl_feature<dim>::coord_type                           const &  avg_radius );                              template                                                                                                                           void cdcl_normalize_matches<dim>( vcl_vector< vbl_smart_ptr< cdcl_match< dim > > >  & matches,                                                                            vnl_vector_fixed< cdcl_feature<dim>::coord_type, dim >                   & center_moving,                                                                      cdcl_feature<dim>::coord_type                                            & avg_radius_moving,                                                                  vnl_vector_fixed< cdcl_feature<dim>::coord_type, dim >                   & center_fixed,                                                                       cdcl_feature<dim>::coord_type                                            & avg_radius_fixed );                               template                                                                                                                           void cdcl_normalize_matches_known<dim>( vcl_vector< vbl_smart_ptr< cdcl_match< dim > > >        & matches,                                                                      vnl_vector_fixed< cdcl_feature<dim>::coord_type, dim >                   const & center_moving,                                                                cdcl_feature<dim>::coord_type                                            const & avg_radius_moving,                                                            vnl_vector_fixed< cdcl_feature<dim>::coord_type, dim >                   const & center_fixed,                                                                 cdcl_feature<dim>::coord_type                                            const & avg_radius_fixed );                                                                                                                                                      template                                                                                                                           void cdcl_unnormalize_matches<dim>( vcl_vector< vbl_smart_ptr< cdcl_match< dim > > >  & matches,                                                                       vnl_vector_fixed< cdcl_feature<dim>::coord_type, dim >             const & center_moving,                                                                      cdcl_feature<dim>::coord_type                                      const & avg_radius_moving,                                                                  vnl_vector_fixed< cdcl_feature<dim>::coord_type, dim >             const & center_fixed,                                                                       cdcl_feature<dim>::coord_type                                      const & avg_radius_fixed );                             template                                                                                                                           double data_spacing_and_bounding_box<dim>( vcl_vector< vbl_smart_ptr< cdcl_feature< dim > > >  const &  feature_points,                                                       vnl_vector_fixed< cdcl_feature<dim>::coord_type, dim >                &  min_pt,                                                   vnl_vector_fixed< cdcl_feature<dim>::coord_type, dim >                &  max_pt );      template                                                                                                                           void location_covariance<dim>( vcl_vector< vbl_smart_ptr< cdcl_feature< dim > > > &  feature_points,                                                              const cdcl_feature<dim>::coord_type &  spacing,                                                                                    vcl_vector< vcl_vector< vbl_smart_ptr< cdcl_feature< dim > > >::size_type > &  indices );           template                                                                                                                           void subsample_data<dim>( vcl_vector< vbl_smart_ptr< cdcl_feature< dim > > >               const &  feature_points,                                               vcl_vector< vcl_vector< vbl_smart_ptr< cdcl_feature< dim > > > >       &  multires_points,                                              vcl_vector< cdcl_feature<dim>::coord_type >                            &  mean_spacing,                                                 const unsigned int                                                        features_at_coarsest );             template                                                                                                                                void subsample_data_fine_covariances<dim>( vcl_vector< vbl_smart_ptr< cdcl_feature< dim > > >               const &  feature_points,                                               vcl_vector< vcl_vector< vbl_smart_ptr< cdcl_feature< dim > > > >       &  multires_points,                                              vcl_vector< cdcl_feature<dim>::coord_type >                            &  mean_spacing   );  template                                                                                                                                void save_data<dim>( vcl_vector< vcl_vector< vbl_smart_ptr< cdcl_feature< dim > > > > const &  multires_points,                                              vcl_string  prefix );                                                                                              template                                                                                                                                bool load_data<dim>( vcl_vector< vcl_vector< vbl_smart_ptr< cdcl_feature< dim > > > > &  multires_points,                                                    vcl_string  prefix );


#endif

