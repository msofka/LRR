
#ifndef cdcl_utils_h_
#define cdcl_utils_h_

#include <vcl_string.h>

#include <vnl/vnl_math.h>

#include "cdcl_feature.h"
#include "cdcl_match.h"
#include "cdcl_trans.h"

#include <vnl/algo/vnl_symmetric_eigensystem.h>

//:
// \file
// \brief  Utility functions for cdcl class.
// \author Michal Sofka
// \date   May 2006


// Note: It would be better if some of these types are imported from the respective classes (e.g. sptr)
//       but gcc complains then.

// Parse raw matrix data into a feature set.
template < unsigned int dim >
int cdcl_parse_raw_data( vnl_matrix< typename cdcl_feature<dim>::coord_type >        const &  raw_matrix, 
                         vcl_vector< vbl_smart_ptr< cdcl_feature< dim > > >       &  feature_set );

// Compute center and average radius of a feature set.
template < unsigned int dim >
void cdcl_radius_center( vcl_vector< vbl_smart_ptr< cdcl_feature< dim > > >      const &  feature_set,
                         vnl_vector_fixed< typename cdcl_feature<dim>::coord_type, dim >        &  center,
                         typename cdcl_feature<dim>::coord_type                                 &  avg_radius );

// Normalize data in a feature set to have location with zero mean and avg. radius of one.
template < unsigned int dim >
void cdcl_normalize_data( vcl_vector< vbl_smart_ptr< cdcl_feature< dim > > >           &  feature_set,
                          vnl_vector_fixed< typename cdcl_feature<dim>::coord_type, dim >       &  center,
                          typename cdcl_feature<dim>::coord_type                                &  avg_radius );

// Normalize data in a feature set to have location with zero mean and avg. radius of one.
// Use given radius.
template < unsigned int dim >
void cdcl_normalize_data_known_radius( vcl_vector< vbl_smart_ptr< cdcl_feature< dim > > >           &  feature_set,
                                       vnl_vector_fixed< typename cdcl_feature<dim>::coord_type, dim >       &  center,
                                       typename cdcl_feature<dim>::coord_type                          const &  avg_radius );

// Normalize data in a feature set to have location with zero mean and avg. radius of one.
// Use given center and radius.
template < unsigned int dim >
void cdcl_normalize_data_known( vcl_vector< vbl_smart_ptr< cdcl_feature< dim > > >            &  feature_set,
                                vnl_vector_fixed< typename cdcl_feature<dim>::coord_type, dim >  const &  center,
                                typename cdcl_feature<dim>::coord_type                           const &  avg_radius );

// Unnormalize data in a feature set.
template < unsigned int dim >
void cdcl_unnormalize_data( vcl_vector< vbl_smart_ptr< cdcl_feature< dim > > >       const &  feature_set,
                            vnl_vector_fixed< typename cdcl_feature<dim>::coord_type, dim >   const &  center,
                            typename cdcl_feature<dim>::coord_type                            const &  avg_radius );

// Normalize moving and fixed points which the matches were formed with.
template < unsigned int dim >
void cdcl_normalize_matches( vcl_vector< vbl_smart_ptr< cdcl_match< dim > > >        & matches, 
                             vnl_vector_fixed< typename cdcl_feature<dim>::coord_type, dim >  & center_moving, 
                             typename cdcl_feature<dim>::coord_type                           & avg_radius_moving, 
                             vnl_vector_fixed< typename cdcl_feature<dim>::coord_type, dim >  & center_fixed, 
                             typename cdcl_feature<dim>::coord_type                           & avg_radius_fixed );

// Normalize moving and fixed points which the matches were formed with
// Use supplied known centers and radii
template < unsigned int dim >
void cdcl_normalize_matches_known( vcl_vector< vbl_smart_ptr< cdcl_match< dim > > >        & matches, 
                                   vnl_vector_fixed< typename cdcl_feature<dim>::coord_type, dim >                   const & center_moving, 
                                   typename cdcl_feature<dim>::coord_type                                            const & avg_radius_moving, 
                                   vnl_vector_fixed< typename cdcl_feature<dim>::coord_type, dim >                   const & center_fixed, 
                                   typename cdcl_feature<dim>::coord_type                                            const & avg_radius_fixed );

// Unormalize moving and fixed points which the matches were formed with.
template < unsigned int dim >
void cdcl_unnormalize_matches( vcl_vector< vbl_smart_ptr< cdcl_match< dim > > >        & matches, 
                               vnl_vector_fixed< typename cdcl_feature<dim>::coord_type, dim >                   const & center_moving, 
                               typename cdcl_feature<dim>::coord_type                                            const & avg_radius_moving, 
                               vnl_vector_fixed< typename cdcl_feature<dim>::coord_type, dim >                   const & center_fixed, 
                               typename cdcl_feature<dim>::coord_type                                            const & avg_radius_fixed );

// Compute location covariance given a set of points.
template < unsigned int dim >
void location_covariance( vcl_vector< vbl_smart_ptr< cdcl_feature< dim > > >                                   &  feature_points,
                          const typename cdcl_feature<dim>::coord_type                                                                         &  spacing,
                          vcl_vector< typename vcl_vector< vbl_smart_ptr< cdcl_feature< dim > > >::size_type > &  indices );

template < unsigned int dim >
double data_spacing_and_bounding_box( vcl_vector< vbl_smart_ptr< cdcl_feature< dim > > >  const &  feature_points,
                                      vnl_vector_fixed< typename cdcl_feature<dim>::coord_type, dim >                           &  min_pt,
                                      vnl_vector_fixed< typename cdcl_feature<dim>::coord_type, dim >                           &  max_pt );

// Subsample data for a multiresolution computation.
// Return mean distance between points in the original dataset before subsampling.
template < unsigned int dim >
void subsample_data( vcl_vector< vbl_smart_ptr< cdcl_feature< dim > > >               const &  feature_points,
                     vcl_vector< vcl_vector< vbl_smart_ptr< cdcl_feature< dim > > > >       &  multires_points,
                     vcl_vector< typename cdcl_feature<dim>::coord_type >                            &  mean_spacing,
                     const unsigned int                                                        features_at_coarsest );
 
// Subsample data for a multiresolution computation.
// Use covariances computed from the finest level (but different radius of the neighborhood).
// Argument indices 
template < unsigned int dim >
void subsample_data_fine_covariances( vcl_vector< vbl_smart_ptr< cdcl_feature< dim > > >               const &  feature_points,
                                      vcl_vector< vcl_vector< vbl_smart_ptr< cdcl_feature< dim > > > >       &  multires_points,
                                      vcl_vector< typename cdcl_feature<dim>::coord_type >                   &  mean_spacing );

template < unsigned int dim >
void save_data( vcl_vector< vcl_vector< vbl_smart_ptr< cdcl_feature< dim > > > > const &  multires_points,
                vcl_string  prefix );

// Load multiresolution data return true if successful, false otherwise.
template < unsigned int dim >
bool load_data( vcl_vector< vcl_vector< vbl_smart_ptr< cdcl_feature< dim > > > > &  multires_points,
                vcl_string  prefix );

namespace {

// Beaton-Tukey robust weight
// Note: residual is assumed to have been squared already
// b ... Beaton-Tukey constant
//inline double weight_BT( double const &  residual )
//{
//  const double b = 4.0;
//  const double b_2 = b*b;
//
//  // Beaton-Tukey robust weight, typical b = 4 .. 5
//  if( residual <= b_2 ) {
//  	double temp = 1-residual/b_2;
//    return temp * temp;
//  }
//  else {
//    return 0;
//  }
//}
inline double weight_BT( double const &  residual )
{
  // Beaton-Tukey robust weight, typical b = 4 .. 5
  if( residual <= 16.0 ) {
  	double temp = 1-residual*0.0625;
    return temp * temp;
  }
  else {
    return 0;
  }
}


// Note: residual is assumed to have been squared already
// derivative of rhod
// rhod ... Modified Beaton-Tukey robust function (multiplied B-T by a constant)
inline double rho( double const &  residual )
{
  const double b = 4.0;
  const double b_2 = b*b;
  //const double f = 2.3931;  // k = 2.0
  //const double f = 2.6466;  // k = 2.5
  const double f = 2.9872;  // k = 3.0

  // Beaton-Tukey robust weight, typical b = 4 .. 5
  const double c = b_2 / 6.0;
  if( residual > b_2 ) {
    return f * c;
  }
  else {
    return f * c * ( 1 - vcl_pow( 1 - residual/b_2, 3.0 ) );
  }

  //// Geman-McClure
  //const double b_2 = 1.0;
  //return residual / ( residual + b_2 );
}


// Note: residual is assumed to have been squared already
// derivative of rhod
// rhod ... Modified Beaton-Tukey robust function (multiplied B-T by a constant)
inline double rho_p( double const &  residual )
{
  const double b = 4.0;
  const double b_2 = b*b;
  //const double f = 2.3931;  // k = 2.0
  //const double f = 2.6466;  // k = 2.5
  const double f = 2.9872;  // k = 3.0

  // Beaton-Tukey robust weight, typical b = 4 .. 5
  if( residual > b_2 ) {
    return 0.0;
  }
  else {
    double temp = 1 - residual/b_2;
    return f * 0.5 * temp * temp;
  }

  //// Geman-McClure
  //const double b_2 = 1.0;
  //return 2.0 / ( 2.0 *residual + b_2*vcl_sqrt( residual ) );
}


// Note: residual is assumed to have been squared already
// b ... Beaton-Tukey constant
inline double weight( double const &  residual )
{
  const double b = 4.0;
  const double b_2 = b*b;
  //const double f = 2.3931;  // k = 2.0
  //const double f = 2.6466;  // k = 2.5
  const double f = 2.9872;  // k = 3.0

  // Beaton-Tukey robust weight, typical b = 4 .. 5
  if( residual <= b_2 ) {
    return f * ( 1-residual/b_2 )*( 1-residual/b_2 );
  }
  else {
    return 0;
  }


  //// Geman-McClure
  //const double b_2 = 1.0;
  //return 2.0 / ( 2.0*residual + b_2*vcl_sqrt( residual ) );
}


inline double weight_Cauchy( double const & residual )
{
  const double b_2 = 9.0;
  //const double b_2 = 16.0;
  return 1.0 / ( 1.0 + residual / b_2 );
}

  
template < unsigned int dim >
double cdcl_gauss_weight( vnl_vector_fixed< typename cdcl_feature<dim>::coord_type, dim > const &  point, double const &  sigma )
{
  const double k = vcl_pow( 2*vnl_math::pi, dim/2.0 ) * vcl_sqrt( sigma );
  return 1.0 / k * vcl_exp( -0.5 * dot_product( point, point ) / (sigma*sigma) );
}


}


namespace {
inline double cdcl_largest_evalue( vnl_matrix_fixed< cdcl_feature<2>::coord_type, 2, 2 > const &  C )
{
  double temp = (( C(0,0) - C(1,1) )*( C(0,0) - C(1,1) ) + 4*C(0,1)*C(0,1) );
  double sqrt_temp = vcl_sqrt( temp );

  double Lpp = (C(0,0) + C(1,1) - sqrt_temp)/2.0;
  double Lqq = (C(0,0) + C(1,1) + sqrt_temp)/2.0;

  // sort the eigenvalues: |lambda1| < |lambda2|
  // return larger
  if( vcl_abs(Lpp) < vcl_abs(Lqq) )
    return Lqq;
  else
    return Lpp;
}

inline double cdcl_smallest_evalue( vnl_matrix_fixed< cdcl_feature<2>::coord_type, 2, 2 > const &  C )
{
  double temp = (( C(0,0) - C(1,1) )*( C(0,0) - C(1,1) ) + 4*C(0,1)*C(0,1) );
  double sqrt_temp = vcl_sqrt( temp );

  double Lpp = (C(0,0) + C(1,1) - sqrt_temp)/2.0;
  double Lqq = (C(0,0) + C(1,1) + sqrt_temp)/2.0;

  // sort the eigenvalues: |lambda1| < |lambda2|
  // return larger
  if( vcl_abs(Lpp) < vcl_abs(Lqq) )
    return Lpp;
  else
    return Lqq;
}

inline double cdcl_largest_evalue( vnl_matrix_fixed< cdcl_feature<3>::coord_type, 3, 3 > const &  C )
{
  double l1, l2, l3;
  vnl_symmetric_eigensystem_compute_eigenvals( C( 0, 0 ), C( 0, 1 ), C( 0, 2 ),
                                                          C( 1, 1 ), C( 1, 2 ),
                                                                     C( 2, 2 ),
                                               l1, l2, l3 );
  return l3;
}


inline double cdcl_smallest_evalue( vnl_matrix_fixed< cdcl_feature<3>::coord_type, 3, 3 > const &  C )
{
  double l1, l2, l3;
  vnl_symmetric_eigensystem_compute_eigenvals( C( 0, 0 ), C( 0, 1 ), C( 0, 2 ),
                                                          C( 1, 1 ), C( 1, 2 ),
                                                                     C( 2, 2 ),
                                               l1, l2, l3 );
  return l1;
}

}


#endif


