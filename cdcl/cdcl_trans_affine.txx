#ifndef cdcl_trans_affine_txx_
#define cdcl_trans_affine_txx_

#include "cdcl_trans_affine.h"

#include <vnl/algo/vnl_svd.h>
#include <vnl/vnl_transpose.h>

//:
// \file
// \brief  Affine transformation in N dmensions.
// \author Michal Sofka
// \date   Oct 2006


// Default construct affine transformation, default is identity.
template < unsigned int dm, unsigned int df >
cdcl_trans_affine<dm,df>::cdcl_trans_affine()
  : cdcl_trans< dm, df >(),
    A_( 0.0 ),
    t_( 0.0 )
{
  // initialize to identity matrix
  for( unsigned d = 0; d < dm; ++d ) A_( d, d ) = 1.0;
}


// Construct given affine matrix and translation.
template < unsigned int dm, unsigned int df >
cdcl_trans_affine<dm,df>::cdcl_trans_affine( vnl_matrix_fixed< coord_type, dm, dm > const &  A,
                                             vnl_vector_fixed< coord_type, dm >     const &  t,
                                             vnl_vector_fixed< coord_type, dm >     const &  center_moving )
  : cdcl_trans< dm, df >(),
    A_( A ),
    t_( t )
{
  this->center_moving_ = center_moving;
}


// Copy construct affine transformation 
template < unsigned int dm, unsigned int df >
cdcl_trans_affine<dm,df>::cdcl_trans_affine( cdcl_trans_affine<dm,df> const &  affine )
  : cdcl_trans< dm, df >( affine )
{
  this->A_   = affine.A_;
  this->t_   = affine.t_;
}


// Return inverse of the current transformation.
template < unsigned int dm, unsigned int df >
vbl_smart_ptr< cdcl_trans< dm, df > >
cdcl_trans_affine<dm,df>::inverse()
{
  vnl_svd< coord_type >  svdA( A_ );
  vnl_matrix< coord_type >  invA = svdA.inverse();
  return new cdcl_trans_affine<dm,df>( invA, -invA * t_ + this->center_moving_, vnl_vector_fixed< coord_type, dm >( 0.0 ) );
}


// Convert transformation to parameterization suitable for optimizer.
template < unsigned int dm, unsigned int df >
vnl_vector< typename cdcl_trans_affine<dm, df>::coord_type >
cdcl_trans_affine<dm,df>::get_parameterization()
{
  // Populate vector with parameterization.
  vnl_vector< coord_type >  param( df, 0.0 );
  unsigned int p = 0;

  // affine matrix
  for( unsigned int d1 = 0; d1 < dm; ++d1 )
    for( unsigned int d2 = 0; d2 < dm; ++d2 ) {
      param[p] = A_( d1, d2 );
      ++p;
    }

  // translation vector
  for( unsigned int d = 0; d < dm; ++d, ++p )
    param[p] = t_[d];

  return param;
}


// Convert transformation from parameterization suitable for optimizer.
template < unsigned int dm, unsigned int df >
void 
cdcl_trans_affine<dm,df>::set_parameterization( vnl_vector< coord_type > const &  param )
{
  assert( param.size() == df );

  // affine matrix
  unsigned int p = 0;
  for( unsigned int d1 = 0; d1 < dm; ++d1 )
    for( unsigned int d2 = 0; d2 < dm; ++d2 ) {
      A_( d1, d2 ) = param[p];
      ++p;
    }

  // translation vector
  for( unsigned int d = 0; d < dm; ++d, ++p ) {
    t_[d] = param[p];
  }

}


// Map location.
template < unsigned int dm, unsigned int df >
vnl_vector_fixed< typename cdcl_trans_affine<dm,df>::coord_type, dm >
cdcl_trans_affine<dm,df>::map_loc( vnl_vector_fixed< coord_type, dm > const &  from ) const
{
  return A_ * ( from - this->center_moving_ ) + t_;
}


// Compute Jacobian w.r.t. location.
template < unsigned int dm, unsigned int df >
vnl_matrix_fixed< typename cdcl_trans_affine<dm,df>::coord_type, dm, dm >
cdcl_trans_affine<dm,df>::jacobian_wrt_loc() const
{
  return A_;
}


// Compute Jacobian w.r.t. parameters.
template < unsigned int dm, unsigned int df >
vnl_matrix_fixed< typename cdcl_trans_affine<dm,df>::coord_type, dm, df > const &
cdcl_trans_affine<dm,df>::jacobian_wrt_par( vnl_vector_fixed< coord_type, dm > const &  location )
{
  this->jacobian_wrt_par_thread( location, this->jacobian_wrt_par_ );

  return this->jacobian_wrt_par_;
}


// Compute Jacobian w.r.t. parameters.
template < unsigned int dm, unsigned int df >
void
cdcl_trans_affine<dm,df>::jacobian_wrt_par_thread( vnl_vector_fixed< coord_type, dm > const &  location,
                                                   vnl_matrix_fixed< coord_type, dm, df >   &  jacobian_wrt_par )
{
  vnl_vector_fixed< coord_type, dm >  const  location_centered = location - this->center_moving_;

  jacobian_wrt_par.fill( 0.0 );

  // derivative w.r.t. each of the affine matrix parameters
  unsigned int p = 0;
  for( unsigned int d1 = 0; d1 < dm; ++d1 ) 
    for( unsigned int d2 = 0; d2 < dm; ++d2 ) {
      jacobian_wrt_par( d1, p ) = location_centered[d2];
      ++p;
    }

  // derivative w.r.t. each of the translation vector parameters
  for( unsigned int d = 0; d < dm; ++d, ++p ) {
    jacobian_wrt_par( d, p ) = 1.0;
  }

}


// Compute Jacobian w.r.t. inverse parameters.
template < unsigned int dm, unsigned int df >
vnl_matrix_fixed< typename cdcl_trans_affine<dm,df>::coord_type, dm, df > const &
cdcl_trans_affine<dm,df>::jacobian_wrt_inv_par( vnl_vector_fixed< coord_type, dm > const &  fixed_location )
{
  // The formula involves A^-1 (inverse of the forward affine component).
  // This function is called on the inverse matrix, so A^-1 is in fact A_.
  
  vnl_vector_fixed< coord_type, dm >  const  location_centered = fixed_location;
  vnl_matrix_fixed< coord_type, dm, dm >  Mkl( 0.0 );

  // derivative w.r.t. each of the affine matrix parameters
  for( unsigned int d1 = 0; d1 < dm; ++d1 )
    for( unsigned int d2 = 0; d2 < dm; ++d2 ) {
      Mkl( d1, d2 ) = 1.0;
      vnl_vector_fixed< coord_type, dm >  temp = -A_ * Mkl * A_ * location_centered;
      this->jacobian_wrt_inv_par_.set_column( d1*dm+d2, temp );
      Mkl( d1, d2 ) = 0.0;
    }

  // derivative w.r.t. each of the translation vector parameters
  // -A_ * [1;0],  -A_ * [0;1]
  // -A_ * [1;0;0],  -A_ * [0;1;0],  -A_ * [0;0;1]
  unsigned int p = dm * dm;
  for( unsigned int d = 0; d < dm; ++d, ++p ) {
    this->jacobian_wrt_inv_par_.set_column(  p, -A_.get_column( d ) );
  }

  return this->jacobian_wrt_inv_par_;
}


// Compute Jacobian w.r.t. parameters of the Jacobian w.r.t. location.
template < unsigned int dm, unsigned int df >
vcl_vector< vnl_matrix_fixed< typename cdcl_trans_affine<dm,df>::coord_type, dm, dm > > const &
cdcl_trans_affine<dm,df>::jacobian_of_Jp()
{
  this->jacobian_of_Jp_.clear();

  vnl_matrix_fixed< coord_type, dm, dm >  Mkl( 0.0 );

  // derivative w.r.t. each of the Jp term
  for( unsigned int d1 = 0; d1 < dm; ++d1 )
    for( unsigned int d2 = 0; d2 < dm; ++d2 ) {
      Mkl( d1, d2 ) = 1.0;
      this->jacobian_of_Jp_.push_back( Mkl );
      Mkl( d1, d2 ) = 0.0;
    }

  return this->jacobian_of_Jp_;
}


// Convert transformation to normalized coordinates.
template < unsigned int dm, unsigned int df >
void
cdcl_trans_affine<dm,df>::normalize( coord_type                         const &  avg_rad_moving,
                                     coord_type                         const &  avg_rad_fixed,
                                     vnl_vector_fixed< coord_type, dm > const &  center_moving,
                                     vnl_vector_fixed< coord_type, dm > const &  center_fixed )
{
  // convert affine matrix
  A_ = avg_rad_moving / avg_rad_fixed * A_;

  // convert translation
  t_ = t_ / avg_rad_fixed + A_ * ( center_moving - this->center_moving_ ) / avg_rad_moving - center_fixed / avg_rad_fixed;


  // normalization matrix
  vnl_matrix_fixed< coord_type, df, df >  B( 0.0 );

  coord_type temp = avg_rad_moving / avg_rad_fixed;
  unsigned int dm2 = dm*dm;
  for( unsigned int d = 0; d < dm2; ++d )
    B( d, d ) = temp;

  const vnl_vector_fixed< coord_type, dm >  factor = ( center_moving - this->center_moving_ ) / avg_rad_fixed;
  for( unsigned int d1 = 0; d1 < dm; ++d1 )
    for( unsigned int d2 = 0; d2 < dm; ++d2 ) {
      B( dm2+d1, d1*dm+d2 ) = factor[d2];
    }

  temp = 1.0 / avg_rad_fixed;
  for( unsigned int d = dm2; d < df; ++d )
    B( d, d ) = temp;

  this->covar_ = B * this->covar_ * B.transpose();

  this->center_moving_.fill( 0.0 );
}


// Convert transformation back from normalized coordinates.
template < unsigned int dm, unsigned int df >
void
cdcl_trans_affine<dm,df>::unnormalize( coord_type                         const &  avg_rad_moving,
                                       coord_type                         const &  avg_rad_fixed,
                                       vnl_vector_fixed< coord_type, dm > const &  center_moving,
                                       vnl_vector_fixed< coord_type, dm > const &  center_fixed )
{
  // convert affine matrix
  A_ = avg_rad_fixed/avg_rad_moving * A_;

  // convert translation
  t_ = avg_rad_fixed * t_ + center_fixed;


  // normalization matrix
  vnl_matrix_fixed< coord_type, df, df >  B( 0.0 );

  coord_type temp = avg_rad_fixed / avg_rad_moving;
  unsigned int dm2 = dm*dm;
  for( unsigned int d = 0; d < dm2; ++d )
    B( d, d ) = temp;

  temp = avg_rad_fixed;
  for( unsigned int d = dm2; d < df; ++d )
    B( d, d ) = temp;

  this->covar_ = B * this->covar_ * B;  // don't need to transpose since B symmetrical
  
  this->center_moving_ = center_moving;
}


template < unsigned int dm, unsigned int df >
void
cdcl_trans_affine<dm,df>::print( vcl_ostream &  ostr )
{
  ostr << "A: " << this->A_ << " t: " << this->t_ << vcl_endl;
  ostr << "from_center: " << this->center_moving_ << vcl_endl;
  ostr << "covar: " << vcl_endl << this->get_covariance() << vcl_endl;
  //ostr << "covarJ: " << vcl_endl << this->get_covarianceJ();
}

#define CDCL_TRANS_AFFINE_INSTANTIATE( dm, df ) \
template class cdcl_trans_affine< dm, df >;

#endif
