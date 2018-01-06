
#include "cdcl_trans_similarity2d.h"

#include <vnl/algo/vnl_svd.h>
#include <vnl/vnl_transpose.h>

//:
// \file
// \brief  2D similarity transformation.
// \author Michal Sofka
// \date   May 2006


// Default construct similarity transformation, default is identity.
cdcl_trans_similarity2d::cdcl_trans_similarity2d()
  : cdcl_trans< dim, dof >(),
    s_( 0.0 ),
    ang_( 0.0 )
{
  t_.fill( 0.0 );
  R_.fill( 0.0 );

  R_( 0, 0 ) = 1.0;
  R_( 1, 1 ) = 1.0;
}


// Construct similarity transformation given scale, angle, and translation.
cdcl_trans_similarity2d::cdcl_trans_similarity2d( coord_type                                                    const &  s,
                                                  coord_type                                                    const &  ang,
                                                  vnl_vector_fixed< coord_type, cdcl_trans_similarity2d::dim >  const &  t )
  : cdcl_trans< dim, dof >(),
    s_( s ),
    ang_( ang ),
    t_( t )
{
  coord_type sina = vcl_sin( ang_ );
  coord_type cosa = vcl_cos( ang_ );

  R_( 0, 0 ) =  cosa;
  R_( 0, 1 ) = -sina;
  R_( 1, 0 ) =  sina;
  R_( 1, 1 ) =  cosa;
}


// Copy construct similarity transformation.
cdcl_trans_similarity2d::cdcl_trans_similarity2d( cdcl_trans_similarity2d const &  similarity2d )
  : cdcl_trans< dim, dof >( similarity2d )
{
  s_   = similarity2d.s_;
  ang_ = similarity2d.ang_;
  R_   = similarity2d.R_;
  t_   = similarity2d.t_;
}


// Return inverse of the current transformation.
cdcl_trans< cdcl_trans_similarity2d::dim, cdcl_trans_similarity2d::dof >::sptr
cdcl_trans_similarity2d::inverse()
{
  vnl_svd< coord_type >  svdR( R_ );
  vnl_matrix< coord_type >  invR = svdR.inverse();

  // Angle.
  coord_type cosa = invR( 0, 0 );
  coord_type sina = invR( 1, 0 );

  // Rotation
  coord_type ang = vcl_atan2( sina, cosa );

  return new cdcl_trans_similarity2d( coord_type(1.0)/s_, ang, coord_type(-1.0)/s_ * invR * t_ );
}


// Convert transformation to parameterization suitable for optimizer.
vnl_vector< cdcl_trans_similarity2d::coord_type >
cdcl_trans_similarity2d::get_parameterization()
{
  coord_type sina = vcl_sin( ang_ );
  coord_type cosa = vcl_cos( ang_ );
//  coord_type sina = R_( 1, 0 );
//  coord_type cosa = R_( 0, 0 );
    
  // Populate vector with parameterization.
  vnl_vector< coord_type >  param( dof, 0.0 );
  param[0] = s_*cosa;
  param[1] = s_*sina;
  param[2] = t_[0];
  param[3] = t_[1];

  return param;
}


// Convert transformation from parameterization suitable for optimizer.
void 
cdcl_trans_similarity2d::set_parameterization( vnl_vector< coord_type > const &  param )
{
  assert( param.size() == dof );

  // Scale.
  s_ = vcl_sqrt( param[0]*param[0] + param[1]*param[1] );

  // Angle.
  coord_type cosa = param[0] / s_;
  coord_type sina = param[1] / s_;

  // Rotation
  ang_ = vcl_atan2( sina, cosa );
  R_( 0, 0 ) =  cosa;
  R_( 0, 1 ) = -sina;
  R_( 1, 0 ) =  sina;
  R_( 1, 1 ) =  cosa;

  // Translation.
  t_[0] = param[2];
  t_[1] = param[3];
}


// Map location.
vnl_vector_fixed< cdcl_trans_similarity2d::coord_type, cdcl_trans_similarity2d::dim >
cdcl_trans_similarity2d::map_loc( vnl_vector_fixed< coord_type, dim > const &  from ) const
{
  return s_ * R_ * from + t_;
}


// Compute Jacobian w.r.t. location.
vnl_matrix_fixed< cdcl_trans_similarity2d::coord_type, cdcl_trans_similarity2d::dim, cdcl_trans_similarity2d::dim >
cdcl_trans_similarity2d::jacobian_wrt_loc() const
{
  return s_ * R_;
}


// Compute Jacobian w.r.t. parameters.
vnl_matrix_fixed< cdcl_trans_similarity2d::coord_type, cdcl_trans_similarity2d::dim, cdcl_trans_similarity2d::dof > const &
cdcl_trans_similarity2d::jacobian_wrt_par( vnl_vector_fixed< coord_type, dim > const &  location )
{
  this->jacobian_wrt_par_thread( location, this->jacobian_wrt_par_ );

  return this->jacobian_wrt_par_;
}


// Compute Jacobian w.r.t. parameters.
void
cdcl_trans_similarity2d::jacobian_wrt_par_thread( vnl_vector_fixed< coord_type, dim > const &  location,
                                                  vnl_matrix_fixed< coord_type, dim, dof > &  jacobian_wrt_par )
{
  jacobian_wrt_par( 0, 0 ) = location[0];
  jacobian_wrt_par( 1, 0 ) = location[1];

  jacobian_wrt_par( 0, 1 ) = -location[1];
  jacobian_wrt_par( 1, 1 ) =  location[0];
    
  jacobian_wrt_par( 0, 2 ) = 1.0;
  jacobian_wrt_par( 1, 2 ) = 0.0;

  jacobian_wrt_par( 0, 3 ) = 0.0;
  jacobian_wrt_par( 1, 3 ) = 1.0;
}


// Compute Jacobian w.r.t. inverse parameters.
vnl_matrix_fixed< cdcl_trans_similarity2d::coord_type, cdcl_trans_similarity2d::dim, cdcl_trans_similarity2d::dof > const &
cdcl_trans_similarity2d::jacobian_wrt_inv_par( vnl_vector_fixed< coord_type, dim > const &  fixed_location )
{
  // The formula involves A^-1 (inverse of the forward rotation*scaling component).
  // This function is called on the inverse matrix, so A^-1 is in fact s_*R_.
  
  vnl_matrix_fixed< coord_type, dim, dim >  const A_1 = s_*R_;
  vnl_vector_fixed< coord_type, dim >  const  location_centered = fixed_location;
  vnl_matrix_fixed< coord_type, dim, dim >  Mkl( 0.0 );

  Mkl( 0, 0 ) = 1.0;
  Mkl( 1, 1 ) = 1.0;
  vnl_vector_fixed< coord_type, dim >  temp = -A_1 * Mkl * A_1 * location_centered;
  jacobian_wrt_inv_par_( 0, 0 ) = temp[0];
  jacobian_wrt_inv_par_( 1, 0 ) = temp[1];

  Mkl.fill( 0.0 );
  Mkl( 0, 1 ) = -1.0;
  Mkl( 1, 0 ) =  1.0;
  temp = - A_1 * Mkl * A_1 * location_centered;
  jacobian_wrt_inv_par_( 0, 1 ) = temp[0];
  jacobian_wrt_inv_par_( 1, 1 ) = temp[1];
    
  // -A_1 * [1;0]
  jacobian_wrt_inv_par_( 0, 2 ) = -A_1( 0, 0 );
  jacobian_wrt_inv_par_( 1, 2 ) = -A_1( 1, 0 );

  // -A_1 * [0;1]
  jacobian_wrt_inv_par_( 0, 3 ) = -A_1( 0, 1 );
  jacobian_wrt_inv_par_( 1, 3 ) = -A_1( 1, 1 );

  return jacobian_wrt_inv_par_;
}



// Compute Jacobian w.r.t. parameters of the Jacobian w.r.t. location. 
vcl_vector< vnl_matrix_fixed< cdcl_trans_similarity2d::coord_type, cdcl_trans_similarity2d::dim, cdcl_trans_similarity2d::dim > > const &
cdcl_trans_similarity2d::jacobian_of_Jp()
{
  jacobian_of_Jp_.clear();

  vnl_matrix_fixed< coord_type, dim, dim >  Mkl( 0.0 );
  Mkl( 0, 0 ) = 1.0;
  Mkl( 1, 1 ) = 1.0;
  jacobian_of_Jp_.push_back( Mkl );

  Mkl( 0, 0 ) =  0.0;
  Mkl( 0, 1 ) = -1.0;
  Mkl( 1, 0 ) =  1.0;
  Mkl( 1, 1 ) =  0.0;
  jacobian_of_Jp_.push_back( Mkl );

  return jacobian_of_Jp_;
}


// Convert transformation to normalized coordinates.
void
cdcl_trans_similarity2d::normalize( coord_type                          const &  avg_rad_moving,
                                    coord_type                          const &  avg_rad_fixed,
                                    vnl_vector_fixed< coord_type, dim > const &  center_moving,
                                    vnl_vector_fixed< coord_type, dim > const &  center_fixed )
{
  // convert scale
  s_ = avg_rad_moving / avg_rad_fixed * s_;

  // convert translation
  //vnl_matrix_fixed< coord_type, dim, dim >  sR = s_*R_;
  coord_type cosa = cos( ang_ );
  coord_type sina = sin( ang_ );
  vnl_matrix_fixed< coord_type, 2, 2 >  sR;
  sR( 0, 0 ) =  s_*cosa;
  sR( 0, 1 ) = -s_*sina;
  sR( 1, 0 ) =  s_*sina;
  sR( 1, 1 ) =  s_*cosa;

  t_ = t_ / avg_rad_fixed + sR * center_moving / avg_rad_moving - center_fixed / avg_rad_fixed;


  // normalization matrix
  vnl_matrix_fixed< coord_type, dof, dof >  B( 0.0 );

  coord_type temp = avg_rad_moving / avg_rad_fixed;
  B( 0, 0 ) = temp;
  B( 1, 1 ) = temp;

  const vnl_vector_fixed< coord_type, dim >  factor = center_moving / avg_rad_fixed;

  B( 2, 0 ) =  factor[0];
  B( 2, 1 ) = -factor[1];
  B( 3, 0 ) =  factor[1];
  B( 3, 1 ) =  factor[0];

  temp = 1.0 / avg_rad_fixed;
  B( 2, 2 ) = temp;
  B( 3, 3 ) = temp;

  this->covar_ = B * this->covar_ * B.transpose();
}


// Convert transformation back from normalized coordinates.
void
cdcl_trans_similarity2d::unnormalize( coord_type                          const &  avg_rad_moving,
                                      coord_type                          const &  avg_rad_fixed,
                                      vnl_vector_fixed< coord_type, dim > const &  center_moving,
                                      vnl_vector_fixed< coord_type, dim > const &  center_fixed )
{
  // convert scale
  s_ = avg_rad_fixed/avg_rad_moving * s_;

  // convert translation
  //vnl_matrix_fixed< coord_type, dim, dim >  sR = s_*R_;  
  coord_type cosa = cos( ang_ );
  coord_type sina = sin( ang_ );
  vnl_matrix_fixed< coord_type, 2, 2 >  sR;
  sR( 0, 0 ) =  s_*cosa;
  sR( 0, 1 ) = -s_*sina;
  sR( 1, 0 ) =  s_*sina;
  sR( 1, 1 ) =  s_*cosa;

  t_ = - sR * center_moving + avg_rad_fixed * t_ + center_fixed;


  // normalization matrix
  vnl_matrix_fixed< coord_type, dof, dof >  B( 0.0 );

  coord_type temp = avg_rad_fixed / avg_rad_moving;
  for( unsigned int d = 0; d < dim; ++d )
    B( d, d ) = temp;

  temp = avg_rad_fixed;
  for( unsigned int d = dim; d < dof; ++d )
    B( d, d ) = temp;

  this->covar_ = B * this->covar_ * B;  // don't need to transpose since B symmetrical
}


void
cdcl_trans_similarity2d::print( vcl_ostream &  ostr )
{
  ostr << "s: " << this->s_ << " ang: " << this->ang_ << " t: " << this->t_[0] << ", " << this->t_[1] << vcl_endl;
  ostr << "covar: " << vcl_endl << this->get_covariance();
  ostr << "covarJ: " << vcl_endl << this->get_covarianceJ();
}

