#ifndef cdcl_keypoint_h_
#define cdcl_keypoint_h_

#include "cdcl_feature.h"

//:
// \file
// \brief  Keypoint derived from the feature for registration with a location and a covariance.
//         Additional members: keypoint direction.
// \author Michal Sofka
// \date   Sep 2007

template < unsigned int dim >
class cdcl_keypoint : public cdcl_feature< dim >
{
public:
  typedef cdcl_feature< dim >   superclass;
  typedef cdcl_keypoint< dim >  self;
  typedef vbl_smart_ptr< self >  sptr;
  typedef typename superclass::coord_type  coord_type;

  // Default constructor.
  cdcl_keypoint()
    : superclass(),
      normal_( 0.0 ) {}

  // Copy constructor.
  cdcl_keypoint( cdcl_keypoint< dim >  const &  other )
    : superclass( other ),
      normal_( other.normal_ ) {}

  // Construct feature given a feature and a direction.
  cdcl_keypoint( superclass                            const &  feature,
                 vnl_vector_fixed< coord_type, dim >   const &  normal )
    : superclass( feature ),
      normal_( normal ) {}

  // Construct feature give location, strength and covariance.
  cdcl_keypoint( vnl_vector_fixed< coord_type, dim >      const &  location,
                 double                                   const &  strength,
                 vnl_matrix_fixed< coord_type, dim, dim > const &  covariance,
                 vnl_vector_fixed< coord_type, dim >      const &  normal )
    : normal_( normal ) {
      this->location_ = location;
      this->strength_ = strength;
      this->covariance_ = covariance;
  }


  // Dynamically create copy of the current feature object.
  virtual typename cdcl_feature< dim >::sptr
  clone() { return new cdcl_keypoint< dim >( *this ); }

//private:
public:
  // Keypoint normal.
  vnl_vector_fixed< coord_type, dim >       normal_;

};


// 3d version has a binormal
VCL_DEFINE_SPECIALIZATION
class cdcl_keypoint<3> : public cdcl_feature< 3 >
{
public:
  typedef cdcl_feature< 3 >   superclass;
  typedef cdcl_keypoint< 3 >  self;
  typedef vbl_smart_ptr< self >  sptr;
  typedef superclass::coord_type  coord_type;

  // Default constructor.
  cdcl_keypoint()
    : superclass(),
      normal_( 0.0 ),
      binormal_( 0.0 ) {}

  // Copy constructor.
  cdcl_keypoint( cdcl_keypoint< 3 >  const &  other )
    : superclass( other ),
      normal_( other.normal_ ),
      binormal_( other.binormal_ ) {}

  // Construct feature given a feature and a direction.
  cdcl_keypoint( superclass                          const &  feature,
                 vnl_vector_fixed< coord_type, 3 >   const &  normal,
                 vnl_vector_fixed< coord_type, 3 >   const &  binormal )
    : superclass( feature ),
      normal_( normal ),
      binormal_( binormal ) {}

  // Construct feature give location, strength and covariance.
  cdcl_keypoint( vnl_vector_fixed< coord_type, 3 >    const &  location,
                 double                               const &  strength,
                 vnl_matrix_fixed< coord_type, 3, 3 > const &  covariance,
                 vnl_vector_fixed< coord_type, 3 >    const &  normal,
                 vnl_vector_fixed< coord_type, 3 >    const &  binormal )
    : normal_( normal ),
      binormal_( binormal ) {
    this->location_ = location;
    this->strength_ = strength;
    this->covariance_ = covariance;
  }


  // Dynamically create copy of the current feature object.
  virtual cdcl_feature< 3 >::sptr
  clone() { return new cdcl_keypoint< 3 >( *this ); }

//private:
public:
  // Keypoint normal.
  vnl_vector_fixed< coord_type, 3 >       normal_;

  // Keypoint binormal.
  vnl_vector_fixed< coord_type, 3 >       binormal_;

};


template < unsigned int dim >
vcl_ostream&
operator<<( vcl_ostream& ostr, const vbl_smart_ptr< cdcl_keypoint<dim> > & feature_sptr )
{
  for( unsigned int d = 0; d < dim; ++d ) {
    ostr << feature_sptr->location_[d] << " ";
  }
  ostr << feature_sptr->normal_ << " ";
  for( unsigned int d1 = 0; d1 < dim; ++d1 ) {
    for( unsigned int d2 = d1; d2 < dim; ++d2 ) {
      ostr << feature_sptr->covariance_(d1, d2);
      if( d1 < dim-1 || d2 < dim-1 ) ostr << " ";
    }
  }
  return ostr;
}


//template < unsigned int dim >
//vcl_istream&
//operator>>( vcl_istream& istr, const vbl_smart_ptr< cdcl_keypoint<dim> > & feature_sptr )
//{
//  for( unsigned int d = 0; d < dim; ++d ) {
//    istr >> feature_sptr->location_[d] << " ";
//  }
//  for( unsigned int d1 = 0; d1 < dim; ++d1 ) {
//    for( unsigned int d2 = d1; d2 < dim; ++d2 ) {
//      istr >> feature_sptr->covariance_(d1, d2);
//      feature_sptr->covariance_(d2, d1) = feature_sptr->covariance_(d1, d2);
//    }
//  }
//  return istr;
//}



#endif
