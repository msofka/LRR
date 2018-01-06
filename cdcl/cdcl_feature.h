#ifndef cdcl_feature_h_
#define cdcl_feature_h_

#include <vcl_ostream.h>
#include <vbl/vbl_ref_count.h>
#include <vbl/vbl_smart_ptr.h>
#include <vnl/vnl_vector_fixed.h>
#include <vnl/vnl_matrix_fixed.h>

//:
// \file
// \brief  Feature for registration with a location and a covariance.
// \author Michal Sofka
// \date   May 2006

template < unsigned int dim >
class cdcl_feature : public vbl_ref_count
{
public:
  typedef vbl_smart_ptr< cdcl_feature< dim > >  sptr;
  typedef double coord_type;

  // Default constructor.
  cdcl_feature()
    : location_( 0.0 ),
      strength_( 0.0 ),
      covariance_( 0.0 ) {}

  // Copy constructor.
  cdcl_feature( cdcl_feature< dim >  const &  other )
    : location_( other.location_ ),
      strength_( other.strength_ ),
      covariance_( other.covariance_ ) {}

  // Construct feature given location, strength and covariance.
  cdcl_feature( vnl_vector_fixed< coord_type, dim >      const &  location,
                double                                   const &  strength,
                vnl_matrix_fixed< coord_type, dim, dim > const &  covariance )
    : location_( location ),
      strength_( strength ),
      covariance_( covariance ) {}

  // Construct feature given location and covariance.
  cdcl_feature( vnl_vector_fixed< coord_type, dim >      const &  location,
                vnl_matrix_fixed< coord_type, dim, dim > const &  covariance )
    : location_( location ),
      strength_( 0.0 ),
      covariance_( covariance ) {}

  // Dynamically create copy of the current feature object.
  virtual typename cdcl_feature< dim >::sptr
  clone() { return new cdcl_feature< dim >( *this ); }

//private:
public:
  // Feature location.
  vnl_vector_fixed< coord_type, dim >       location_;

  // Strength of a feature.
  double strength_;

  // Feature covariance.
  vnl_matrix_fixed< coord_type, dim, dim >  covariance_;
};


template < unsigned int dim >
vcl_ostream&
operator<<( vcl_ostream& ostr, const vbl_smart_ptr< cdcl_feature<dim> > & feature_sptr )
{
  for( unsigned int d = 0; d < dim; ++d ) {
    ostr << feature_sptr->location_[d] << " ";
  }
  for( unsigned int d1 = 0; d1 < dim; ++d1 ) {
    for( unsigned int d2 = d1; d2 < dim; ++d2 ) {
      ostr << feature_sptr->covariance_(d1, d2);
      if( d1 < dim-1 || d2 < dim-1 ) ostr << " ";
    }
  }
  return ostr;
}


template < unsigned int dim >
vcl_istream&
operator>>( vcl_istream& istr, const vbl_smart_ptr< cdcl_feature<dim> > & feature_sptr )
{
  for( unsigned int d = 0; d < dim; ++d ) {
    istr >> feature_sptr->location_[d] << " ";
  }
  for( unsigned int d1 = 0; d1 < dim; ++d1 ) {
    for( unsigned int d2 = d1; d2 < dim; ++d2 ) {
      istr >> feature_sptr->covariance_(d1, d2);
      feature_sptr->covariance_(d2, d1) = feature_sptr->covariance_(d1, d2);
    }
  }
  return istr;
}



#endif
