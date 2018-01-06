#ifndef cdcl_feature_ICP_h_
#define cdcl_feature_ICP_h_

#include <vcl_ostream.h>
#include <vcl_vector.h>
#include <vbl/vbl_ref_count.h>
#include <vbl/vbl_smart_ptr.h>
#include <vnl/vnl_vector_fixed.h>
#include <vnl/vnl_matrix_fixed.h>

//#include "cdcl_feature.h"

//:
// \file
// \brief  Feature for registration with a location, strength,
//         shape (corner, tube, sheet) and error projector (depending on the shape).
// \author Michal Sofka
// \date   May 2008

#define HAS_DIRECTIONS 1

template < unsigned int dim >
class cdcl_feature_ICP : public vbl_ref_count
{
public:
  typedef vbl_smart_ptr< cdcl_feature_ICP< dim > >  sptr;
  typedef double coord_type;
  typedef vnl_vector_fixed< coord_type, dim >  location_type;


  enum shape_type { CORNER, TUBE, SHEET };

  // Default constructor.
  cdcl_feature_ICP()
    : location_( 0.0 ),
      strength_( 0.0 ),
      shape_( CORNER ),
      error_projector_( 0.0 ) {}

  // Copy constructor.
  cdcl_feature_ICP( cdcl_feature_ICP< dim >  const &  other )
    : location_( other.location_ ),
      strength_( other.strength_ ),
      shape_( other.shape_ ),
      error_projector_( other.error_projector_ )
      #if HAS_DIRECTIONS
      , directions_( other.directions_ )
      #endif
      {}

  // Construct feature given location, strength and covariance.
  cdcl_feature_ICP( vnl_vector_fixed< coord_type, dim >      const &  location,
                    double                                   const &  strength,
                    vnl_matrix_fixed< coord_type, dim, dim > const &  error_projector,
                    shape_type                                        shape
                    #if HAS_DIRECTIONS
                    , vcl_vector< vnl_vector_fixed< coord_type, dim > > const &  directions
                    #endif
                    )
    : location_( location ),
      strength_( strength ),
      error_projector_( error_projector ),
      shape_( shape )
      #if HAS_DIRECTIONS
      , directions_( directions )
      #endif
      {}

  // Construct feature given location and covariance.
  cdcl_feature_ICP( vnl_vector_fixed< coord_type, dim >      const &  location )
    : location_( location ) {}

  // Dynamically create copy of the current feature object.
  virtual typename cdcl_feature_ICP< dim >::sptr
  clone() { return new cdcl_feature_ICP< dim >( *this ); }

//private:
public:
  // Feature location.
  location_type       location_;

  // Strength of a feature.
  float strength_;

  // Feature shape type (corner, sheet, or tube).  
  shape_type  shape_;

  // Error projector, in 2D, this is  I - tt^T, in 3D, it's  nn^T .
  // Residual is computed as: e^T  error_projector_  e .
  vnl_matrix_fixed< coord_type, dim, dim >  error_projector_;

  // Directions of the feature. Depending on the shape, the vector can hold one, two, or three directions.
  #if HAS_DIRECTIONS
  vcl_vector< vnl_vector_fixed< coord_type, dim > >  directions_;
  #endif
};


//// Define a postfix increment operator for feature shape.
//template < unsigned int dim >
//inline typename cdcl_feature_ICP<dim>::shape_type operator++( typename cdcl_feature_ICP<dim>::shape_type &rs, int ) {
//   cdcl_feature_ICP<dim>::shape_type  temp = rs;
//   ++rs;
//   return temp;
//}
//
//// Define a prefix increment operator for feature shape.
//template < unsigned int dim >
//inline typename cdcl_feature_ICP<dim>::shape_type operator++( typename cdcl_feature_ICP<dim>::shape_type &rs ) {
//   rs = (cdcl_feature_ICP<dim>::shape_type)(rs + 1);
//   return rs;
//}


template < unsigned int dim >
vcl_ostream&
operator<<( vcl_ostream& ostr, const vbl_smart_ptr< cdcl_feature_ICP<dim> > & feature_sptr )
{
  for( unsigned int d = 0; d < dim; ++d ) {
    ostr << feature_sptr->location_[d] << " ";
  }

  return ostr;
}


template < unsigned int dim >
vcl_istream&
operator>>( vcl_istream& istr, const vbl_smart_ptr< cdcl_feature_ICP<dim> > & feature_sptr )
{
  for( unsigned int d = 0; d < dim; ++d ) {
    istr >> feature_sptr->location_[d] << " ";
  }

  return istr;
}



#endif
