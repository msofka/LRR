#ifndef cdcl_utils_VTK_h_
#define cdcl_utils_VTK_h_

#include "cdcl_feature.h"
#include "cdcl_feature_with_shape.h"
#include "cdcl_feature_ICP.h"
#include "cdcl_keypoint.h"
#include "cdcl_match.h"
#include "cdcl_trans.h"

#include <vcl_string.h>

#include "vtkPolyData.h"
#include "vtkSmartPointer.h"



#include "vtkPointData.h"
#include "vtkDoubleArray.h"

//:
// \file
// \brief  Writing matches, points, and covariances (axes of the ellipse or ellipsoid).
// \author Michal Sofka
// \date   May 2006

// Write matches in a VTK format
template < unsigned int dim, unsigned int dof >
void
cdcl_write_matches_VTK( vcl_vector< vbl_smart_ptr< cdcl_feature< dim > > > const &  moving,
                        vcl_vector< vbl_smart_ptr< cdcl_feature< dim > > > const &  fixed,
                        vcl_vector< vbl_smart_ptr< cdcl_match< dim > > >   const &  matches,
                        vbl_smart_ptr< cdcl_trans< dim, dof > >            const &  trans,
                        unsigned int                                       const    iteration );

// THE FOLLOWING 3 SHOULD BE DONE IN TERMS OF ITERATOR (templated over iterator), then only one function would be enough (but with template specialization)
//
// Convert cdcl features to VTK polydata for display purposes
template < unsigned int dim >
void
cdcl_features_to_poly_data( vcl_vector< vbl_smart_ptr< cdcl_feature< dim > > > const &  points,
                            vtkSmartPointer< vtkPolyData >                           &  poly_data );

// Convert cdcl features to VTK polydata for display purposes
template < unsigned int dim >
void
cdcl_features_to_poly_data( vcl_vector< vbl_smart_ptr< cdcl_feature_with_shape< dim > > > const &  points,
                            vtkSmartPointer< vtkPolyData >                                      &  poly_data );

template < unsigned int dim >
void
cdcl_features_to_poly_data( vcl_vector< vbl_smart_ptr< cdcl_feature_ICP< dim > > > const &  points,
                            vtkSmartPointer< vtkPolyData >                               &  poly_data );

template < unsigned int dim >
bool
cdcl_read_keypoint_descriptors_VTK( vcl_vector< vbl_smart_ptr< cdcl_keypoint< dim > > >        &  keypoints,
                                    vcl_vector< vnl_vector< float > >                          &  descriptors,
                                    vcl_string                                           const    filename );

template < unsigned int dim >
void
cdcl_write_keypoint_descriptors_VTK( vcl_vector< vbl_smart_ptr< cdcl_keypoint< dim > > >  const &  keypoints,
                                     vcl_vector< vnl_vector< float > >                    const &  descriptors,
                                     vcl_string                                           const    filename );



#endif
