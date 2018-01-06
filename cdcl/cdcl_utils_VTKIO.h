#ifndef cdcl_utils_VTKIO_h_
#define cdcl_utils_VTKIO_h_

#include "cdcl_feature.h"
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
// \brief  Converting keypoints to/from poly_data.
// \author Michal Sofka
// \date   Sep 2007


template < unsigned int dim >
void
cdcl_keypoints_to_poly_data_VTK( vcl_vector< vbl_smart_ptr< cdcl_keypoint< dim > > >  const &  keypoints,
                                 vtkSmartPointer< vtkPolyData >                             &  poly_data );

template < unsigned int dim >
void
cdcl_poly_data_to_keypoints_VTK( vtkSmartPointer< vtkPolyData >                       const &  poly_data,
                                 vcl_vector< vbl_smart_ptr< cdcl_keypoint< dim > > >        &  keypoints );





#endif
