#ifndef cdcl_extract_data_h_
#define cdcl_extract_data_h_

#include <cdcl/cdcl_feature.h>
#include <cdcl/cdcl_match.h>
#include <cdcl/cdcl_trans.h>

#include <vnl/vnl_transpose.h>

#include "vtkPolyData.h"
#include "vtkSmartPointer.h"

#include "vtkVoxelModeller.h"

#include "itkImage.h"
#include "itkImageToVTKImageFilter.h"
#include "itkVTKImageToImageFilter.h"
#include "itkRGBPixel.h"

//:
// \file
// \brief  Displaying points and matches with a display callback function.
//         The dependence of the estimation classes on the display is avoided
//         by giving a choice of supplying callback or not.
//         The other advantage is that data members do not neeed to be exposed
//         in order to be displayed.
//         In case callback is not available, estimation can still run
//         without displying progress.
//         Previously, callback generated display (when dependent on vgui)
//         but now it only prepares values (poly data) to render with VTK.
// \author Michal Sofka
// \date   Aug 2007


// Callback which displays points is passed to the registration object.
// This way the registration object does not depend on vgui.
template < unsigned int dim, unsigned int dof, class feature_typeT = cdcl_feature< dim > >
class cdcl_extract_data_callback3d
{
public:
  cdcl_extract_data_callback3d();

  typedef feature_typeT                                   feature_type;
  typedef typename feature_type::sptr                     feature_sptr_type;
  typedef typename cdcl_match< dim, feature_type >::sptr  match_sptr_type;
  typedef typename cdcl_trans< dim, dof >::sptr   trans_sptr_type;

  typedef itk::RGBPixel< unsigned char > RGBPixelType;
  typedef itk::Image<signed short, 2>                         InputImageType;
  typedef itk::Image<unsigned char, 3>                        ImageType;
  typedef itk::Image<RGBPixelType, 2>                         OutputImageType;
  typedef itk::ImageToVTKImageFilter< InputImageType >        ToVTKAdaptorFilterType;
  typedef itk::VTKImageToImageFilter< OutputImageType >       FromVTKAdaptorFilterType;
  typedef FromVTKAdaptorFilterType::OutputImageType     OutputImageType;
  typedef FromVTKAdaptorFilterType::OutputImagePointer  OutputImagePointer;


  // Given moving points, transform them and initilialize PolyDataMovingROIMapped.
  void initialize_moving_mapped_features( vcl_vector< feature_sptr_type > const &  moving,
                                          trans_sptr_type                 const &  trans );
  
  // Display 2d points.
  virtual void display_points( vcl_vector< feature_sptr_type > const &  moving,
                               vcl_vector< feature_sptr_type > const &  fixed,
                               vcl_vector< match_sptr_type >   const &  matches,
                               trans_sptr_type                 const &  trans,
                               unsigned int                             iteration );

  // Callback which calls member display_points to draw points on tableau.
  // It is done this way because 1) we need to know the caller, 2) callback must be static.
  static void display_points( void * caller,   
                              vcl_vector< feature_sptr_type > const &  moving,
                              vcl_vector< feature_sptr_type > const &  fixed,
                              vcl_vector< match_sptr_type >   const &  matches,
                              trans_sptr_type                 const &  transform,
                              unsigned int                             iteration )
  {
    cdcl_extract_data_callback3d* mySelf = (cdcl_extract_data_callback3d*)  caller;
    mySelf->display_points( moving, fixed, matches, transform, iteration );
    mySelf->finalize_display( iteration );
  }

  // Finalize display callback tableau (post redraw and save).
  virtual void finalize_display( unsigned int  iteration );

  // Poly data for the covariances of the correspondences.
  vtkSmartPointer< vtkPolyData >  GetPolyDataCovariances() { return m_PolyDataCovariances; }

  // Poly data for the mapped moving points.
  vtkSmartPointer< vtkPolyData >  GetPolyDataMovingROIMapped() { return m_PolyDataMovingROIMapped; }

protected:
  // Transfer error covariance.
  inline virtual vnl_matrix_fixed< double, dim, dim >  transfer_covar( vnl_matrix_fixed< double, dim, dof > const &  Jth,
                                                                       trans_sptr_type                      const &  trans ) {
    return Jth * trans->get_covariance() * vnl_transpose( Jth );
  }

private:
  // Poly data for the covariances of the correspondences.
  vtkSmartPointer< vtkPolyData >  m_PolyDataCovariances;

  // Poly data for the mapped moving points.
  vtkSmartPointer< vtkPolyData >  m_PolyDataMovingROIMapped;
  
};


#endif
