#ifndef __ImageToImageFilter_h
#define __ImageToImageFilter_h


// \brief  Off screen rendering.
// \author Michal Sofka, portions (GenerateInputRequestedRegion and GenerateOutputInformation)
//                       taken from itkProjectionImageFilter.
// \date   Oct 2007


#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#pragma warning ( disable : 4996 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif


#include "itkImageToImageFilter.h"

#include "vtkFileOutputWindow.h"
#include "vtkPolyData.h"
#include "vtkLookupTable.h"
#include "vtkPolyDataMapper.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkInteractorStyleTrackballCamera.h"
#include "vtkInteractorStyleImage.h"
#include "vtkPolyDataReader.h"
#include "vtkXMLPolyDataReader.h"
#include "vtkSphereSource.h"
#include "vtkGenericGlyph3DFilter.h"
#include "vtkGlyph3D.h"
#include "vtkTensorGlyph.h"
//#include "vtkTensorGlyphScaled.h"
#include <cassert>
#include <iostream>
#include <string>
#include "vtkCylinderSource.h"
#include "vtkCamera.h"
#include "vtkPNGReader.h"
#include "vtkPNGWriter.h"

#include "vtkContourFilter.h"
#include "vtkCleanPolyData.h"
#include "vtkPolyDataNormals.h"
#include "vtkGeometryFilter.h"


#include "vtkImageData.h"
#include "vtkUnstructuredGrid.h"
#include "vtkSmartPointer.h"
#include "vtkOpenGLRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkOpenGLRenderWindow.h"
#include "vtkTensorGlyph.h"
#include "vtkWindowToImageFilter.h"
#include "vtkActor.h"
#include "vtkPolyDataMapper.h"
#include "vtkJPEGWriter.h"
#include "vtkPNGReader.h"
#include "vtkGlyphSource2D.h"
#include "vtkProperty.h"
#include "vtkGraphicsFactory.h"
#include "vtkImagingFactory.h"
#include "vtkImageActor.h"
#include "vtkWindowLevelLookupTable.h"
#include "vtkImageMapToColors.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkImageViewer.h"
#include "vtkImageViewer2.h"
//#include "vtkMesaRenderer.h"
//#include "vtkMesaRenderWindow.h"
//#include "vtkMesaActor.h"
//#include "vtkMesaPolyDataMapper.h" 
#include "vtkXMLPolyDataReader.h"
#include "vtkExtractPolyDataGeometry.h"
#include "vtkExtractGeometry.h"
#include "vtkBox.h"
#include "vtkClipPolyData.h"


#include "itkImageSeriesReader.h"
#include "itkImageSeriesWriter.h"
#include "itkDICOMImageIO2.h"
#include "itkDICOMSeriesFileNames.h"
#include "itkNumericSeriesFileNames.h"

#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"

#include "itkImage.h"
#include "itkRegionOfInterestImageFilter.h"
#include "itkImageToVTKImageFilter.h"
#include "itkVTKImageToImageFilter.h"


namespace itk {

template <class TInputImageType, class TOutputImageType>
class ITK_EXPORT ImageSlicesWithGeometryFilter :
    public ImageToImageFilter<TInputImageType, TOutputImageType>
{
public:

  typedef ImageSlicesWithGeometryFilter               Self;
  typedef ImageToImageFilter<TInputImageType,TOutputImageType> Superclass;
  typedef SmartPointer<Self>                        Pointer;
  typedef SmartPointer<const Self>                  ConstPointer;

  /** Method for creation through object factory */
  itkNewMacro(Self);

  /** Run-time type information */
  itkTypeMacro(ImageSlicesWithGeometryFilter, ImageToImageFilter);

  /** Display */
  void PrintSelf( std::ostream& os, Indent indent ) const;

  typedef typename TInputImageType::PixelType  InputPixelType;
  typedef typename TOutputImageType::PixelType  OutputPixelType;

  typedef typename itk::ImageRegion<3>  RegionOfInterestType;

  itkGetMacro( ROI, RegionOfInterestType );
  itkSetMacro( ROI, RegionOfInterestType );

  typedef vtkSmartPointer< vtkPolyData >   PolyDataPointer;

  itkGetMacro( PolyData, PolyDataPointer );
  itkSetMacro( PolyData, PolyDataPointer );

  itkGetMacro( EnlargeBy, double );
  itkSetMacro( EnlargeBy, double );

  // We cannot use ImageViewer2's SetSliceOrientation, because that already implies the rendering is set up
  enum
  {
    SLICE_ORIENTATION_YZ = 0,
    SLICE_ORIENTATION_XZ = 1,
    SLICE_ORIENTATION_XY = 2
  };
  itkGetMacro( SliceOrientation, int );
  itkSetMacro( SliceOrientation, int );

protected:

  ImageSlicesWithGeometryFilter();

protected:

  typedef itk::RegionOfInterestImageFilter< TInputImageType, TInputImageType >  RegionInterestType;

  typedef itk::ImageToVTKImageFilter< TInputImageType >        ToVTKAdaptorFilterType;

  typedef vtkSmartPointer< vtkExtractGeometry >  ExtractGeometryFilterPointer;

  typedef itk::VTKImageToImageFilter< TOutputImageType >              ToITKAdaptorFilterType;

  typedef vtkSmartPointer< vtkLookupTable >  LookupTablePointer;
  typedef vtkSmartPointer< vtkSphereSource >  SphereSourecePointer;
  typedef vtkSmartPointer< vtkGlyph3D >  Glyph3DPointer;
  typedef vtkSmartPointer< vtkTensorGlyph/*Scaled*/ >  TensorGlyphPointer;
  typedef vtkSmartPointer< vtkPolyDataMapper >  PolyDataMapperPointer;
  typedef vtkSmartPointer< vtkActor >  ActorPointer;
  typedef vtkSmartPointer< vtkImageViewer2 >  ImageViewer2Pointer;
  typedef vtkSmartPointer< vtkWindowToImageFilter >  WindowImageFilterPointer;

  void GenerateData();

  void GenerateOutputInformation();

  void GenerateInputRequestedRegion();

private:

  ImageSlicesWithGeometryFilter(Self&);   // intentionally not implemented
  void operator=(const Self&);          // intentionally not implemented

  typename RegionInterestType::Pointer     m_RegionOfInterestFilter;
  typename ToVTKAdaptorFilterType::Pointer  m_ToVTKAdaptorFilter;
  typename ToITKAdaptorFilterType::Pointer  m_ToITKAdaptorFilter;
  ExtractGeometryFilterPointer  m_ExtractGeometryFilter;
  PolyDataPointer  m_PolyData;
  LookupTablePointer  m_LookupTable;
  SphereSourecePointer  m_SphereSource;
  Glyph3DPointer  m_Glyph3D;
  TensorGlyphPointer  m_TensorGlyph;
  PolyDataMapperPointer  m_PolyDataMapper;
  ActorPointer  m_Actor;
  ImageViewer2Pointer  m_ImageViewer2;
  WindowImageFilterPointer  m_WindowImageFilter;

  RegionOfInterestType  m_ROI;

  int m_SliceOrientation;

  double m_EnlargeBy;
};

}


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkImageSlicesWithGeometryFilter.txx"
#endif

#endif
