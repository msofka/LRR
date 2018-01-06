#ifndef __ImageToImageFilter_txx
#define __ImageToImageFilter_txx

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

#include "itkImageSlicesWithGeometryFilter.h"

namespace itk 
{

//  Software Guide : BeginCodeSnippet
template <class TInputImageType, class TOutputImageType>
ImageSlicesWithGeometryFilter<TInputImageType, TOutputImageType>
::ImageSlicesWithGeometryFilter()
{
  m_EnlargeBy = 10.0;



  m_SliceOrientation = SLICE_ORIENTATION_XY;

  m_RegionOfInterestFilter = RegionInterestType::New();
  m_ToVTKAdaptorFilter = ToVTKAdaptorFilterType::New();
  m_ToITKAdaptorFilter = ToITKAdaptorFilterType::New();
  m_ExtractGeometryFilter = ExtractGeometryFilterPointer::New();
  m_WindowImageFilter = WindowImageFilterPointer::New();

  // initialize ROI
  RegionOfInterestType::SizeType  size;
  RegionOfInterestType::IndexType  index;

  index[0] = 0;
  index[1] = 0;
  index[2] = 0;

  size[0] = 0;
  size[1] = 0;
  size[2] = 0;

  m_ROI.SetSize( size );
  m_ROI.SetIndex( index );

  m_RegionOfInterestFilter->SetRegionOfInterest( m_ROI );

  m_ToVTKAdaptorFilter->SetInput( m_RegionOfInterestFilter->GetOutput() );

  m_ExtractGeometryFilter->ExtractInsideOn();


  m_LookupTable = LookupTablePointer::New();
  double scalarRange[2] = { 0.0, 1.0 };
  m_LookupTable->SetTableRange( scalarRange );
  m_LookupTable->SetHueRange( 0.667, 0.0 );
  m_LookupTable->SetSaturationRange( 1, 1 );
  m_LookupTable->SetValueRange( 1, 1 );

  m_SphereSource = SphereSourecePointer::New();
  m_SphereSource->SetRadius( 5.0 );
  m_SphereSource->SetThetaResolution( 12 );
  m_SphereSource->SetPhiResolution( 12 );

  m_PolyDataMapper = PolyDataMapperPointer::New();
  
#define RENDER_COVARIANCES 0

#if !RENDER_COVARIANCES
  m_Glyph3D = Glyph3DPointer::New();
  m_Glyph3D->SetInput( m_ExtractGeometryFilter->GetOutput() );
  m_Glyph3D->SetSource( m_SphereSource->GetOutput() );
  m_Glyph3D->SetVectorModeToUseNormal();
  m_Glyph3D->SetScaleModeToScaleByVector();
  m_Glyph3D->SetScaleFactor( 0.1 );
  
  m_PolyDataMapper->SetInput( m_Glyph3D->GetOutput() );
#else
  vtkSmartPointer< vtkGlyphSource2D >  DashSource = vtkSmartPointer< vtkGlyphSource2D >::New();
  DashSource->SetGlyphTypeToDash();
  DashSource->SetScale( 30.0 );
  DashSource->SetScale2( 30.0 );
  DashSource->FilledOff();

  // warning, this is not TensorGlyphScaled, but only TensorGlyph
  m_TensorGlyph = TensorGlyphPointer::New();
  m_TensorGlyph->SetInput( m_ExtractGeometryFilter->GetOutput() );
  //m_TensorGlyph->SetSource( m_SphereSource->GetOutput() );
  m_TensorGlyph->SetSource( DashSource->GetOutput() );
  m_TensorGlyph->ExtractEigenvaluesOn();  // eigen values extracted and used
  m_TensorGlyph->ThreeGlyphsOn();  // not set -> one glyph produced and eigenvalues scale the glyph
  m_TensorGlyph->SymmetricOn(); // set -> each glyph is mirrored
  m_TensorGlyph->ScalingOn(); // set -> use additional scaling by a factor
  m_TensorGlyph->SetScaleFactor( 0.1 );
  m_TensorGlyph->ClampScalingOn();
  m_TensorGlyph->ColorGlyphsOff();
  
  m_PolyDataMapper->SetInput( m_TensorGlyph->GetOutput() );
#endif
  
  m_PolyDataMapper->SetLookupTable( m_LookupTable );
  m_PolyDataMapper->SetScalarRange( scalarRange );

  m_Actor = ActorPointer::New();
  m_Actor->SetMapper( m_PolyDataMapper );
  m_Actor->GetProperty()->SetLineWidth( 3 );
  m_Actor->GetProperty()->SetPointSize( 3 );
  m_Actor->GetProperty()->SetColor( 0.0, 1.0, 0.0 );


  m_ImageViewer2 = ImageViewer2Pointer::New();
  vtkSmartPointer< vtkRenderWindowInteractor >  renderWindowInteractor = vtkSmartPointer< vtkRenderWindowInteractor >::New();
  m_ImageViewer2->SetupInteractor( renderWindowInteractor );
  m_ImageViewer2->SetInput( m_ToVTKAdaptorFilter->GetOutput() );
  m_ImageViewer2->SetColorWindow( 255 );
  m_ImageViewer2->SetColorLevel( 128 );


}


template <class TInputImage, class TOutputImage>
void
ImageSlicesWithGeometryFilter<TInputImage,TOutputImage>
::GenerateOutputInformation()
{
  //Superclass::GenerateOutputInformation();


  typename TOutputImage::RegionType outputRegion;
  typename TInputImage::RegionType inputRegion;
  typename TInputImage::IndexType inputIndex;
  typename TInputImage::SizeType  inputSize;
  typename TOutputImage::SizeType  outputSize;
  typename TOutputImage::IndexType outputIndex;
  typename TInputImage::SpacingType inSpacing;
  typename TInputImage::PointType inOrigin;
  typename TOutputImage::SpacingType outSpacing;
  typename TOutputImage::PointType outOrigin;


  // the other two dimensions
  int d1 = (m_SliceOrientation + 1 ) % 3;
  int d2 = (m_SliceOrientation + 2 ) % 3;

  // Get pointers to the input and output
  typename Superclass::OutputImagePointer outputImage = this->GetOutput();
  typename Superclass::InputImagePointer inputImage = 
    const_cast< TInputImage * >( this->GetInput() );

  inputRegion = inputImage->GetLargestPossibleRegion ();

  inputIndex  = inputRegion.GetIndex();
  inputSize   = inputRegion.GetSize();

  inSpacing = inputImage->GetSpacing();

  outputIndex[0] = 0;
  outputIndex[1] = 0;


  bool canCrop = m_ROI.Crop( inputImage->GetLargestPossibleRegion() );
  if( !canCrop || !inputImage->GetLargestPossibleRegion().IsInside( m_ROI ) ) {
    std::cout << "Error: ROI outside the image." << std::endl;
    return;
  }


  // size of the window is reversed, enlarge the size
  typedef typename TOutputImage::SizeValueType  OutputSizeValueType;
  outputSize[0] = OutputSizeValueType( m_ROI.GetSize()[d2] * inSpacing[d2] * m_EnlargeBy );
  outputSize[1] = OutputSizeValueType( m_ROI.GetSize()[d1] * inSpacing[d1] * m_EnlargeBy );

  outSpacing[0] = inSpacing[d2];
  outSpacing[1] = inSpacing[d1];

  outOrigin[0] = inOrigin[d2];
  outOrigin[1] = inOrigin[d1];


  outputRegion.SetIndex( outputIndex );
  outputRegion.SetSize( outputSize );

  outputImage->SetLargestPossibleRegion( outputRegion );
  outputImage->SetSpacing( outSpacing );

  outputImage->SetOrigin( outOrigin);

  //outputImage->SetDirection( inputImage->GetDirection() );



  //itkDebugMacro("GenerateOutputInformation Start");

  //unsigned int projectionDimension = this->m_ImageViewer2->GetSliceOrientation();

  //typename TOutputImage::RegionType outputRegion;
  //typename TInputImage::IndexType inputIndex;
  //typename TInputImage::SizeType  inputSize;
  //typename TOutputImage::SizeType  outputSize;
  //typename TOutputImage::IndexType outputIndex;
  //typename TInputImage::SpacingType inSpacing;
  //typename TInputImage::PointType inOrigin;
  //typename TOutputImage::SpacingType outSpacing;
  //typename TOutputImage::PointType outOrigin;

  //// Get pointers to the input and output
  //typename Superclass::OutputImagePointer output = this->GetOutput();
  //typename Superclass::InputImagePointer input = 
  //  const_cast< TInputImage * >( this->GetInput() );

  //inputIndex = m_ROI.GetIndex();
  //inputSize = m_ROI.GetSize();
  //inSpacing = input->GetSpacing();
  //inOrigin = input->GetOrigin();

  //// Set the LargestPossibleRegion of the output.
  //// Reduce the size of the accumulated dimension.

  //if( static_cast< unsigned int >( InputImageDimension ) == 
  //    static_cast< unsigned int >( OutputImageDimension )    )
  //  {
  //  for(unsigned int i = 0; i<InputImageDimension; i++)
  //    {
  //    if (i != projectionDimension)
  //      {
  //      outputSize[i]  = inputSize[i];
  //      outputIndex[i] = inputIndex[i];
  //      outSpacing[i] = inSpacing[i];
  //      outOrigin[i]  = inOrigin[i];
  //      }
  //    else
  //      {
  //      outputSize[i]  = 1;
  //      outputIndex[i] = 0;
  //      outSpacing[i] = inSpacing[i]*inputSize[i];
  //      outOrigin[i]  = inOrigin[i] + (i-1)*inSpacing[i]/2;
  //      }
  //    }
  //  }
  //else
  //  {
  //  // Then OutputImageDimension = InputImageDimension - 1
  //  for(unsigned int i = 0; i<OutputImageDimension; i++)
  //    {
  //    if (i != projectionDimension)
  //      {
  //      outputSize[i]  = inputSize[i];
  //      outputIndex[i] = inputIndex[i];
  //      outSpacing[i] = inSpacing[i];
  //      outOrigin[i]  = inOrigin[i];
  //      }
  //    else
  //      {
  //      outputSize[i]  = inputSize[InputImageDimension - 1];
  //      outputIndex[i] = inputIndex[InputImageDimension - 1];
  //      outSpacing[i] = inSpacing[InputImageDimension - 1];
  //      outOrigin[i]  = inOrigin[InputImageDimension - 1];
  //      }
  //    }
  //  }

  //outputRegion.SetSize(outputSize);
  //outputRegion.SetIndex(outputIndex);
  //output->SetOrigin(outOrigin);
  //output->SetSpacing(outSpacing);
  //output->SetLargestPossibleRegion(outputRegion);

  //itkDebugMacro("GenerateOutputInformation End");
}


template <class TInputImage, class  TOutputImage>
void
ImageSlicesWithGeometryFilter<TInputImage,TOutputImage>
::GenerateInputRequestedRegion()
{
  itkDebugMacro("GenerateInputRequestedRegion Start");
  Superclass::GenerateInputRequestedRegion();


  typename TOutputImage::RegionType outputRegion;
  typename TInputImage::RegionType inputRegion;
  typename TInputImage::IndexType inputIndex;
  typename TInputImage::SizeType  inputSize;
  //typename TOutputImage::SizeType  outputSize;
  //typename TOutputImage::IndexType outputIndex;
  typename TInputImage::SpacingType inSpacing;
  typename TInputImage::PointType inOrigin;
  typename TOutputImage::SpacingType outSpacing;
  typename TOutputImage::PointType outOrigin;





  // Get pointers to the input and output
  typename Superclass::OutputImagePointer outputImage = this->GetOutput();
  typename Superclass::InputImagePointer inputImage = 
    const_cast< TInputImage * >( this->GetInput() );

  inputRegion = m_ROI;

  inputIndex  = inputRegion.GetIndex();
  inputSize   = inputRegion.GetSize();



  typename TInputImage::RegionType RequestedRegion;

  RequestedRegion.SetSize(inputSize);
  RequestedRegion.SetIndex(inputIndex);
  typename Superclass::InputImagePointer input = const_cast< TInputImage * > ( this->GetInput() );
  input->SetRequestedRegion (RequestedRegion);



  //unsigned int projectionDimension = this->m_ImageViewer2->GetSliceOrientation();

  //if ( this->GetInput() )
  //  {
  //  typename TInputImage::RegionType RequestedRegion;
  //  typename TInputImage::SizeType  inputSize;
  //  typename TInputImage::IndexType inputIndex;
  //  typename TInputImage::SizeType  inputLargSize;
  //  typename TInputImage::IndexType inputLargIndex;
  //  typename TOutputImage::SizeType  outputSize;
  //  typename TOutputImage::IndexType outputIndex;

  //  outputIndex = this->GetOutput()->GetRequestedRegion().GetIndex();
  //  outputSize = this->GetOutput()->GetRequestedRegion().GetSize();
  //  inputLargSize = this->GetInput()->GetLargestPossibleRegion().GetSize();
  //  inputLargIndex = this->GetInput()->GetLargestPossibleRegion().GetIndex();

  //  if( static_cast< unsigned int >( InputImageDimension ) == 
  //      static_cast< unsigned int >( OutputImageDimension )    )
  //    {
  //    for(unsigned int i=0; i<TInputImage::ImageDimension; i++)
  //      {
  //      if(i!=projectionDimension)
  //        {
  //        inputSize[i] = outputSize[i];
  //        inputIndex[i] = outputIndex[i];
  //        }
  //      else
  //        {
  //        inputSize[i]=inputLargSize[i];
  //        inputIndex[i]=inputLargIndex[i];
  //        }
  //      }
  //    }
  //  else
  //    {
  //    for(unsigned int i=0; i<OutputImageDimension; i++)
  //      {
  //      if(i!=projectionDimension)
  //        {
  //        inputSize[i] = outputSize[i];
  //        inputIndex[i] = outputIndex[i];
  //        }
  //      else
  //        {
  //        // the size of the output image on this dimension is the size
  //        // of the input image on the removed dimension
  //        inputSize[InputImageDimension - 1] = outputSize[i];
  //        inputIndex[InputImageDimension - 1] = outputIndex[i];
  //        }
  //      }
  //      inputSize[projectionDimension] = 
  //                   inputLargSize[projectionDimension];
  //      inputIndex[projectionDimension] = 
  //                   inputLargIndex[projectionDimension];
  //    }

  //  RequestedRegion.SetSize(inputSize);
  //  RequestedRegion.SetIndex(inputIndex);
  //  InputImagePointer input = const_cast< TInputImage * > ( this->GetInput() );
  //  input->SetRequestedRegion (RequestedRegion);
  //  }

  itkDebugMacro("GenerateInputRequestedRegion End");
}


template <class TInputImageType, class TOutputImageType>
void
ImageSlicesWithGeometryFilter<TInputImageType, TOutputImageType>::
GenerateData()
{
  // there must be a bug in this adaptor
  // unless it is created here again (creation in the constructor is not considered), the pipeline is not propagated properly
  // and the adaptor seems to have no output
  m_ToVTKAdaptorFilter = ToVTKAdaptorFilterType::New();
  m_ToVTKAdaptorFilter->SetInput( m_RegionOfInterestFilter->GetOutput() );


  // extract the image in the region of interest
  m_RegionOfInterestFilter->SetInput( this->GetInput() );
  m_RegionOfInterestFilter->SetRegionOfInterest( this->m_ROI );

  // convert the ITK image to VTK
  m_ToVTKAdaptorFilter->Update();
  m_ToVTKAdaptorFilter->GetOutput()->Update();

  
  // has to be before the extract geometry filter is initialized
  m_ImageViewer2->SetInput( m_ToVTKAdaptorFilter->GetOutput() );

  double *range = m_ToVTKAdaptorFilter->GetOutput()->GetScalarRange();
  int window = int( range[1] - range[0] );
  int level = int( 0.5 * ( range[1] + range[0] ) );
  m_ImageViewer2->SetColorLevel( level );
  m_ImageViewer2->SetColorWindow( window );

  //m_ImageViewer2->GetRenderer()->ResetCameraClippingRange();
  //m_ImageViewer2->GetRenderer()->ResetCamera( Box->GetBounds() );
  //viewer->GetRenderer()->ResetCameraClippingRange( bounds );
  //viewer->UpdateDisplayExtent();

  m_ImageViewer2->OffScreenRenderingOn();





  // prepare box for extract geometry filter
  //
  itk::ImageRegion<3>::IndexType  Start = m_ROI.GetIndex();
  itk::ImageRegion<3>::IndexType  End = m_ROI.GetIndex() + m_ROI.GetSize();

  // shrink size along which we are projecting to only a few slices, so that there is not too many points shown
  itk::ImageRegion<3>::SizeValueType  halfSize = itk::ImageRegion<3>::SizeValueType( m_ROI.GetSize()[m_SliceOrientation] / 2.0 );
  //itk::ImageRegion<3>::SizeValueType  shrinkFactor = 0;//m_ROI.GetSize()[m_SliceOrientation] / 2.0 - 3.0;
  Start[m_SliceOrientation] = m_ROI.GetIndex()[m_SliceOrientation] + halfSize - 1;
  End[m_SliceOrientation] = m_ROI.GetIndex()[m_SliceOrientation] + halfSize + 1;

  itk::Point< double, 3 >  StartPhysical;
  this->GetInput()->TransformIndexToPhysicalPoint( Start, StartPhysical );
  itk::Point< double, 3 >  EndPhysical;
  this->GetInput()->TransformIndexToPhysicalPoint( End, EndPhysical );

  vtkSmartPointer< vtkBox >  Box = vtkSmartPointer< vtkBox >::New();
  Box->SetBounds( StartPhysical[0], EndPhysical[0],
                  StartPhysical[1], EndPhysical[1],
                  StartPhysical[2], EndPhysical[2] );

  itkDebugMacro( "Box Bounds: " << StartPhysical[0] << ".." << EndPhysical[0] << ", "
                                << StartPhysical[1] << ".." << EndPhysical[1] << ", "
                                << StartPhysical[2] << ".." << EndPhysical[2] << std::endl );

  typename TInputImageType::SpacingType  spacing = this->GetInput()->GetSpacing();


  // if poly data supplied, add this actor and render them
  if( m_PolyData ) {
    // extract points inside a box
    m_ExtractGeometryFilter->SetInput( m_PolyData );
    m_ExtractGeometryFilter->SetImplicitFunction( Box );

    double scalarRange[2];
    m_PolyData->GetScalarRange( scalarRange );
    m_LookupTable->SetTableRange( scalarRange );
    m_PolyDataMapper->SetScalarRange( scalarRange );

    m_ImageViewer2->GetRenderer()->AddActor( m_Actor );
  }
  else {
    std::cout << "Warning: No polydata supplied." << std::endl;
  }

  // the other two dimensions
  int d1 = (m_SliceOrientation + 1) % 3;
  int d2 = (m_SliceOrientation + 2) % 3;


  // set the size -- this is in screen coordinates, but we use world coordinates (it does not matter for now)
  // remember, vertical direction is first
  //int winHeight = dimensions[d2];
  //int winWidth  = dimensions[d1];
  int winHeight = int( m_ROI.GetSize()[d2] * spacing[d2] );
  int winWidth  = int( m_ROI.GetSize()[d1] * spacing[d1] );
  m_ImageViewer2->SetSize( int( winWidth * m_EnlargeBy ), int( winHeight * m_EnlargeBy ) );  // size of the window is larger than the image
  //m_ImageViewer2->SetSize( winHeight * m_EnlargeBy, winWidth * m_EnlargeBy );  // size of the window is larger than the image

  // dimension of the orientation
  this->m_ImageViewer2->SetSliceOrientation( m_SliceOrientation );

  int sliceRange[2];
  m_ImageViewer2->GetSliceRange( sliceRange );

  int middleSlice = ( sliceRange[1] - sliceRange[0] + 1 ) / 2;
  m_ImageViewer2->SetSlice( middleSlice );

  // adjust the camera
  //
  vtkSmartPointer< vtkCamera >  camera = m_ImageViewer2->GetRenderer()->GetActiveCamera();
m_ImageViewer2->GetRenderer()->ResetCamera();
  //camera->ParallelProjectionOn();

  //int *dimensions = m_ToVTKAdaptorFilter->GetOutput()->GetDimensions();
  //double *origin = m_ToVTKAdaptorFilter->GetOutput()->GetOrigin();

  //double focalPoint[3];
  //double position[3];

  ////std::cout << "position: " << std::endl;
  //for ( unsigned int cc = 0; cc < 3; cc++)
  //  {
  //  focalPoint[cc] = origin[cc] + ( spacing[cc] * dimensions[cc] ) / 2.0;
  //  position[cc]   = focalPoint[cc];
  //  //std::cout << position[cc] << "  ";
  //  }
  ////std::cout << std::endl;

  //const double distanceToFocalPoint = -1000;
  //position[m_SliceOrientation] += distanceToFocalPoint;

  ////camera->SetPosition( position );
  //camera->SetFocalPoint( focalPoint );


 
  //double max = std::max( spacing[d1] * dimensions[d1], spacing[d2] * dimensions[d2] );

  //itkDebugMacro( "Max dim: " << max << " img size: " << m_ROI.GetSize()[d1] << " " << m_ROI.GetSize()[d2]
  //                           << " size*spacing: " << m_ROI.GetSize()[d1]*spacing[d1] << " " << m_ROI.GetSize()[d2]*spacing[d2] << std::endl );
  // parallel scale is w.r.t. height
  camera->SetParallelScale( winHeight / 2 );
  //camera->SetParallelScale( max / 2 );

  // FOR SOME REASON, THE WINDOW SIZE CHANGES AFTER RENDER EVEN THOUGH IMAGE EXTENT AND GEOMETRY BOUNDING BOX ARE SMALLER
  m_ImageViewer2->Render();
  //m_ImageViewer2->GetRenderWindow()->GetInteractor()->Start();

  m_WindowImageFilter->SetInput( m_ImageViewer2->GetRenderWindow() );
  vtkSmartPointer< vtkImageData >  imagePointsImage = m_WindowImageFilter->GetOutput();

  imagePointsImage->SetSpacing( spacing[0], spacing[1], spacing[2] );
	//imagePointsImage->SetOrigin( ROIOrigin[0], ROIOrigin[1], ROIOrigin[2] );
	imagePointsImage->SetScalarType( VTK_UNSIGNED_CHAR );
	imagePointsImage->SetNumberOfScalarComponents( 3 );

  // for debugging
  //vtkSmartPointer< vtkPNGWriter >  png = vtkSmartPointer< vtkPNGWriter >::New();
  //png->SetFileName( "frame.png" );
  //png->SetInput( imagePointsImage );
  //png->Write();

  m_ToITKAdaptorFilter->SetInput( imagePointsImage );
  m_ToITKAdaptorFilter->Update();

  typedef typename ToITKAdaptorFilterType::OutputImageType  ToITKAdaptorOutputType;
  ToITKAdaptorOutputType * output = const_cast<ToITKAdaptorOutputType*>( m_ToITKAdaptorFilter->GetOutput() );
  this->GraftOutput( output );
}


template <class TInputImageType, class TOutputImageType>
void
ImageSlicesWithGeometryFilter<TInputImageType, TOutputImageType>::
PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os,indent);

  os
    << indent << "ROI:" << this->m_ROI
    << std::endl;
}

} /* end namespace itk */


#endif
