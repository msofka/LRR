/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: ExtractFeatures.cxx,v $
  Language:  C++
  Date:      $Date: 2008/05/19 19:03:03 $
  Version:   $Revision: 1.27 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#pragma warning ( disable : 4996 )
#endif


//  Software Guide : BeginLatex
//
// \brief  Executable for extracting features from a Dicom image.
// \author Michal Sofka, sofka at cs dot rpi dot edu
// \date   Nov 2007
//
//  Software Guide : EndLatex 

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
//#include <direct.h>
#include <fstream>


// Software Guide : BeginCodeSnippet
#include "itkImageSeriesReader.h"
#include "itkImageSeriesWriter.h"
#include "itkNumericSeriesFileNames.h"

#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"

#include "itkImageFileReader.h"

#include "itkImage.h"

#include "itkTimeProbe.h"


//#include "itkVTKPolyDataWriter.h"  // there are some drawbacks in this class (under Review for ITK), e.g. input must be a triangle mesh
#include "itkPointSet.h"
#include "itkFeatureImageFilter.h"
#include "itkRegionOfInterestImageFilter.h"
#include "itkResampleImageFilter.h"

#include "vtkXMLPolyDataWriter.h"
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"

#include "vtkITKPointSetToPolyDataFilter.h"
//#include "vtkFeaturesToPolyDataFilter.h"

#include "itkMeshSpatialFilterClean.h"


int main( int argc, char * argv[] )
{

  if( argc < 3 )
    {
    std::cerr << "Missing Parameters " << std::endl;
    std::cerr << "Usage: " << argv[0];
    std::cerr << " ImageFile ";
    std::cerr << " FeaturesFile";
    std::cerr << " [SpacingForResampling (default is 1.0)]" << std::endl;
    return 1;
    }

  const     unsigned int   ImageDimension = 3;

  typedef   short  InputPixelType;
  typedef   itk::Image< InputPixelType, ImageDimension >  ImageType;


  // Try using image reader first (for mhd file type, for example).
  //
  typedef itk::ImageFileReader< ImageType > ImageReaderType;

  ImageReaderType::Pointer imageReader = ImageReaderType::New();
  imageReader->SetFileName( argv[1] );

  try
  {
    imageReader->Update();
  }
  catch( itk::ExceptionObject & excp )
  {
    std::cerr << "Exception thrown " << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
  }

  ImageType::Pointer inputImage = imageReader->GetOutput();



  itk::TimeProbe timer;
  timer.Start();

  if( inputImage.IsNull() ) {
    // Try using Dicom reader

    typedef   itk::ImageSeriesReader< ImageType >   ReaderType;

    // Read the input image.
    // USE GDCMImageIO, DICOMImageIO2 is OLD
    typedef itk::GDCMImageIO                        ImageIOType;
    typedef itk::GDCMSeriesFileNames                NamesGeneratorType;
    
    ImageIOType::Pointer io = ImageIOType::New();

    // Get the DICOM filenames from the directory
    NamesGeneratorType::Pointer  names = NamesGeneratorType::New();
    names->SetDirectory( argv[1] );
    const ReaderType::FileNamesContainer & fixedFileNames = 
                              names->GetInputFileNames();

    ReaderType::Pointer  reader = ReaderType::New();
    reader->SetFileNames( fixedFileNames );
    reader->SetImageIO( io );
    std::cout << "FIXED NAMES" << std::endl << names;

    try
      {
      reader->Update();
      }
    catch( itk::ExceptionObject & excp )
      {
      std::cerr << "Exception thrown " << std::endl;
      std::cerr << excp << std::endl;
      return EXIT_FAILURE;
      }


  // use region of interest
  #if 0
      // Setup region of interest filter.
      typedef  itk::RegionOfInterestImageFilter< ImageType, ImageType >  RegionInterestType;

      RegionInterestType::Pointer  regionInterest = RegionInterestType::New();

      //ImageType::SizeType  ROIsize = {{50, 50, 30}};
      //ImageType::IndexType  ROIindex = {{100, 100, 20}};
      ImageType::SizeType  ROIsize = {{100, 100, 40}};
      ImageType::IndexType  ROIindex = {{100, 100, 10}};
      ImageType::RegionType  region;
      region.SetSize( ROIsize );
      region.SetIndex( ROIindex );

      regionInterest->SetInput( reader->GetOutput() );
      regionInterest->SetRegionOfInterest( region );


      inputImage = regionInterest->GetOutput();
  #else
      inputImage = reader->GetOutput();
  #endif

  }


  // compute features
  //typedef itk::PointSet< InputPixelType,  ImageDimension>   PointsetType;
  typedef itk::PointAttribute  PointAttributeType;
  typedef itk::PointSet< PointAttributeType,  ImageDimension >   PointSetType;
  typedef itk::FeatureImageFilter< ImageType, PointSetType >  FeatureImageFilterType;
  FeatureImageFilterType::Pointer  featureFilter = FeatureImageFilterType::New();


#define RESAMPLE_TO_ISOTROPIC 1

#if RESAMPLE_TO_ISOTROPIC



// Software Guide : BeginCodeSnippet
  typedef itk::RecursiveGaussianImageFilter< 
                                ImageType,
                                ImageType > GaussianFilterType;
// Software Guide : EndCodeSnippet




// Software Guide : BeginLatex
//
// We create two instances of the smoothing filter, one will smooth along the
// $X$ direction while the other will smooth along the $Y$ direction. They are
// connected in a cascade in the pipeline, while taking their input from the
// intensity windowing filter. Note that you may want to skip the intensity
// windowing scale and simply take the input directly from the reader.
//
// Software Guide : EndLatex 

// Software Guide : BeginCodeSnippet
  GaussianFilterType::Pointer smootherX = GaussianFilterType::New();
  GaussianFilterType::Pointer smootherY = GaussianFilterType::New();

  smootherX->SetInput( inputImage );
  smootherY->SetInput( smootherX->GetOutput() );
// Software Guide : EndCodeSnippet


// Software Guide : BeginLatex
//
// We must now provide the settings for the resampling itself. This is done by
// searching for a value of isotropic resolution that will provide a trade-off
// between the evil of subsampling and the evil of supersampling. We advance
// here the conjecture that the geometrical mean between the in-plane and the
// inter-slice resolutions should be a convenient isotropic resolution to use.
// This conjecture is supported on nothing else than intuition and common
// sense. You can rightfully argue that this choice deserves a more technical
// consideration, but then, if you are so inclined to the technical correctness
// of the image sampling process, you should not be using this code, and should
// rather we talking about such technical correctness to the radiologist who
// acquired this ugly anisotropic dataset.
//
// We take the image from the input and then request its array of pixel spacing
// values.
//
// Software Guide : EndLatex 

// Software Guide : BeginCodeSnippet

  const ImageType::SpacingType& inputSpacing = inputImage->GetSpacing();
// Software Guide : EndCodeSnippet




// Software Guide : BeginLatex
//
// and apply our ad-hoc conjecture that the correct anisotropic resolution
// to use is the geometrical mean of the in-plane and inter-slice resolutions.
// Then set this spacing as the Sigma value to be used for the Gaussian
// smoothing at the preprocessing stage.
//
// Software Guide : EndLatex 

// Software Guide : BeginCodeSnippet

//  const double isoSpacing = sqrt( inputSpacing[2] * inputSpacing[0] );

  double isoSpacing = 1.0;
  if( argc > 3 ) {
    isoSpacing = atof( argv[3] );
  }

  
  smootherX->SetSigma( isoSpacing );
  smootherY->SetSigma( isoSpacing );
// Software Guide : EndCodeSnippet




// Software Guide : BeginLatex
//
// We instruct the smoothing filters to act along the $X$ and $Y$ direction
// respectively. And define the settings for avoiding the loss of intensity as
// a result of the diffusion process that is inherited from the use of a
// Gaussian filter.
//
// \index{RecursiveGaussianImageFilter!SetNormalizeAcrossScale}
// Software Guide : EndLatex 

 
// Software Guide : BeginCodeSnippet
  smootherX->SetDirection( 0 );
  smootherY->SetDirection( 1 );

  smootherX->SetNormalizeAcrossScale( true );
  smootherY->SetNormalizeAcrossScale( true );
// Software Guide : EndCodeSnippet


// Software Guide : BeginLatex
//
// Now that we have taken care of the smoothing in-plane, we proceed to
// instantiate the resampling filter that will reconstruct an isotropic image.
// We start by declaring the pixel type to be use at the output of such filter,
// then instantiate the image type and the type for the resampling filter.
// Finally we construct an instantiation of such a filter.
//
// Software Guide : EndLatex 


// Software Guide : BeginCodeSnippet
  typedef ImageType   OutputImageType;

  typedef itk::ResampleImageFilter<
                ImageType, OutputImageType >  ResampleFilterType;

  ResampleFilterType::Pointer resampler = ResampleFilterType::New();
// Software Guide : EndCodeSnippet




// Software Guide : BeginLatex
//
// The resampling filter requires that we provide a Transform, that in this
// particular case can simply be an identity transform.
//
// Software Guide : EndLatex 

// Software Guide : BeginCodeSnippet
  const     unsigned int    Dimension = 3;

  typedef itk::IdentityTransform< double, Dimension >  TransformType;

  TransformType::Pointer transform = TransformType::New();
  transform->SetIdentity();

  resampler->SetTransform( transform );
// Software Guide : EndCodeSnippet




// Software Guide : BeginLatex
//
// The filter also requires an interpolator to be passed to it. In this case we
// chose to use a linear interpolator.
//
// Software Guide : EndLatex 

// Software Guide : BeginCodeSnippet
  typedef itk::LinearInterpolateImageFunction< 
                          ImageType, double >  InterpolatorType;

  InterpolatorType::Pointer interpolator = InterpolatorType::New();

  resampler->SetInterpolator( interpolator );
// Software Guide : EndCodeSnippet




  resampler->SetDefaultPixelValue( 255 ); // highlight regions without source




// Software Guide : BeginLatex
//
// The pixel spacing of the resampled dataset is loaded in a \code{SpacingType}
// and passed to the resampling filter.
//
// Software Guide : EndLatex 

// Software Guide : BeginCodeSnippet
  OutputImageType::SpacingType spacing;

  spacing[0] = isoSpacing;
  spacing[1] = isoSpacing;
  spacing[2] = isoSpacing;

  resampler->SetOutputSpacing( spacing );
// Software Guide : EndCodeSnippet


// Software Guide : BeginLatex
//
// The origin of the output image is maintained, since we decided to resample
// the image in the same physical extent of the input anisotropic image.
//
// Software Guide : EndLatex 

// Software Guide : BeginCodeSnippet
  resampler->SetOutputOrigin( inputImage->GetOrigin() );
// Software Guide : EndCodeSnippet




// Software Guide : BeginLatex
//
// The number of pixels to use along each dimension in the grid of the
// resampled image is computed using the ratio between the pixel spacings of the
// input image and those of the output image. Note that the computation of the
// number of pixels along the $Z$ direction is slightly different with the
// purpose of making sure that we don't attempt to compute pixels that are
// outside of the original anisotropic dataset.
//
// Software Guide : EndLatex 

// Software Guide : BeginCodeSnippet
  ImageType::SizeType   inputSize = 
                    inputImage->GetLargestPossibleRegion().GetSize();
  
  typedef ImageType::SizeType::SizeValueType SizeValueType;

  const double dx = inputSize[0] * inputSpacing[0] / isoSpacing;
  const double dy = inputSize[1] * inputSpacing[1] / isoSpacing;

  const double dz = (inputSize[2] - 1 ) * inputSpacing[2] / isoSpacing;
// Software Guide : EndCodeSnippet



// Software Guide : BeginLatex
//
// Finally the values are stored in a \code{SizeType} and passed to the
// resampling filter. Note that this process requires a casting since the
// computation are performed in \code{double}, while the elements of the
// \code{SizeType} are integers.
//
// Software Guide : EndLatex 

// Software Guide : BeginCodeSnippet
  ImageType::SizeType   size;

  size[0] = static_cast<SizeValueType>( dx );
  size[1] = static_cast<SizeValueType>( dy );
  size[2] = static_cast<SizeValueType>( dz );

  resampler->SetSize( size );
// Software Guide : EndCodeSnippet


  std::cout << "Input spacing: " << inputSpacing[0] << ", " << inputSpacing[1] << ", " << inputSpacing[2] << std::endl;
  std::cout << "Resampling the image from: " << inputSize[0] << ", " << inputSize[1] << ", " << inputSize[2] << " to: "
                                             << size[0] << ", " << size[1] << ", " << size[2] << std::endl;



// Software Guide : BeginLatex
//
// Our last action is to take the input for the resampling image filter from
// the output of the cascade of smoothing filters, and then to trigger the
// execution of the pipeline by invoking the \code{Update()} method on the
// resampling filter.
//
// Software Guide : EndLatex 


// Software Guide : BeginCodeSnippet
  resampler->SetInput( smootherY->GetOutput() );




  ImageType::Pointer inputToFeatureDetector = resampler->GetOutput();

  // The resampler->SetOutputSpacing above does not seem to work -> set spacing on the image.
  inputToFeatureDetector->SetSpacing( spacing );

#else
  ImageType::Pointer inputToFeatureDetector = inputImage;
#endif

  featureFilter->SetInput( inputToFeatureDetector );

  ImageType::SpacingType  spacingForFeatureDet = inputToFeatureDetector->GetSpacing();
  float filteringRadiusMM[3] = {30.0, 30.0, 30.0}; // in milimeters
  FeatureImageFilterType::InputSizeType  filteringRadius;
  std::cout << "Spacing: " << spacingForFeatureDet << std::endl;
  filteringRadius[0] = filteringRadiusMM[0] / spacingForFeatureDet[0];
  filteringRadius[1] = filteringRadiusMM[1] / spacingForFeatureDet[1];
  filteringRadius[2] = filteringRadiusMM[2] / spacingForFeatureDet[2];

  featureFilter->SetFilteringRadius( filteringRadius ); // in pixels
  FeatureImageFilterType::OutputMeshPointer  features = featureFilter->GetOutput();
  featureFilter->Update();


  itk::TimeProbe timer1;
  timer1.Start();
  typedef itk::MeshSpatialFilterClean< PointSetType, PointSetType >  SpatialFilterType;
  SpatialFilterType::Pointer  spatialFilter = SpatialFilterType::New();
  spatialFilter->SetInput( featureFilter->GetOutput() );
  spatialFilter->SetMinDistanceBetweenPoints( 2.0f );
  spatialFilter->Update();
  features = spatialFilter->GetOutput();
//  features = featureFilter->GetOutput();
  timer1.Stop();
  std::cout << "Spatial filtering took " << timer1.GetMeanTime() << " seconds.\n";

  


  FeatureImageFilterType::PointDataContainerPointer  featurePointData = features->GetPointData();
  std::cout << "Extracted features: " << featurePointData->size() << std::endl;

  // count the number of each feature type
  unsigned int numCorners = 0;
  unsigned int numTubes = 0;
  unsigned int numSheets = 0;
  for( FeatureImageFilterType::PointDataContainer::VectorContainerSizeType  i = 0; i < featurePointData->size(); ++i ) {
    itk::PointAttribute::FeatureShape  featureShape = (*featurePointData)[i].m_Shape;
    if( featureShape == itk::PointAttribute::Corner ) ++numCorners;
    if( featureShape == itk::PointAttribute::Tube ) ++numTubes;
    if( featureShape == itk::PointAttribute::Sheet ) ++numSheets;
  }
  std::cout << "Extracted: " << numCorners << " corners, " << numTubes << " tubes, " << numSheets << " sheets." << std::endl;


//// THIS VTKFEATURE SET COULD BE ITK::POINTSET INSTEAD, IT CAN HAVE ATTRIBUTES, WHICH WOULD BE A CDCL FEATURE W/OUT LOCATION
//// THIS WOULD HELP TO DEFINE ITK::POINTSETREGISTRATION FILTER
//  vtkSmartPointer< vtkFeatureSet >  featureSet = vtkSmartPointer< vtkFeatureSet >::New();
//  //FeatureImageFilterType::PointDataContainerPointer  pointData = features->GetPointData();
//  FeatureImageFilterType::PointType  point;// = FeatureImageFilterType::PointType::New();
//  for( unsigned int i = 0; i < features->GetNumberOfPoints(); ++i ) {
//    //FeatureImageFilterType::PointType  point = pointData[i];
//    //itk::Point< InputPixelType, 3 >*  
//    features->GetPoint( i, &point );
//    featureSet->AddPoint( point );
//  }
//  


  timer.Stop();
  std::cout << "Feature extraction took " << timer.GetMeanTime() << " seconds.\n";


  // create ITK PointSet to VTK object connector
  vtkSmartPointer< vtkPointAttributeSet >  pointAttrubuteSetVTK = vtkSmartPointer< vtkPointAttributeSet >::New();
  pointAttrubuteSetVTK->SetInput( features );

  // create the converting filter
  vtkSmartPointer< vtkITKPointSetToPolyDataFilter >  pointSetToPolyDataFilter = vtkSmartPointer< vtkITKPointSetToPolyDataFilter >::New();
  pointSetToPolyDataFilter->SetInput( pointAttrubuteSetVTK->GetOutput() );
  pointSetToPolyDataFilter->Update();

  // write the features
  vtkSmartPointer< vtkXMLPolyDataWriter >  featureWriter = vtkSmartPointer< vtkXMLPolyDataWriter >::New();

  vtkSmartPointer< vtkPolyData >  poly_data = pointSetToPolyDataFilter->GetOutput();

  featureWriter->SetInput( poly_data );
  featureWriter->SetFileName( argv[2] );
  featureWriter->Write();




  try
    {
    featureWriter->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << "Exception thrown " << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
    }



  //// Write Debug image
  ////
  //typedef   short  OutputPixelType;
  //typedef   itk::Image< OutputPixelType, 2 >              SliceImageType;
  //
  //typedef   itk::ImageSeriesWriter< FeatureImageFilterType::InputImageType, SliceImageType >  WriterType;
  //WriterType::Pointer writer = WriterType::New();
  //std::string  debugFileName = "scores";
  ////_mkdir(argv[3]);
  //itksys::SystemTools::MakeDirectory( debugFileName.c_str() );
  //
  //names->SetOutputDirectory( debugFileName.c_str() );
  //writer->SetImageIO( io );
  // 
  ////const WriterType::FileNamesContainer &  movingMappedFileNames = 
  ////                          croppedNames->GetOutputFileNames();
  //WriterType::FileNamesContainer  writerFileNames = 
  //                          names->GetOutputFileNames();
  //// resize the names container (there can be more moving names than fixed names)
  //WriterType::FileNamesContainer::size_type  newSize = featureFilter->GetDebugImage()->GetLargestPossibleRegion().GetSize()[2];
  //writerFileNames.resize( newSize );
  //
  //writer->SetFileNames( writerFileNames );
  //
  //writer->SetMetaDataDictionaryArray( 
  //                      reader->GetMetaDataDictionaryArray() );
  //
  //typedef itk::RescaleIntensityImageFilter< FeatureImageFilterType::InternalImageType, FeatureImageFilterType::InputImageType >  RescalerType;
  ////  Rescaler to set the final image to 0..255
  //RescalerType::Pointer  rescaler = RescalerType::New();
  //rescaler->SetOutputMinimum(-10000);
  //rescaler->SetOutputMaximum(10000);
  //rescaler->SetInput( featureFilter->GetDebugImage() );
  //
  //
  //writer->SetInput( rescaler->GetOutput() );
  //
  //writer->Update();



  return EXIT_SUCCESS;
}

