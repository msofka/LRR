/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: ExtractKeypoints.cxx,v $
  Language:  C++
  Date:      $Date: 2008/07/22 22:26:45 $
  Version:   $Revision: 1.8 $

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
// \brief  Executable for extracting keypoints (corners) from a Dicom image.
//         This uses feature extraction filter, but only corners are extracted as keypoints.
//         Note the orientation of a keypoint is modulo \pi and unique orientation is assigned
//         during descriptor computation.
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

#include "itkImage.h"

//#include "itkVTKPolyDataWriter.h"  // there are some drawbacks in this class (under Review for ITK), e.g. input must be a triangle mesh
#include "itkPointSet.h"
#include "itkFeatureImageFilter.h"
#include "itkRegionOfInterestImageFilter.h"

#include "vtkPolyData.h"
#include "vtkPointData.h"
#include "vtkSmartPointer.h"
#include "vtkXMLPolyDataWriter.h"
#include "vtkXMLPolyDataReader.h"
//#include "vtkPolyDataToFeaturesWithShapeFilter.h"

#include "vtkITKPointSetToPolyDataFilter.h"
//#include "vtkFeaturesToPolyDataFilter.h"

#include "itkMeshSpatialFilterClean.h"
#include "itkVTKPolyDataToMeshFilter.h"

// Indicate whether to extract only corners during feature extraction or whether to extract all feature types
// and then only use corners as keypoints. Spatial filtering will be different in these cases.
#define EXTRACT_ONLY_CORNERS 0

int main( int argc, char * argv[] )
{

  if( argc < 3 )
    {
    std::cerr << "Missing Parameters " << std::endl;
    std::cerr << "Usage: " << argv[0];
    std::cerr << " FeatureFile ";
    std::cerr << " KeypointFile" << std::endl;
    return 1;
    }

  const     unsigned int   ImageDimension = 3;

#if EXTRACT_ONLY_CORNERS


  typedef   short  InputPixelType;
  typedef   itk::Image< InputPixelType, ImageDimension >  ImageType;

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

  itk::TimeProbe timer;
  timer.Start();

  // compute features
  //typedef itk::PointSet< InputPixelType,  ImageDimension>   PointsetType;
  typedef itk::PointAttribute  PointAttributeType;
  typedef itk::PointSet< PointAttributeType,  ImageDimension >   PointSetType;
  typedef itk::FeatureImageFilter< ImageType, PointSetType >  FeatureImageFilterType;
  FeatureImageFilterType::Pointer  featureFilter = FeatureImageFilterType::New();

// use region of interest
#if 0
    // Setup region of interest filter.
    typedef  itk::RegionOfInterestImageFilter< ImageType, ImageType >  RegionInterestType;

    RegionInterestType::Pointer  regionInterest = RegionInterestType::New();

    //ImageType::SizeType  size = {{50, 50, 30}};
    //ImageType::IndexType  index = {{100, 100, 20}};
    ImageType::SizeType  size = {{100, 100, 40}};
    ImageType::IndexType  index = {{100, 100, 10}};
    ImageType::RegionType  region;
    region.SetSize( size );
    region.SetIndex( index );

    regionInterest->SetInput( reader->GetOutput() );
    regionInterest->SetRegionOfInterest( region );


    featureFilter->SetInput( regionInterest->GetOutput() );
#else
    featureFilter->SetInput( reader->GetOutput() );
#endif

  ImageType::SpacingType  spacing = reader->GetOutput()->GetSpacing();
  float filteringRadiusMM[3] = {30.0, 30.0, 30.0}; // in milimeters
  FeatureImageFilterType::InputSizeType  filteringRadius;
  filteringRadius[0] = filteringRadiusMM[0] / spacing[0];
  filteringRadius[1] = filteringRadiusMM[1] / spacing[1];
  filteringRadius[2] = filteringRadiusMM[2] / spacing[2];

  featureFilter->SetOnlyCorners( true );
  featureFilter->SetFilteringRadius( filteringRadius ); // in pixels
  FeatureImageFilterType::OutputMeshPointer  features = featureFilter->GetOutput();
  featureFilter->Update();


  typedef itk::MeshSpatialFilterClean< PointSetType, PointSetType >  SpatialFilterType;
  SpatialFilterType::Pointer  spatialFilter = SpatialFilterType::New();
  spatialFilter->SetInput( featureFilter->GetOutput() );
  spatialFilter->SetMinDistanceBetweenPoints( 10.0f );

  spatialFilter->Update();
  features = spatialFilter->GetOutput();
//  features = featureFilter->GetOutput();


  timer.Stop();
  FeatureImageFilterType::PointDataContainerPointer  featurePointData = features->GetPointData();
  std::cout << "Extracted keypoints: " << featurePointData->size() << " took: " << timer.GetMeanTime() << std::endl;

#else
  vtkSmartPointer< vtkXMLPolyDataReader >  featureReader = vtkSmartPointer< vtkXMLPolyDataReader >::New();
  std::string  filename = std::string( argv[1] );
  featureReader->SetFileName( filename.c_str() );
  featureReader->Update();  // need this update (bug in the converting filter)


  typedef itk::PointAttribute  PointAttributeType;
  typedef itk::PointSet< PointAttributeType,  ImageDimension >   PointSetType;

  // create the feature converting filter
  typedef itk::VTKPolyDataToMeshFilter< PointSetType >  FeatureConvertorType;
  FeatureConvertorType::Pointer  featureConverter = FeatureConvertorType::New();
  featureConverter->SetInput( featureReader->GetOutput() );
  featureConverter->Update();

  FeatureConvertorType::OutputMeshPointer  featuresAll = featureConverter->GetOutput();

  // Go through all features and get only corners.
  // First create the new mesh.
  //
  typedef FeatureConvertorType::PointsContainer  PointsContainer;
  typedef FeatureConvertorType::PointDataContainer  PointDataContainer;
  typedef FeatureConvertorType::PointsContainerPointer  PointsContainerPointer;
  typedef FeatureConvertorType::PointDataContainerPointer  PointDataContainerPointer;

  PointsContainerPointer      points    = featuresAll->GetPoints();
  PointDataContainerPointer   pointData = featuresAll->GetPointData();

  PointsContainerPointer     corners     = PointsContainer::New();
  PointDataContainerPointer  cornerData  = PointDataContainer::New();

  typedef itk::PointSet< PointAttributeType,  ImageDimension >   PointSetType;
  PointSetType::Pointer  cornerFeatures = PointSetType::New();
  cornerFeatures->SetPoints( corners );
  cornerFeatures->SetPointData( cornerData );

  assert( points->size() == pointData->size() );


  PointDataContainer::ConstIterator  itData = pointData->Begin();
  for( PointsContainer::ConstIterator  it = points->Begin(); it != points->End(); ++it, ++itData ) {
    if( itData->Value().m_Shape == PointAttributeType::Corner ) {
      corners->push_back( it->Value() );
      cornerData->push_back( itData->Value() );
    }
  }


  // Spatial filtering with radius larger than for feature extraction.
  typedef itk::MeshSpatialFilterClean< PointSetType, PointSetType >  SpatialFilterType;
  SpatialFilterType::Pointer  spatialFilter = SpatialFilterType::New();
  spatialFilter->SetInput( cornerFeatures );
  spatialFilter->SetMinDistanceBetweenPoints( 6.0f );
  spatialFilter->Update();
  PointSetType::Pointer  features = spatialFilter->GetOutput();



#endif

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

  //typedef   itk::ImageSeriesWriter< FeatureImageFilterType::InternalImageType, SliceImageType >  WriterType;
  //WriterType::Pointer writer = WriterType::New();
  //std::string  debugFileName = "scores";
  ////_mkdir(argv[3]);
  //itksys::SystemTools::MakeDirectory( debugFileName.c_str() );

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

  //writer->SetFileNames( writerFileNames );

  //writer->SetMetaDataDictionaryArray( 
  //                      reader->GetMetaDataDictionaryArray() );

  //typedef itk::RescaleIntensityImageFilter< FeatureImageFilterType::InternalImageType, FeatureImageFilterType::InternalImageType >  RescalerType;
  ////  Rescaler to set the final image to 0..255
  //RescalerType::Pointer  rescaler = RescalerType::New();
  //rescaler->SetOutputMinimum(-10000);
  //rescaler->SetOutputMaximum(10000);
  //rescaler->SetInput( featureFilter->GetDebugImage() );


  //writer->SetInput( rescaler->GetOutput() );

  //writer->Update();



  return EXIT_SUCCESS;
}

