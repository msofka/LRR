/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: CropDicom.cxx,v $
  Language:  C++
  Date:      $Date: 2008/08/22 22:58:42 $
  Version:   $Revision: 1.21 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#if defined(_MSC_VER)
#pragma warning ( disable : 4996 )
#pragma warning ( disable : 4786 )
#endif


//  Software Guide : BeginLatex
//
//  \brief
//  This example illustrates how to crop Dicom lung CT volume to remove air around the body.
//  Ideas from: "Lung CT segmentation for image retrieval using the Insight Toolkit (ITK)",
//  Joris Heuberger, Antoine Geissbuhler, Henning Muller
//  Medical Imaging and Telemedicine (MIT 2005), pages 57-62, WuyiShan, China, August 2005.
//
// \author Michal Sofka, sofka at cs dot rpi dot edu
// \date   Nov 2007
//
//  Software Guide : EndLatex 

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
//#include <direct.h>


// Software Guide : BeginCodeSnippet
#include "itkImageSeriesReader.h"
#include "itkImageSeriesWriter.h"
#include "itkDICOMImageIO2.h"
#include "itkDICOMSeriesFileNames.h"
#include "itkNumericSeriesFileNames.h"

#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkMetaDataObject.h"
#include "itkMetaDataObjectBase.h"

#include "itkImage.h"
#include "itkResampleImageFilter.h"
#include "itkLinearInterpolateImageFunction.h"

#include "itkBSplineDeformableTransform.h"
#include "itkAffineTransform.h"
#include "itkTransformFileWriter.h"
#include "itkTransformFileReader.h"
#include "itkTransformFactory.h"

#include "itkBinaryThresholdImageFilter.h"
#include "itkNeighborhoodConnectedImageFilter.h"
#include "itkStatisticsImageFilter.h"
#include "itkImageRegionConstIterator.h"

#include "itkExtractImageFilter.h"
#include "itkRegionOfInterestImageFilter.h"
#include "itkBoundingBox.h"

#include <fstream>


int main( int argc, char * argv[] )
{

  if( argc < 3 )
    {
    std::cerr << "Missing Parameters " << std::endl;
    std::cerr << "Usage: " << argv[0];
    std::cerr << " Image ";
    std::cerr << " croppedImage";
    std::cerr << " [noDecrease]" << std::endl;
    return 1;
    }

  const     unsigned int   ImageDimension = 3;

  typedef   short  InputPixelType;
  typedef   short  OutputPixelType;

  typedef   itk::Image< InputPixelType, ImageDimension >  ImageType;
  typedef   itk::Image< OutputPixelType, 2 >              SliceImageType;

  typedef   itk::ImageSeriesReader< ImageType >   ReaderType;

  typedef   itk::ImageSeriesWriter< ImageType, SliceImageType >  WriterType;


  // Read the input image.
  // USE GDCMImageIO, DICOMImageIO2 is OLD
  typedef itk::GDCMImageIO                        ImageIOType;
  typedef itk::GDCMSeriesFileNames                NamesGeneratorType;
  //typedef itk::DICOMImageIO2                      ImageIOType;
  //typedef itk::DICOMSeriesFileNames               NamesGeneratorType;
  
  ImageIOType::Pointer io = ImageIOType::New();

  // Get the DICOM filenames from the directory
  NamesGeneratorType::Pointer  names = NamesGeneratorType::New();
  names->SetDirectory( argv[1] );
  const ReaderType::FileNamesContainer & fixedFileNames = 
                            names->GetInputFileNames();

  ReaderType::Pointer  reader = ReaderType::New();
  reader->SetFileNames( fixedFileNames );
  reader->SetImageIO( io );
  // WRITER DOES NOT HAVE ReverseOrderOn FUNCTION! -> It would get flipped
  //reader->ReverseOrderOn();
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

  ImageType::RegionType fixedImageRegion = reader->GetOutput()->GetBufferedRegion();
  ImageType::SizeType   fixedImageSize = fixedImageRegion.GetSize();

  WriterType::Pointer writer = WriterType::New();

  //_mkdir(argv[3]);
  itksys::SystemTools::MakeDirectory( argv[2] );

  names->SetOutputDirectory( argv[2] );
  writer->SetImageIO( io );
  
  
//  Software Guide : BeginLatex
//
//  The following line of code is extremely important for this process to work
//  correctly.  The line is taking the MetaDataDictionary from the input reader
//  and passing it to the output writer. The reason why this step is so
//  important is that the MetaDataDictionary contains all the entries of the
//  input DICOM header.
//
//  \index{itk::ImageSeriesReader!GetMetaDataDictionaryArray()}
//  \index{itk::ImageSeriesWriter!SetMetaDataDictionaryArray()}
//
//  Software Guide : EndLatex 

// Software Guide : BeginCodeSnippet
  //writer->SetMetaDataDictionaryArray( 
  //                      reader->GetMetaDataDictionaryArray() );
// Software Guide : EndCodeSnippet


  ImageType::ConstPointer image = reader->GetOutput();

  std::cout << "Read image origin: " << image->GetOrigin() << std::endl;

  // find optimal threshold for the background / tissue separation
  //
  typedef itk::BinaryThresholdImageFilter< ImageType, ImageType >  ThresholderType;
  ThresholderType::Pointer  thresholder = ThresholderType::New();
  thresholder->SetInput( image );

  //typedef itk::StatisticsImageFilter< ImageType >  StatisticsFilterType;
  //StatisticsFilterType::Pointer  statisticsFilter = StatisticsFilterType::New();
  //statisticsFilter->SetInput( thresholder->GetOutput() );

  typedef itk::ImageRegionConstIterator< ImageType >  ConstIteratorType; 

  // find a threshold which separates background and forground means the most
  ThresholderType::InputPixelType  threshold = 128;
  ThresholderType::InputPixelType  thresholdPrevious = 128;
  bool converged = false;
  while( !converged ) {

    ConstIteratorType  it( image, image->GetLargestPossibleRegion() );

    double backgroundMean = 0.0;
    double tissueMean = 0.0;
    unsigned int backgroundCount = 0;
    unsigned int tissueCount = 0;
    for( it.GoToBegin(); ! it.IsAtEnd(); ++it ) {

      ImageType::PixelType  pixel = it.Value();
// common background value in the scans appearing as a tube
if( pixel == -3024 ) continue;
      if( pixel < threshold ) {
        backgroundMean += pixel;
        ++backgroundCount;
      }
      else {
        tissueMean += pixel;
        ++tissueCount;
      }
    }
    backgroundMean /= double( backgroundCount );
    tissueMean /= double( tissueCount );

    thresholdPrevious = threshold;
    threshold = ThresholderType::InputPixelType( ( tissueMean + backgroundMean ) / 2 );

    if( vcl_abs( thresholdPrevious - threshold ) < 1 ) converged = true;

    vcl_cout << "tissue mean: " << tissueMean << "  back mean: " << backgroundMean << "  thresh: " << threshold << vcl_endl;
  }

  thresholder->SetLowerThreshold( threshold );
  thresholder->SetUpperThreshold( itk::NumericTraits<ThresholderType::InputPixelType>::max() );
  //thresholder->SetLowerThreshold( itk::NumericTraits<ThresholderType::InputPixelType>::NonpositiveMin() );
  //thresholder->SetUpperThreshold( threshold );
  thresholder->SetOutsideValue( itk::NumericTraits<ThresholderType::InputPixelType>::NonpositiveMin() );


  // find background connected component
  //
  typedef itk::NeighborhoodConnectedImageFilter< ImageType, ImageType >  ConnectedFilterType;
  ConnectedFilterType::Pointer  connectedFilter = ConnectedFilterType::New();

  ThresholderType::OutputImageType::Pointer  thresholdedImage = thresholder->GetOutput();

  ImageType::RegionType  imageRegion = image->GetLargestPossibleRegion();
  ImageType::SizeType  imageSize = imageRegion.GetSize();
  ConnectedFilterType::IndexType  startIndex = imageRegion.GetIndex();
  ConnectedFilterType::IndexType  endIndex_plusOne = startIndex + imageSize;
  for( ConnectedFilterType::IndexType::IndexValueType  x = startIndex[0]; x < endIndex_plusOne[0]; ++x ) {
    for( ConnectedFilterType::IndexType::IndexValueType  z = startIndex[2]; z < endIndex_plusOne[2]; ++z ) {
  // add seeds as the bounding box of the volume (-> seeds will start on the outside)
  //for( ConnectedFilterType::IndexType::IndexValueType  x = startIndex[0]+0; x < startIndex[0]+11; ++x ) {
  //  for( ConnectedFilterType::IndexType::IndexValueType  z = startIndex[2]+0; z < startIndex[2]+11; ++z ) {
      ConnectedFilterType::IndexType  seed;
      seed[0] = x;
      seed[1] = startIndex[1];
      seed[2] = z;
      connectedFilter->AddSeed( seed );
      //vcl_cout << "Seed: " << seed << vcl_endl;
    }
  }


  ConnectedFilterType::InputImageSizeType  radius;
  for( unsigned int d = 0; d < ImageDimension; ++d ) {
    radius[d] = 1;
  }
  connectedFilter->SetRadius( radius );
  connectedFilter->SetLower( itk::NumericTraits<ThresholderType::InputPixelType>::NonpositiveMin() );
  connectedFilter->SetUpper( threshold );
  ImageType::PixelType  replaceValue = 100;
  connectedFilter->SetReplaceValue( replaceValue );
  connectedFilter->SetInput( image );
  
  connectedFilter->Update();
  ImageType::Pointer  connectedComponent = connectedFilter->GetOutput();

  // Find bounding box of the foreground region (not background component found above).
  // Do this by: for each slice step from the center towards the boundaries (on a cross), collect points
  // until reaching background and compute bounding box using those points.
  //
  typedef  itk::BoundingBox< >  BoundingBoxType;
  BoundingBoxType::Pointer  boundingBox = BoundingBoxType::New();
  typedef BoundingBoxType::PointType  PointType;
  BoundingBoxType::PointsContainerPointer  points = BoundingBoxType::PointsContainer::New();
  // for each slice
  for( ImageType::IndexType::IndexValueType  z = startIndex[2]; z < endIndex_plusOne[2]; ++z ) {
    // x center to max
    for( ImageType::IndexType::IndexValueType  x = startIndex[0] + imageSize[0] / 2; x < endIndex_plusOne[0]; ++x ) {
      ImageType::IndexType  index;
      index[0] = x;
      index[1] = startIndex[1] + imageSize[1] / 2;
      index[2] = z;
      PointType  point;
      //std::cout << "x: " << x << "  " << connectedComponent->GetPixel( index ) << std::endl;
      if( connectedComponent->GetPixel( index ) != replaceValue ) {
        connectedComponent->TransformIndexToPhysicalPoint( index, point );
        points->push_back( point );
      }
      else break;
    }

    // x center to min
    for( ImageType::IndexType::IndexValueType  x = startIndex[0] + imageSize[0] / 2; x >= startIndex[0]; --x ) {
      ImageType::IndexType  index;
      index[0] = x;
      index[1] = startIndex[1] + imageSize[1] / 2;
      index[2] = z;
      PointType  point;
      if( connectedComponent->GetPixel( index ) != replaceValue ) {
        connectedComponent->TransformIndexToPhysicalPoint( index, point );
        points->push_back( point );
      }
      else break;
    }

    // y center to max
    for( ImageType::IndexType::IndexValueType  y = startIndex[1] + imageSize[1] / 2; y < endIndex_plusOne[1]; ++y ) {
      ImageType::IndexType  index;
      index[0] = startIndex[0] + imageSize[0] / 2;
      index[1] = y;
      index[2] = z;
      PointType  point;
      if( connectedComponent->GetPixel( index ) != replaceValue ) {
        connectedComponent->TransformIndexToPhysicalPoint( index, point );
        points->push_back( point );
      }
      else break;
    }

    // y center to min
    for( ImageType::IndexType::IndexValueType  y = startIndex[1] + imageSize[1] / 2; y >= startIndex[1]; --y ) {
      ImageType::IndexType  index;
      index[0] = startIndex[0] + imageSize[0] / 2;
      index[1] = y;
      index[2] = z;
      PointType  point;
      if( connectedComponent->GetPixel( index ) != replaceValue ) {
        connectedComponent->TransformIndexToPhysicalPoint( index, point );
        points->push_back( point );
      }
      else break;
    }

  }
  boundingBox->SetPoints( points );

  PointType  minPoint = boundingBox->GetMinimum();
  PointType  maxPoint = boundingBox->GetMaximum();
  ImageType::IndexType  minIndex;
  image->TransformPhysicalPointToIndex( minPoint, minIndex );
  ImageType::IndexType  maxIndex;
  image->TransformPhysicalPointToIndex( maxPoint, maxIndex );

  if( !(argc == 4 && std::string( argv[3] ) == "noDecrease") ) {
    std::cout << "Decreasing the size of the bounding box .... " << std::endl;
    // decrease the size of the bounding box
    minIndex[0] += 100;
    minIndex[1] += 100;
    minIndex[2] += 20;

    maxIndex[0] -= 100;
    maxIndex[1] -= 100;
    maxIndex[2] -= 20;
  }

  ImageType::SizeType  boundingBoxSize;
  boundingBoxSize[0] = maxIndex[0] - minIndex[0];
  boundingBoxSize[1] = maxIndex[1] - minIndex[1];
  boundingBoxSize[2] = maxIndex[2] - minIndex[2];

#if 0
  typedef itk::RegionOfInterestImageFilter< ImageType, ImageType >  ROIFilterType;
  ROIFilterType::Pointer  extractFilter = ROIFilterType::New();

  ROIFilterType::InputImageRegionType  extractionRegion;

	extractFilter->SetInput( image );

  extractionRegion.SetSize( boundingBoxSize );
  extractionRegion.SetIndex( minIndex );

	extractFilter->SetRegionOfInterest( extractionRegion );
	extractFilter->Update();
#else

  typedef itk::ExtractImageFilter< ImageType, ImageType >  ExtractFilterType;

  ExtractFilterType::Pointer  extractFilter = ExtractFilterType::New();
  ExtractFilterType::InputImageRegionType  extractionRegion;

  extractFilter->SetInput( image );

  extractionRegion.SetSize( boundingBoxSize );
  extractionRegion.SetIndex( minIndex );

  extractFilter->SetExtractionRegion( extractionRegion );
  extractFilter->Update();
#endif


  WriterType::FileNamesContainer  croppedFileNames = 
                            names->GetOutputFileNames();
  croppedFileNames.resize( boundingBoxSize[2] );
  writer->SetFileNames( croppedFileNames );

  std::cout << "Cropped names, size: " << croppedFileNames.size() << std::endl;
  for( WriterType::FileNamesContainer::size_type  i = 0; i < croppedFileNames.size(); ++i ) {
    std::cout << croppedFileNames[i] << std::endl;
  }

  ImageType::Pointer  extractedImage = extractFilter->GetOutput();
  extractedImage->DisconnectPipeline();


  std::cout << "Extracted image origin: " << extractedImage->GetOrigin() << std::endl;
  std::cout << "Extracted image size: " << extractedImage->GetLargestPossibleRegion().GetSize() << std::endl;

	std::cout << "WARNING!!! Make sure to build this executable with itkGDCMImageIO.cxx revision 1.141 (otherwise Image Origin gets corrupted). " << std::endl;
	std::cout << "  There was a bug introduced in a later version. See: http://public.kitware.com/Bug/view.php?id=0006637" << std::endl;


  PointType  newOriginPoint;
  image->TransformIndexToPhysicalPoint( minIndex, newOriginPoint );

  //itk::Point< double, 3 >  newOriginPoint = extractedImage->GetOrigin();


  double newOrigin[3] = { newOriginPoint[0], newOriginPoint[1], newOriginPoint[2] };
  //extractedImage->SetOrigin( newOrigin );

  writer->SetInput( extractedImage );


//  Software Guide : EndCodeSnippet
//
//  std::cout << "Origin should be (convert min index from extracted region to a point): " << newOriginPoint << std::endl;
//  std::cout << "(If it isn't -> ITK bug -- origin should be modified in ExtractImageFilter) -- however, largest possible region is properly modified and origin is properly saved." << std::endl;
//
//
//
//  //std::cout << "Warning, this probably does not modify DICOM -- need to change the tags, doing that now... (but maybe this bug was fixed)" << std::endl;
//  //Update: the problem was in writing series (reported as a bug). Keeping this code only as a reference.
//
//
//  //typedef itk::MetaDataDictionary   DictionaryType;
//  //DictionaryType & dictionary = io->GetMetaDataDictionary();
//
//  //std::string  tagkey = "0020|0032";
//  //
//  //std::ostringstream  valueStream;
//  //valueStream << newOrigin[0] << "\\" << newOrigin[1] << "\\" << newOrigin[2];
//  //std::string  value = valueStream.str();
//
//
//  ////MetaDataObjectBase::Pointer  metaObject = dictionary[tagkey];
//  ////metaObject->SetMetaDataObjectValue( value );
//
//  //DictionaryType::Iterator  keyToChange = dictionary.Find( tagkey );
//  //*keyToChange = value;
//
//
//  //itk::EncapsulateMetaData<std::string>( dictionary, tagkey, value );
//  //std::cout << "FIXMEFIXME here this thing is not changed. I think that we need tagkey w/out | and vaule separated by backslashes (use ExtractSlice to see the tags, see DicomImageReadChangeHeaderWrite)." << std::endl;
//
//
//  //io->SetMetaDataDictionary( dictionary );
//  //writer->SetImageIO( io );
//
//
//
//
//  std::string  patientDir = itksys::SystemTools::GetParentDirectory( argv[1] );
//  std::string  nodulePath = patientDir + "/nodules.txt";
//  std::string  newNodulePath = patientDir + "/nodules_whole.txt";
//
//  std::cout << "Reading nodule file: " << nodulePath << std::endl;
//  std::ifstream  noduleFile( nodulePath.c_str() );
//  std::ofstream  newNoduleFile( newNodulePath.c_str() );
//  if( !noduleFile ) {
//    std::cout << "Specified nodule file cannot be opened: " << noduleFile << std::endl;
//    return 1;
//  }
//  else {
//
//    itk::Point< double, 3 >  pointRAS;
//    while( noduleFile >> pointRAS[0] >> pointRAS[1] >> pointRAS[2] ) {
//      std::cout << "Read nodule location: " << pointRAS << std::endl;
//
//      pointRAS[0] += ( image->GetOrigin()[0] - newOrigin[0] );
//      pointRAS[1] += ( image->GetOrigin()[1] - newOrigin[1]);
//      pointRAS[2] += ( image->GetOrigin()[2] - newOrigin[2] );
//
//      newNoduleFile << pointRAS[0] << " " << pointRAS[1] << " " << pointRAS[2] << std::endl;
//    }
//
//    newNoduleFile.close();
//  }
//  noduleFile.close();

  
  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << "Exception thrown " << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
    }
//  Software Guide : EndCodeSnippet




  return EXIT_SUCCESS;
}

