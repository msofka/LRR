/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: ExtractSlice.cxx,v $
  Language:  C++
  Date:      $Date: 2008/03/13 23:32:41 $
  Version:   $Revision: 1.3 $

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
// \brief  Executable for extracting a slice (ROI) around a given volume.
// \author Michal Sofka, sofka at cs dot rpi dot edu
// \date   Mar 2008
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

#include "itkRegionOfInterestImageFilter.h"

#include "itkRescaleIntensityImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkMaximumProjectionImageFilter.h"
#include "itkResampleImageFilter.h"
#include "itkTileImageFilter.h"
#include "itkPNGImageIO.h"

#include "itkOrientedImage.h"



namespace itk {

template <class ImageType>
void WriteROIImage( typename ImageType::Pointer const &  Image,
                    itk::Point< double, 3 > const &  Point,
                    std::string  FileName )
{
  typedef typename itk::ImageRegion<3>  RegionOfInterestType;
  typedef itk::Point< double, 3 >  PointType;
  
  RegionOfInterestType  m_ROI;

  RegionOfInterestType::SizeType  size;
  RegionOfInterestType::IndexType  index;
  Image->TransformPhysicalPointToIndex( Point, index );

  std::cout << "Index of Point: " << index << std::endl;

  typename ImageType::SpacingType  spacing = Image->GetSpacing();

  size[0] = 100 / spacing[0];
  size[1] = 100 / spacing[1];
  size[2] = 100 / spacing[2];

  index[0] = index[0] - size[0] / 2;
  index[1] = index[1] - size[1] / 2;
  index[2] = index[2] - size[2] / 2;

  m_ROI.SetSize( size );
  m_ROI.SetIndex( index );

  bool canCrop = m_ROI.Crop( Image->GetLargestPossibleRegion() );
  if( !canCrop || !Image->GetLargestPossibleRegion().IsInside( index ) ) {
    std::cout << "Error: point outside the image. Or the ROI outside of image." << std::endl;
    return;
  }


  // GenerateSlices
  //
  typedef itk::Image<signed short, 2>                           SliceImageType;
  typedef itk::ExtractImageFilter< ImageType, SliceImageType >  ExtractFilterType;
  typedef itk::ExtractImageFilter< ImageType, ImageType >       ExtractFilter3DType;
  typedef itk::MaximumProjectionImageFilter< ImageType, SliceImageType >  MaximumProjectionFilterType;

  // for some reason, having one filter and just changing the region didn't work here
  //
  // XY plane
  //
  SliceImageType::Pointer  movingSliceXY;
  {
  //ExtractFilter3DAdaptorType::Pointer  extractFilterXY = ExtractFilter3DAdaptorType::New();
  //ExtractFilter3DAdaptorType::InputImageRegionType  extractionRegion;

  //extractFilterXY->SetInput( m_MovingMappedImage );
  //extractFilterXY->SetInput( m_PointsImage );

  //// using m_ROI -- since m_MovingMappedImage is in the fixed coordinate system
  //ExtractFilter3DAdaptorType::InputImageRegionType::SizeType  size = m_ROI.GetSize();
  //ExtractFilter3DAdaptorType::InputImageRegionType::IndexType  index = m_ROI.GetIndex();
  //// add half size to get the slice in the middle of the ROI
  //index[2] += size[2] / 2 - 3;
  //// extracting Z slice
  //size[2] = 6;
  typename ExtractFilter3DType::Pointer  extractFilterXY = ExtractFilter3DType::New();
  typename ExtractFilter3DType::InputImageRegionType  extractionRegion;

  extractFilterXY->SetInput( Image );

  // using m_ROI -- since m_MovingMappedImage is in the fixed coordinate system
  typename ExtractFilter3DType::InputImageRegionType::SizeType  size = m_ROI.GetSize();
  typename ExtractFilter3DType::InputImageRegionType::IndexType  index = m_ROI.GetIndex();
  // add half size to get the slice in the middle of the ROI
  index[2] += size[2] / 2;
  // extracting Z slice
  size[2] = 1;

  extractionRegion.SetSize( size );
  extractionRegion.SetIndex( index );
  bool canCrop = extractionRegion.Crop( Image->GetLargestPossibleRegion() );
  if( !canCrop || !Image->GetLargestPossibleRegion().IsInside( index ) ) {
    std::cout << "Error: point outside the image. Or the ROI outside of image." << std::endl;
    return;
  }


  extractFilterXY->SetExtractionRegion( extractionRegion );

  typename MaximumProjectionFilterType::Pointer  maximumProjectionFilter = MaximumProjectionFilterType::New();
  maximumProjectionFilter->SetInput( extractFilterXY->GetOutput() );
  maximumProjectionFilter->Update();
  movingSliceXY = maximumProjectionFilter->GetOutput();
  movingSliceXY->DisconnectPipeline();
  }

  // factor by which the slice in XZ or YZ plane will be resampled along Z direction
  double resamplingFactor = spacing[2] / spacing[0];

  // YZ plane
  //
  SliceImageType::Pointer  movingSliceYZ;
  {
  typename ExtractFilterType::Pointer  extractFilterYZ = ExtractFilterType::New();
  typename ExtractFilterType::InputImageRegionType  extractionRegion;

  extractFilterYZ->SetInput( Image );

  // using m_ROI -- since m_MovingMappedImage is in the fixed coordinate system
  typename ExtractFilterType::InputImageRegionType::SizeType  size = m_ROI.GetSize();
  typename ExtractFilterType::InputImageRegionType::IndexType  index = m_ROI.GetIndex();
  // add half size to get the slice in the middle of the ROI
  index[0] += size[0] / 2;
  // extracting X slice
  size[0] = 0;

  extractionRegion.SetSize( size );
  extractionRegion.SetIndex( index );
  bool canCrop = extractionRegion.Crop( Image->GetLargestPossibleRegion() );
  if( !canCrop || !Image->GetLargestPossibleRegion().IsInside( index ) ) {
    std::cout << "Error: point outside the image. Or the ROI outside of image." << std::endl;
    return;
  }

  extractFilterYZ->SetExtractionRegion( extractionRegion );
  extractFilterYZ->Update();

  SliceImageType::Pointer  extractedYZ = extractFilterYZ->GetOutput();
  extractedYZ->DisconnectPipeline();

  // The resolution in the Z direction is small, so resample the image
  typedef itk::ResampleImageFilter< SliceImageType, SliceImageType >  ResamplerType;
  ResamplerType::Pointer movingMappedResampler = ResamplerType::New();
  movingMappedResampler->SetInput( extractedYZ );

  movingMappedResampler->SetOutputOrigin( extractedYZ->GetOrigin() );
  movingMappedResampler->SetOutputDirection( extractedYZ->GetDirection() );

  // Trick: Need to change the size of the index as well
  typename SliceImageType::IndexType    movingMappedIndex   = extractedYZ->GetLargestPossibleRegion().GetIndex();
  movingMappedIndex[1] = movingMappedIndex[1] * resamplingFactor;
  movingMappedResampler->SetOutputStartIndex( movingMappedIndex);

  typename SliceImageType::SizeType     movingMappedSize    = extractedYZ->GetLargestPossibleRegion().GetSize();
  movingMappedSize[1] = movingMappedSize[1] * resamplingFactor;
  movingMappedResampler->SetSize( movingMappedSize );

  SliceImageType::SpacingType  movingMappedSpacing = extractedYZ->GetSpacing();
  movingMappedSpacing[1] = movingMappedSpacing[1] / resamplingFactor;
  movingMappedResampler->SetOutputSpacing( movingMappedSpacing);

  movingMappedResampler->Update();

  movingSliceYZ = movingMappedResampler->GetOutput();
  movingSliceYZ->DisconnectPipeline();
  }


  // XZ plane
  //
  SliceImageType::Pointer  movingSliceXZ;
  {
  typename ExtractFilterType::Pointer  extractFilterXZ = ExtractFilterType::New();
  typename ExtractFilterType::InputImageRegionType  extractionRegion;

  extractFilterXZ->SetInput( Image );

  // using m_ROI -- since m_MovingMappedImage is in the fixed coordinate system
  typename ExtractFilterType::InputImageRegionType::SizeType  size = m_ROI.GetSize();
  typename ExtractFilterType::InputImageRegionType::IndexType  index = m_ROI.GetIndex();
  // add half size to get the slice in the middle of the ROI
  index[1] += size[1] / 2;
  // extracting Y slice
  size[1] = 0;

  extractionRegion.SetSize( size );
  extractionRegion.SetIndex( index );
  bool canCrop = extractionRegion.Crop( Image->GetLargestPossibleRegion() );
  if( !canCrop || !Image->GetLargestPossibleRegion().IsInside( index ) ) {
    std::cout << "Error: point outside the image. Or the ROI outside of image." << std::endl;
    return;
  }

  extractFilterXZ->SetExtractionRegion( extractionRegion );
  extractFilterXZ->Update();

  SliceImageType::Pointer  extractedXZ = extractFilterXZ->GetOutput();
  extractedXZ->DisconnectPipeline();

  // The resolution in the Z direction is small, so resample the image
  typedef itk::ResampleImageFilter< SliceImageType, SliceImageType >  ResamplerType;
  ResamplerType::Pointer movingMappedResampler = ResamplerType::New();
  movingMappedResampler->SetInput( extractedXZ );

  movingMappedResampler->SetOutputOrigin( extractedXZ->GetOrigin() );
  movingMappedResampler->SetOutputDirection( extractedXZ->GetDirection() );

  // Trick: Need to change the size of the index as well
  typename SliceImageType::IndexType    movingMappedIndex   = extractedXZ->GetLargestPossibleRegion().GetIndex();
  movingMappedIndex[1] = movingMappedIndex[1] * resamplingFactor;
  movingMappedResampler->SetOutputStartIndex( movingMappedIndex);

  typename SliceImageType::SizeType     movingMappedSize    = extractedXZ->GetLargestPossibleRegion().GetSize();
  movingMappedSize[1] = movingMappedSize[1] * resamplingFactor;
  movingMappedResampler->SetSize( movingMappedSize );

  typename SliceImageType::SpacingType  movingMappedSpacing = extractedXZ->GetSpacing();
  movingMappedSpacing[1] = movingMappedSpacing[1] / resamplingFactor;
  movingMappedResampler->SetOutputSpacing( movingMappedSpacing);

  movingMappedResampler->Update();

  movingSliceXZ = movingMappedResampler->GetOutput();
  movingSliceXZ->DisconnectPipeline();
  }


  // the trick is to tile 2d images into a 3d volume, so that it can be passed into vtk display in the base class
  typedef itk::TileImageFilter< SliceImageType, SliceImageType >  TileFilterType;
  //typedef itk::TileImageFilter< typename cdcl_extract_data_callback3d< dim, dof >::OutputImageType, typename cdcl_extract_data_callback3d< dim, dof >::OutputImageType >  TileFilterType;
  //typedef itk::TileImageFilter< typename cdcl_extract_data_callback3d< dim, dof >::OutputImageType, SliceImageType >  TileFilterType;
  TileFilterType::Pointer  tileFilter = TileFilterType::New();
  tileFilter->SetInput( 0, movingSliceXY );
  tileFilter->SetInput( 1, movingSliceYZ );
  tileFilter->SetInput( 2, movingSliceXZ );
  TileFilterType::LayoutArrayType  layout;
  layout[0] = 3;
  layout[1] = 1;
  tileFilter->SetLayout( layout );
  tileFilter->Update();


  SliceImageType::Pointer                                     m_SlicesImage;
  m_SlicesImage = tileFilter->GetOutput();



  // Write the checkerboard
  //
  typedef SliceImageType                  WriteImageType;
  typedef itk::Image<vxl_byte, 2>                                            ByteImageType;
  typedef itk::ImageFileWriter<ByteImageType>                               WriterType;

  typedef itk::RescaleIntensityImageFilter<WriteImageType, ByteImageType>  RescalerType;
  //typedef itk::RescaleIntensityImageFilter<SliceImageType, WriteImageType>  RescalerType;
  //  Rescaler to set the final image to 0..255
  RescalerType::Pointer rescaler = RescalerType::New();
  rescaler->SetOutputMinimum( 0 );
  rescaler->SetOutputMaximum( 255 );
  rescaler->SetInput( m_SlicesImage );

  // need to set IO because the filename might include periods from decimals
  // and the factory does not know which writer to use (based on extension)
  typedef itk::PNGImageIO  ImageIOType;
  ImageIOType::Pointer  io = ImageIOType::New();

  //  Writer
  WriterType::Pointer writer = WriterType::New();
  writer->SetImageIO( io );
  writer->SetInput( rescaler->GetOutput() );

  writer->SetFileName( FileName.c_str() );
  writer->Update();

}

}


int main( int argc, char * argv[] )
{

  if( argc != 5 && argc != 6 && argc != 3 )
    {
    std::cerr << "Missing Parameters " << std::endl;
    std::cerr << "Usage: " << argv[0];
    std::cerr << " Image ";
    std::cerr << " PhysicalX PhysicalY PhysicalZ [prefix]" << std::endl;
    std::cerr << "OR: " << argv[0];
    std::cerr << " Image ";
    std::cerr << " LocationsFile" << std::endl;
    return 1;
    }

  const     unsigned int   ImageDimension = 3;

  typedef   short  InputPixelType;
  typedef   itk::Image< InputPixelType, ImageDimension >  ImageType;
  //typedef   itk::OrientedImage< InputPixelType, ImageDimension > ImageType; 

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
  //std::cout << "FIXED NAMES" << std::endl << names;
  std::cout << "Reading volume..." << std::endl;

  //reader->ReverseOrderOn();


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

  ImageType::Pointer  image = reader->GetOutput();


   typedef itk::MetaDataDictionary   DictionaryType;

   const  DictionaryType & dictionary = io->GetMetaDataDictionary();


   // Since we are interested only in the DICOM tags that can be expressed in
   // strings, we declare a MetaDataObject suitable for managing strings.
  typedef itk::MetaDataObject< std::string > MetaDataStringType;


  itk::Point< double, 3 >  patientPosition;

  // We instantiate the iterators that will make possible to walk through all the
  // entries of the MetaDataDictionary.
  DictionaryType::ConstIterator itr = dictionary.Begin();
  DictionaryType::ConstIterator end = dictionary.End();

  // For each one of the entries in the dictionary, we check first if its element
  // can be converted to a string, a \code{dynamic\_cast} is used for this purpose.
  while( itr != end )
    {
    itk::MetaDataObjectBase::Pointer  entry = itr->second;

    MetaDataStringType::Pointer entryvalue = 
      dynamic_cast<MetaDataStringType *>( entry.GetPointer() ) ;

    // For those entries that can be converted, we take their DICOM tag and pass it
    // to the \code{GetLabelFromTag()} method of the GDCMImageIO class. This method
    // checks the DICOM dictionary and returns the string label associated to the
    // tag that we are providing in the \code{tagkey} variable. If the label is
    // found, it is returned in \code{labelId} variable. The method itself return
    // false if the tagkey is not found in the dictionary.  For example "0010|0010"
    // in \code{tagkey} becomes "Patient's Name" in \code{labelId}.
    if( entryvalue )
      {
      std::string tagkey   = itr->first;
      std::string labelId;
      bool found =  itk::GDCMImageIO::GetLabelFromTag( tagkey, labelId );

      // The actual value of the dictionary entry is obtained as a string with the
      // \code{GetMetaDataObjectValue()} method.
      std::string tagvalue = entryvalue->GetMetaDataObjectValue();

      // At this point we can print out an entry by concatenating the DICOM Name or
      // label, the numeric tag and its actual value.
      if( found )
        {
        std::cout << "(" << tagkey << ") " << labelId;
        std::cout << " = " << tagvalue.c_str() << std::endl;

        if( tagkey == "0020|0032" ) // Patient Position
          {
            //std::string::iterator  backslash = find( tagkey.begin(), tagkey.end(), '\\' );
            //std::cout << "pos" << tagvalue.substr( 0, backslash - tagkey.begin() ) << std::endl;

            // Note: this could be done by ExposeMetaData from itk::MetaDataObject

            // the value is separated by slashes, so we replace them by spaces
            unsigned int first_backslash = tagvalue.find( '\\' );
            unsigned int second_backslash = tagvalue.find( '\\', first_backslash + 1 );

            std::string  spacesNotSlashes = tagvalue.replace( first_backslash, 1, " " );
            spacesNotSlashes = spacesNotSlashes.replace( second_backslash, 1, " " );

            std::istringstream  nameStream;
            nameStream.str( spacesNotSlashes );
            nameStream >> patientPosition[0] >> patientPosition[1] >> patientPosition[2];

          }
        }
      else
        {
        std::cout << "(" << tagkey <<  ") " << "Unknown";
        std::cout << " = " << tagvalue.c_str() << std::endl;
        }
      }

    // Finally we just close the loop that will walk through all the Dictionary
    // entries.
    ++itr;
    }

  std::cout << "Extracted Image Position (Patient): " << patientPosition << std::endl;


  std::cout << image << std::endl;



  std::string  prefix = "";
  if( argc == 6 ) prefix = argv[5];


  std::vector< itk::Point< double, 3 > >  locations;
  // if nodule file specified, read locations from it otherwise look for locations on the command line
  std::ifstream  noduleFile;
  if( argc == 3 ) {
    noduleFile.open( argv[2] );
    if( !noduleFile ) {
      std::cout << "Specified nodule file cannot be opened: " << argv[2] << std::endl;
      return 1;
    }
    prefix = itksys::SystemTools::GetFilenameWithoutExtension( argv[2] );

    itk::Point< double, 3 >  nodulePoint;
    while( ( noduleFile >> nodulePoint[0] >> nodulePoint[1] >> nodulePoint[2] ) ) {
      locations.push_back( nodulePoint );
    }

  }
  else {
    itk::Point< double, 3 >  pointRAS;
    pointRAS[0] = atof( argv[2] );
    pointRAS[1] = atof( argv[3] );
    pointRAS[2] = atof( argv[4] );

    locations.push_back( pointRAS );
  }



  for( std::vector< itk::Point< double, 3 > >::size_type  n = 0; n < locations.size(); ++n ) {
    itk::Point< double, 3 >  const &  pointRAS = locations[n];

    ImageType::SpacingType  spacing = image->GetSpacing();

    // This is not RAS but rather RAI
    // RAI (indicating Right to Left, Anterior to Posterior, Inferior to Superior) has the equivalent direction cosines [ 1, 0, 0 ], [ 0, 1, 0 ], [ 0, 0, 1].

    itk::ContinuousIndex< double, ImageDimension >  indexXYZ;

    #define RIDER 1

    #if RIDER

    indexXYZ[0] = pointRAS[0];
    indexXYZ[1] = pointRAS[1];
    indexXYZ[2] = pointRAS[2];

    #else   // LungCancerCases

    indexXYZ[0] = ( pointRAS[0] - patientPosition[0] ) / spacing[0];
    indexXYZ[1] = ( pointRAS[1] - patientPosition[1] ) / spacing[1];
    //indexXYZ[2] = ( patientPosition[2] - pointRAS[2] ) / spacing[2];
    indexXYZ[2] = ( patientPosition[2] - pointRAS[2] ) / spacing[2];

    // for some reason the slice number is from the last one
    indexXYZ[2] = image->GetLargestPossibleRegion().GetSize()[2] - indexXYZ[2];
    #endif


    std::cout << "Image Size: " << image->GetLargestPossibleRegion().GetSize() << std::endl;
    std::cout << "Image Origin: " << image->GetOrigin() << std::endl;


    itk::Point< double, 3 >  pointXYZ;
    image->TransformContinuousIndexToPhysicalPoint( indexXYZ, pointXYZ );

    std::cout << "RAI point: " << pointRAS << std::endl;
    std::cout << "XYZ index: " << indexXYZ << std::endl;
    std::cout << "XYZ point: " << pointXYZ << std::endl;


    //// If the volumes were cropped an meta data tags not modified, need to change the origin manually.
    ////
    //double newOrigin[3] = { -235.4, -184.769, -362.59 };
    ////either do this
    ////image->SetOrigin( newOrigin );
    ////or do this
    //pointRAS[0] += ( image->GetOrigin()[0] - newOrigin[0] );
    //pointRAS[1] += ( image->GetOrigin()[1] - newOrigin[1]);
    //pointRAS[2] += ( image->GetOrigin()[2] - newOrigin[2] );

    std::ostringstream  sliceName;
    sliceName << argv[1] << prefix << "_" << int( pointRAS[0] ) << "_" << int( pointRAS[1] ) << "_" << int( pointRAS[2] ) << ".png";
    std::cout << "Writing image slice: " << sliceName.str().c_str() << std::endl;

    itk::WriteROIImage<ImageType>( image, pointRAS, sliceName.str() );
    //itk::WriteROIImage<ImageType>( image, pointXYZ, "slice.png" );   // this one is off probably due to rounding in TransformContinuousIndexToPhysicalPoint
    //itk::WriteROIImage<ImageType>( image, pointRAS, "slice2.png" );  // this is the correct one

    std::cout << pointXYZ[0] << " " << pointXYZ[1] << " " << pointXYZ[2] << std::endl;
  }

  return EXIT_SUCCESS;
}

