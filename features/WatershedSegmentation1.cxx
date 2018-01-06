/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: WatershedSegmentation1.cxx,v $
  Language:  C++
  Date:      $Date: 2008/05/19 19:03:03 $
  Version:   $Revision: 1.11 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifdef _MSC_VER
#pragma warning ( disable : 4786 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif

#if defined(_MSC_VER)
#pragma warning ( disable : 4996 )
#pragma warning ( disable : 4786 )
#endif

//  Software Guide : BeginCommandLineArgs
//  INPUTS: {VisibleWomanEyeSlice.png}
//  OUTPUTS: {WatershedSegmentation1Output1.png}
//  2 10 0 0.05 1
//  Software Guide : EndCommandLineArgs
//  Software Guide : BeginCommandLineArgs
//  INPUTS: {VisibleWomanEyeSlice.png}
//  OUTPUTS: {WatershedSegmentation1Output2.png}
//  2 10 0.001 0.15 0
//  Software Guide : EndCommandLineArgs



// Software Guide : BeginLatex
//
// The following example illustrates how to preprocess and segment images
// using the \doxygen{WatershedImageFilter}. Note that the care with which
// the data is preprocessed will greatly affect the quality of your result.
// Typically, the best results are obtained by preprocessing the original
// image with an edge-preserving diffusion filter, such as one of the
// anisotropic diffusion filters, or with the bilateral image filter.  As
// noted in Section~\ref{sec:AboutWatersheds}, the height function used as
// input should be created such that higher positive values correspond to
// object boundaries.  A suitable height function for many applications can
// be generated as the gradient magnitude of the image to be segmented.
//
// The \doxygen{VectorGradientMagnitudeAnisotropicDiffusionImageFilter} class
// is used to smooth the image and the
// \doxygen{VectorGradientMagnitudeImageFilter} is used to generate the
// height function.  We begin by including all preprocessing filter header
// files and the header file for the WatershedImageFilter.  We
// use the vector versions of these filters because the input data is a color
// image.
//
// 
// Software Guide : EndLatex
#include <iostream>

#include "itkImageSeriesReader.h"
#include "itkImageSeriesWriter.h"
#include "itkNumericSeriesFileNames.h"

#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"

// Software Guide : BeginCodeSnippet
#include "itkGradientAnisotropicDiffusionImageFilter.h"
#include "itkGradientMagnitudeImageFilter.h"
#include "itkWatershedImageFilter.h"
//#include "itkMorphologicalWatershedImageFilter.h"
#include "itkGrayscaleMorphologicalClosingImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkThresholdLabelerImageFilter.h"
#include "itkMultiResolutionPyramidImageFilter.h"

// Software Guide : EndCodeSnippet

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkVectorCastImageFilter.h"
#include "itkUnaryFunctorImageFilter.h"
#include "itkScalarToRGBPixelFunctor.h"



// WatershedSegmentation1.exe 000002whole 000002wholewatershed 2.0 10 0.001 0.10


int main( int argc, char *argv[] )
{
  if (argc < 7 )
    {
    std::cerr << "Missing Parameters " << std::endl;
    std::cerr << "Usage: " << argv[0];
    std::cerr << " inputImage outputImage conductanceTerm diffusionIterations lowerThreshold outputScaleLevel " << std::endl;
    return 1;
    }
  
  // Software Guide : BeginLatex
  //
  // We now declare the image and pixel types to use for instantiation of the
  // filters.  All of these filters expect real-valued pixel types in order to
  // work properly.  The preprocessing stages are done directly on the
  // vector-valued data and the segmentation is done using floating point
  // scalar data.  Images are converted from RGB pixel type to
  // numerical vector type using \doxygen{VectorCastImageFilter}.
  //
  // Software Guide : EndLatex



  const     unsigned int   ImageDimension = 3;

  typedef   float  InputPixelType;

  typedef   itk::Image< InputPixelType, ImageDimension >  ImageType;

  typedef   itk::ImageSeriesReader< ImageType >   ReaderType;

  typedef itk::RGBPixel<unsigned char>   RGBPixelType;

  typedef itk::Image<RGBPixelType, ImageDimension>    RGBImageType;

  typedef itk::Image<unsigned long, ImageDimension>   LabeledImageType;

  typedef   itk::Image< RGBPixelType, 2 >              SliceImageType;

  typedef   itk::ImageSeriesWriter< RGBImageType, SliceImageType >  WriterType;

  typedef   itk::ImageSeriesWriter< LabeledImageType, SliceImageType >  UnsignedLongWriterType;




  // Software Guide : BeginLatex
  //
  // The various image processing filters are declared using the types created
  // above and eventually used in the pipeline.
  //
  // Software Guide : EndLatex

  // Software Guide : BeginCodeSnippet
  typedef itk::GradientAnisotropicDiffusionImageFilter<ImageType,
    ImageType>  DiffusionFilterType;

  typedef itk::GradientMagnitudeImageFilter<ImageType,ImageType>
    GradientMagnitudeFilterType; 

  typedef itk::WatershedImageFilter<ImageType> WatershedFilterType;

  //typedef itk::MorphologicalWatershedImageFilter<ImageType,LabeledImageType> MorphologicalWatershedFilterType;
  // Software Guide : EndCodeSnippet




  //typedef itk::BinaryBallStructuringElement< float, 3 >  StructuringElementType;
  //typedef itk::GrayscaleMorphologicalClosingImageFilter<ImageType, ImageType, StructuringElementType>  ClosingFilterType;

  //ClosingFilterType::Pointer  closingFilter = ClosingFilterType::New();

  //StructuringElementType  element;
  //element.SetRadius( 5 );
  //element.CreateStructuringElement();

  //closingFilter->SetKernel( element );


  std::cout << "Reading image..." << std::endl;

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

  ImageType::Pointer readImage = imageReader->GetOutput();



  if( readImage.IsNull() ) {
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

    readImage = reader->GetOutput();
  }

 
  // Software Guide : BeginLatex
  //
  // Next we instantiate the filters and set their parameters.  The first
  // step in the image processing pipeline is diffusion of the color input
  // image using an anisotropic diffusion filter.  For this class of filters,
  // the CFL condition requires that the time step be no more than 0.25 for
  // two-dimensional images, and no more than 0.125 for three-dimensional
  // images.  The number of iterations and the conductance term will be taken
  // from the command line. See
  // Section~\ref{sec:EdgePreservingSmoothingFilters} for more information on
  // the ITK anisotropic diffusion filters.
  //
  // Software Guide : EndLatex
  
  // Software Guide : BeginCodeSnippet
  DiffusionFilterType::Pointer diffusion = DiffusionFilterType::New();
  diffusion->SetNumberOfIterations( atoi(argv[4]) );
  diffusion->SetConductanceParameter( atof(argv[3]) );
  diffusion->SetTimeStep(0.05);
  //diffusion->SetTimeStep(0.125);
  // Software Guide : EndCodeSnippet

  // Software Guide : BeginLatex
  //
  // The ITK gradient magnitude filter for vector-valued images can optionally
  // take several parameters.  Here we allow only enabling or disabling 
  // of principal component analysis.
  //
  // Software Guide : EndLatex

  // Software Guide : BeginCodeSnippet
  GradientMagnitudeFilterType::Pointer 
  gradient = GradientMagnitudeFilterType::New();
  // Software Guide : EndCodeSnippet


  // Software Guide : BeginLatex
  //
  // Finally we set up the watershed filter.  There are two parameters.
  // \code{Level} controls watershed depth, and \code{Threshold} controls the
  // lower thresholding of the input.  Both parameters are set as a
  // percentage (0.0 - 1.0) of the maximum depth in the input image.
  //
  // Software Guide : EndLatex

  // Software Guide : BeginCodeSnippet
  WatershedFilterType::Pointer watershed = WatershedFilterType::New();
  watershed->SetLevel( atof(argv[6]) );
  watershed->SetThreshold( atof(argv[5]) );





  //typedef itk::ThresholdLabelerImageFilter< ImageType, LabeledImageType >  ThresholdLabelerType;
  //ThresholdLabelerType::Pointer  thresholdLabeler = ThresholdLabelerType::New();

  //ThresholdLabelerType::ThresholdVector  thresholds;
  //thresholds.push_back( -1000 );
  //thresholds.push_back( -120 );
  //thresholds.push_back( 0 );
  //thresholds.push_back( 40 );
  //thresholds.push_back( 400 );

  //thresholdLabeler->SetThresholds( thresholds );





  //MorphologicalWatershedFilterType::Pointer watershed = MorphologicalWatershedFilterType::New();
  //watershed->SetLevel( atof(argv[6]) );
  // Software Guide : EndCodeSnippet

  // Software Guide : BeginLatex
  //
  // The output of WatershedImageFilter is an image of unsigned long integer
  // labels, where a label denotes membership of a pixel in a particular
  // segmented region.  This format is not practical for visualization, so
  // for the purposes of this example, we will convert it to RGB pixels.  RGB
  // images have the advantage that they can be saved as a simple png file
  // and viewed using any standard image viewer software.  The
  // \subdoxygen{Functor}{ScalarToRGBPixelFunctor} class is a special
  // function object designed to hash a scalar value into an
  // \doxygen{RGBPixel}. Plugging this functor into the
  // \doxygen{UnaryFunctorImageFilter} creates an image filter for that
  // converts scalar images to RGB images.
  //
  // Software Guide : EndLatex

  // Software Guide : BeginCodeSnippet
  typedef itk::Functor::ScalarToRGBPixelFunctor<unsigned long>
    ColorMapFunctorType;

  typedef itk::UnaryFunctorImageFilter<LabeledImageType,
    RGBImageType, ColorMapFunctorType> ColorMapFilterType;

  ColorMapFilterType::Pointer colormapper = ColorMapFilterType::New();
  // Software Guide : EndCodeSnippet


  std::string  watershedDir = std::string( argv[2] );

  itk::NumericSeriesFileNames::Pointer fit = itk::NumericSeriesFileNames::New();
  char format[4096];
  sprintf (format, "%s/image%%04d.dcm", watershedDir.c_str());

  std::cout << "Format = " << format << std::endl;

  fit->SetSeriesFormat( format );
  fit->SetStartIndex( 0 );
  fit->SetIncrementIndex( 1 );


  std::string  watershedColorDir = watershedDir + "color";

  itk::NumericSeriesFileNames::Pointer fitColor = itk::NumericSeriesFileNames::New();
  char formatColor[4096];
  sprintf (formatColor, "%s/image%%04d.dcm", watershedColorDir.c_str());

  std::cout << "Format = " << formatColor << std::endl;

  fitColor->SetSeriesFormat( formatColor );
  fitColor->SetStartIndex( 0 );
  fitColor->SetIncrementIndex( 1 );


  WriterType::Pointer writer = WriterType::New();

  itksys::SystemTools::MakeDirectory( watershedColorDir.c_str() );

  //writer->SetMetaDataDictionaryArray( 
  //                      reader->GetMetaDataDictionaryArray() );


  UnsignedLongWriterType::Pointer unsignedLongWriter = UnsignedLongWriterType::New();

  itksys::SystemTools::MakeDirectory( watershedDir.c_str() );

  //unsignedLongWriter->SetMetaDataDictionaryArray( 
  //                      reader->GetMetaDataDictionaryArray() );



  typedef   itk::Image< unsigned int, 3 >  UnsignedIntImageType;
  typedef   itk::CastImageFilter< LabeledImageType, UnsignedIntImageType >  CasterType;

  CasterType::Pointer  caster = CasterType::New();


  typedef   itk::ImageFileWriter< UnsignedIntImageType >  LabeledFileWriterType;
  LabeledFileWriterType::Pointer  labeledFileWriter = LabeledFileWriterType::New();
  std::string  fileName = watershedDir + ".mhd";
  labeledFileWriter->SetFileName( fileName );


  // output filenames the same as the input filenames
  //
  //names->SetOutputDirectory( watershedColorDir.c_str() );

  //WriterType::FileNamesContainer  outputFileNamesColor = 
  //                          names->GetOutputFileNames();
  //writer->SetFileNames( outputFileNamesColor );


  //names->SetOutputDirectory( watershedDir.c_str() );

  //WriterType::FileNamesContainer  outputFileNames = 
  //                          names->GetOutputFileNames();

  //unsignedLongWriter->SetFileNames( outputFileNames );




  typedef itk::MultiResolutionPyramidImageFilter< ImageType, ImageType >  MultiResolutionFilterType;
  MultiResolutionFilterType::Pointer  multiResolutionFilter = MultiResolutionFilterType::New();

  unsigned int numberOfLevels = 2;
  multiResolutionFilter->SetInput( readImage );
  multiResolutionFilter->SetNumberOfLevels( numberOfLevels );
  multiResolutionFilter->Update();

  unsigned int level = 0;
  std::cout << "Processing level " << level << "..." << std::endl;

  MultiResolutionFilterType::OutputImagePointer  inputImage = multiResolutionFilter->GetOutput( level );

  //ImageType::Pointer inputImage = reader->GetOutput();



  //writer->SetFileName(argv[2]);

  // Software Guide : BeginLatex
  //
  // The filters are connected into a single pipeline, with readers and
  // writers at each end.
  //
  // Software Guide : EndLatex

  //  Software Guide : BeginCodeSnippet
  diffusion->SetInput(inputImage);

  //closingFilter->SetInput(diffusion->GetOutput());
  //gradient->SetInput(closingFilter->GetOutput());

  //thresholdLabeler->SetInput(diffusion->GetOutput());
  //colormapper->SetInput(thresholdLabeler->GetOutput());

  gradient->SetInput(diffusion->GetOutput());
  watershed->SetInput(gradient->GetOutput());
  colormapper->SetInput(watershed->GetOutput());

  writer->SetInput(colormapper->GetOutput());
  unsignedLongWriter->SetInput(watershed->GetOutput());

  caster->SetInput(watershed->GetOutput());
  labeledFileWriter->SetInput(caster->GetOutput());

  //unsignedLongWriter->SetInput(thresholdLabeler->GetOutput());

  // Software Guide : EndCodeSnippet

  try 
    {
    std::cout << "Diffusion..." << std::endl;
    diffusion->Update();
    std::cout << "Gradient..." << std::endl;
    gradient->Update();
    std::cout << "Watershed..." << std::endl;
    watershed->Update();
    std::cout << "Colormap..." << std::endl;
    colormapper->Update();


    //for( unsigned int i = 0; i < 30; ++i ) {
    //  LabeledImageType::IndexType index = {30, 30, 30};
    //  index[2] = i;
    //  std::cout << watershed->GetOutput()->GetPixel( index ) << std::endl;
    //}


    ImageType::SizeType  watershedSize = watershed->GetOutput()->GetLargestPossibleRegion().GetSize();
    fit->SetEndIndex( watershedSize[2]-1 );  // The number of slices to write
    fitColor->SetEndIndex( watershedSize[2]-1 );  // The number of slices to write

    unsignedLongWriter->SetFileNames( fit->GetFileNames() );
    writer->SetFileNames( fitColor->GetFileNames() );


    std::cout << "Writing image..." << std::endl;

    writer->Update();
    unsignedLongWriter->Update();

    labeledFileWriter->Update();
    }
  catch (itk::ExceptionObject &e)
    {
    std::cerr << e << std::endl;
    }
    
  return 0;
}

//
// Software Guide : BeginLatex
//
// \begin{figure} \center
// \includegraphics[width=0.32\textwidth]{VisibleWomanEyeSlice.eps}
// \includegraphics[width=0.32\textwidth]{WatershedSegmentation1Output1.eps}
// \includegraphics[width=0.32\textwidth]{WatershedSegmentation1Output2.eps}
// \itkcaption[Watershed segmentation output]{Segmented section of Visible Human
// female head and neck cryosection data.  At left is the original image.  The
// image in the middle was generated with parameters: conductance = 2.0,
// iterations = 10, threshold = 0.0, level = 0.05, principal components = on.
// The image on the right was generated with parameters: conductance = 2.0,
// iterations = 10, threshold = 0.001, level = 0.15, principal components =
// off. } \label{fig:outputWatersheds} \end{figure}
//
//
// Tuning the filter parameters for any particular application is a process
// of trial and error.  The \emph{threshold} parameter can be used to great
// effect in controlling oversegmentation of the image.  Raising the
// threshold will generally reduce computation time and produce output with
// fewer and larger regions.  The trick in tuning parameters is to consider
// the scale level of the objects that you are trying to segment in the
// image.  The best time/quality trade-off will be achieved when the image is
// smoothed and thresholded to eliminate features just below the desired
// scale.
//
// Figure~\ref{fig:outputWatersheds} shows output from the example code.  The
// input image is taken from the Visible Human female data around the right
// eye.  The images on the right are colorized watershed segmentations with
// parameters set to capture objects such as the optic nerve and
// lateral rectus muscles, which can be seen just above and to the left and
// right of the eyeball.  Note that a critical difference between the two
// segmentations is the mode of the gradient magnitude calculation.
//
// A note on the computational complexity of the watershed algorithm is
// warranted.  Most of the complexity of the ITK implementation lies in
// generating the hierarchy. Processing times for this stage are non-linear
// with respect to the number of catchment basins in the initial segmentation.
// This means that the amount of information contained in an image is more
// significant than the number of pixels in the image.  A very large, but very
// flat input take less time to segment than a very small, but very detailed
// input.
// 
// Software Guide : EndLatex

