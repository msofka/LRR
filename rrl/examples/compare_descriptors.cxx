#include <vul/vul_file.h>
#include <vul/vul_arg.h>

#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include "itkBSplineDeformableTransform.h"
#include "itkTransformFileReader.h"
#include "itkAffineTransform.h"
#include "itkTransformFactory.h"
#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkMaximumProjectionImageFilter.h"
#include "itkResampleImageFilter.h"
#include "itkTileImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkBoundingBox.h"

#include "itkImageSeriesReader.h"
#include "itkImageSeriesWriter.h"

#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"

#include "itkKdTree.h"
#include "itkKdTreeGenerator.h"
#include "itkListSample.h"

#include <cdcl/cdcl_feature.h>
#include <cdcl/cdcl_utils_VTK.h>

#include <vnl/algo/vnl_symmetric_eigensystem.h>
#include <vnl/algo/vnl_generalized_eigensystem.h>
#include <vnl/vnl_random.h>

#include <vcl_utility.h>
#include <vcl_algorithm.h>


//:
// \file
// \brief  Comparing descriptors for two sets of keypoints -- testing changes.
// \author Michal Sofka
// \date   September 2007


// write debugging image chips for each descriptor
#define OUTPUT_IMAGES 0


namespace {

// Kd-tree for descriptor vectors
//static const unsigned int descriptorSize = 771;
//typedef itk::Vector< float, descriptorSize >                      DescriptorVectorType;
typedef itk::Array< float >                      DescriptorVectorType;
typedef itk::Statistics::ListSample< DescriptorVectorType >       DescriptorSampleType;
typedef itk::Statistics::KdTreeGenerator< DescriptorSampleType >  DescriptorTreeGeneratorType;
typedef DescriptorTreeGeneratorType::KdTreeType                   DescriptorTreeType;


// Function object for sorting vector of descriptor indices based on distance
class DescriptorAtIndexIsGreater {
public:
  DescriptorAtIndexIsGreater( DescriptorVectorType &  FixedDescriptor,
                              DescriptorTreeType::Pointer const &  MovingDescriptorTree  )
    : m_FixedDescriptor( FixedDescriptor ),
      m_MovingDescriptorTree( MovingDescriptorTree ) {}
  bool operator()( int i, int j )
    {
      DescriptorVectorType  movingDescriptor_i = m_MovingDescriptorTree->GetMeasurementVector( i );
      vnl_vector< float >  distancei = m_FixedDescriptor - movingDescriptor_i;
      DescriptorVectorType  movingDescriptor_j = m_MovingDescriptorTree->GetMeasurementVector( j );
      vnl_vector< float >  distancej = m_FixedDescriptor - movingDescriptor_j;
      return distancei.magnitude() < distancej.magnitude();
    }
private:
  DescriptorVectorType const &  m_FixedDescriptor;
  DescriptorTreeType::Pointer  const &  m_MovingDescriptorTree;
};



typedef std::pair< DescriptorTreeType::InstanceIdentifierVectorType, double >  distance_ratio_type;
bool DistanceRatioIsSmaller( distance_ratio_type const &  left,
                     distance_ratio_type const &  right )
{
  return left.second < right.second;
}

}


int
main( int argc, char* argv[] )
{
  // deformation field case:
  // suppose we have a deformation field which warps image1 to image2
  // then warper in itk warps image image1 (output) to image2 input, s.t. input = output + deformation
  // therefore movingDescriptorFile is the one of image2 and fixedDescriptorFile is the one of image1
  //
  // e.g. indexing_shape_context3dITK.exe 241/4804/000002croppeddesc0.vtk 241/4802/000002croppeddesc0.vtk -trans 241/4802_000002cropped-4804_000002cropped.vtk -mov 241/4804/000002cropped -fix 241/4802/000002cropped
  vul_arg< const char* >        movingDescriptorFile     ( 0, "VTK moving descriptor file (of a synthetically deformed image)" );
  vul_arg< const char* >        fixedDescriptorFile      ( 0, "VTK fixed descriptor file" );
  vul_arg< const char* >        transformFile            ( "-trans", "ITK transform file (transform used to warp points)", 0 );

  vul_arg_parse( argc, argv );

  vcl_cout << "Command: ";
  for( int i = 0; i < argc; ++i )
    vcl_cout << argv[i] << " ";
  vcl_cout << vcl_endl;



  // read moving descriptors
  vcl_vector< cdcl_keypoint< 3 >::sptr >                movingKeypoints;
  vcl_vector< vnl_vector< float > >                    movingDescriptors;
  cdcl_read_keypoint_descriptors_VTK( movingKeypoints, movingDescriptors, movingDescriptorFile() );


  // read fixed descriptors
  //
  vcl_vector< cdcl_keypoint< 3 >::sptr >                fixedKeypoints;
  vcl_vector< vnl_vector< float > >                    fixedDescriptors;
  cdcl_read_keypoint_descriptors_VTK( fixedKeypoints, fixedDescriptors, fixedDescriptorFile() );


  // setup transforms

  typedef itk::AffineTransform<double,3> AffineTransformType;

  const unsigned int ImageDimension = 3;
	const unsigned int SpaceDimension = ImageDimension;
	const unsigned int SplineOrder = 3;
	typedef double CoordinateRepType;

	typedef itk::BSplineDeformableTransform<
													CoordinateRepType,
													SpaceDimension,
													SplineOrder >     BSplineTransformType;

  BSplineTransformType::Pointer  bsplineTransform;

  // Deformation field
  typedef itk::Vector< float, ImageDimension >           VectorPixelType;
  typedef itk::Image<  VectorPixelType, ImageDimension > DeformationFieldType;
  DeformationFieldType::Pointer  deformationField;

  enum deformationTypeEnum { NONE, BSPLINE, FIELD };
  deformationTypeEnum  deformationType = NONE;

  if( transformFile.set() ) {
    vcl_string  transformExtension = itksys::SystemTools::GetFilenameLastExtension( transformFile() );
    //std::cout << "Transform extension: " << transformExtension << std::endl;
    // if warp field supplied
    if( transformExtension == ".vtk" || transformExtension == ".mhd" ) {
      deformationType = FIELD;

      typedef   itk::ImageFileReader< DeformationFieldType >  FieldReaderType;

      std::cout << "Reading deformation field: " << transformFile();

      FieldReaderType::Pointer fieldReader = FieldReaderType::New();
      fieldReader->SetFileName( transformFile() );
      fieldReader->Update();
      deformationField = fieldReader->GetOutput();

      std::cout << " of size " << deformationField->GetLargestPossibleRegion().GetSize() << std::endl;
    }
    else { // assume BSpline transform supplied
      deformationType = BSPLINE;

      // In order to read a transform file, we instantiate a TransformFileReader. 
      // Like the writer, the reader is not templated.
      itk::TransformFileReader::Pointer reader;
      reader = itk::TransformFileReader::New();

      // Some transforms (like the BSpline transform) might not be registered 
      // with the factory so we add them manually. 
      itk::TransformFactory<BSplineTransformType>::RegisterTransform();

      // We then set the name of the file we want to read, and call the
      // Update() function.
      reader->SetFileName( transformFile() );
    
      try
        {
        reader->Update();
        }
      catch( itk::ExceptionObject & excp )
        {
        std::cerr << "Error while reading the transform file" << std::endl;
        std::cerr << excp << std::endl;
        std::cerr << "[FAILED]" << std::endl;
        return EXIT_FAILURE;
        }

      // The transform reader is not template and therefore it retunrs a list
      // of \doxygen{Transform}. However, the reader instantiate the appropriate
      // transform class when reading the file but it is up to the user to
      // do the approriate cast.
      // To get the output list of transform we use the GetTransformList() function.
      typedef itk::TransformFileReader::TransformListType * TransformListType;
      TransformListType transforms = reader->GetTransformList();
      std::cout << "Number of transforms = " << transforms->size() << std::endl;

      // We then use an STL iterator to go trought the list of transforms. We show here
      // how to do the proper casting of the resulting transform.
      for( itk::TransformFileReader::TransformListType::const_iterator it = transforms->begin(); it != transforms->end(); ++it ) {
        if(!strcmp((*it)->GetNameOfClass(),"BSplineDeformableTransform"))
          {
          bsplineTransform = static_cast<BSplineTransformType*>((*it).GetPointer());
          bsplineTransform->Print(std::cout);
          }
        else
          {
            std::cout << "Do not know how to read: " << (*it)->GetNameOfClass() << std::endl;
            return -1;
          }
      }
    }
  }
  else { // set to identity
    std::cout << "No transform specified, errors will be printed without transforming keypoints." << std::endl;
    bsplineTransform = BSplineTransformType::New();

    BSplineTransformType::RegionType  region;
    BSplineTransformType::SizeType  size;
    size.Fill(10);
    region.SetSize(size);
    bsplineTransform->SetGridRegion( region );
    BSplineTransformType::OriginType  origin;
    origin.Fill ( 100 );
    bsplineTransform->SetGridOrigin ( origin );
    BSplineTransformType::SpacingType  spacing;
    spacing.Fill ( 1.5 );
    bsplineTransform->SetGridSpacing ( spacing );

    BSplineTransformType::ParametersType  parameters( bsplineTransform->GetNumberOfParameters() );
    bsplineTransform->SetParameters( parameters );

    bsplineTransform->SetIdentity();
    deformationType = BSPLINE;
  }


  if( movingDescriptors.size() != fixedDescriptors.size() ) {
    std::cout << "Error: The sizes are different: " << movingDescriptors.size() << " " << fixedDescriptors.size() << std::endl;
    return -1;
  }

  std::cout << "Size of descriptor: " << movingDescriptors[0].size() << std::endl;


  unsigned int numDifferent = 0;
  double difference = 0.0;
  double maxDifference = 0.0;
  double avgMagnitudeMoving = 0.0;
  for( vcl_vector< vnl_vector< float > >::size_type  n = 0; n < movingDescriptors.size(); ++n ) {
    if( movingDescriptors[n] != fixedDescriptors[n] ) {
      ++numDifferent;
      double currDifference = ( movingDescriptors[n] - fixedDescriptors[n] ).magnitude();
      difference += currDifference;
      if( currDifference > maxDifference ) maxDifference = currDifference;

      std::cout << movingDescriptors[n] << std::endl;
      std::cerr << fixedDescriptors[n] << std::endl;
      //std::cout << std::endl;
    }
    avgMagnitudeMoving += movingDescriptors[n].magnitude();
  }
  avgMagnitudeMoving /= double( movingDescriptors.size() );

  std::cout << "Num different: " << numDifferent << std::endl;
  std::cout << "Total difference: " << difference << std::endl;
  std::cout << "Average difference: " << difference / double( numDifferent ) << std::endl;
  std::cout << "Max difference: " << maxDifference << std::endl;
  std::cout << "Average moving magnitude: " << avgMagnitudeMoving << std::endl;

  return 0;
}
