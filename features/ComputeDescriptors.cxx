#include "vtkXMLPolyDataWriter.h"
#include "vtkXMLPolyDataReader.h"
#include "vtkSmartPointer.h"
#include "vtkPolyData.h"

#include "vtkITKDescriptorPointSetToPolyDataFilter.h"
#include "itkDescriptorMeshFilter.h"
#include "itkFeatureImageFilter.h"
#include "itkVTKPolyDataToMeshFilter.h"

//:
// \file
// \brief  Computing shape context in 3D.
// \author Michal Sofka
// \date   June 2008



int
main( int argc, char* argv[] )
{
  if( argc < 4 )
    {
    std::cerr << "Missing Parameters " << std::endl;
    std::cerr << "Usage: " << argv[0];
    std::cerr << " KeypointFile ";
    std::cerr << " FeatureFile ";
    std::cerr << " DescriptorFile" << std::endl;
    return 1;
    }

  const     unsigned int   Dimension = 3;

  vtkSmartPointer< vtkXMLPolyDataReader >  featureReader = vtkSmartPointer< vtkXMLPolyDataReader >::New();
  std::string  featureFilename = std::string( argv[2] );
  featureReader->SetFileName( featureFilename.c_str() );
  featureReader->Update();

  vtkSmartPointer< vtkXMLPolyDataReader >  keypointReader = vtkSmartPointer< vtkXMLPolyDataReader >::New();
  keypointReader->SetFileName( argv[1] );
  keypointReader->Update();

  typedef itk::PointAttribute  PointAttributeType;
  typedef itk::PointSet< PointAttributeType,  Dimension >   PointSetType;

  // create the feature converting filter
  typedef itk::VTKPolyDataToMeshFilter< PointSetType >  FeatureConvertorType;
  FeatureConvertorType::Pointer  featureConverter = FeatureConvertorType::New();
  featureConverter->SetInput( featureReader->GetOutput() );

  // create the keypoint converting filter
  typedef itk::VTKPolyDataToMeshFilter< PointSetType >  KeypointConvertorType;
  KeypointConvertorType::Pointer  keypointConverter = KeypointConvertorType::New();
  keypointConverter->SetInput( keypointReader->GetOutput() );


  typedef itk::DescriptorAttribute  DescriptorAttributeType;
  typedef itk::PointSet< DescriptorAttributeType,  Dimension >   DescriptorSetType;

  itk::TimeProbe timer;
  timer.Start();

  typedef itk::DescriptorMeshFilter< PointSetType, PointSetType, DescriptorSetType >  DescriptorFilterType;
  DescriptorFilterType::Pointer  descriptorFilter = DescriptorFilterType::New();
  descriptorFilter->SetInputKeypoints( keypointConverter->GetOutput() );
  descriptorFilter->SetInputFeatures( featureConverter->GetOutput() );

  descriptorFilter->Update();
  DescriptorFilterType::OutputMeshPointer  descriptors = descriptorFilter->GetOutput();

  
  timer.Stop();
  std::cout << "Descriptor computation took " << timer.GetMeanTime() << " seconds.\n";


  // create ITK PointSet to VTK object connector
  vtkSmartPointer< vtkDescriptorAttributeSet >  descriptorAttrubuteSetVTK = vtkSmartPointer< vtkDescriptorAttributeSet >::New();
  descriptorAttrubuteSetVTK->SetInput( descriptors );

  // create the converting filter
  vtkSmartPointer< vtkITKDescriptorPointSetToPolyDataFilter >  descriptorSetToPolyDataFilter = vtkSmartPointer< vtkITKDescriptorPointSetToPolyDataFilter >::New();
  descriptorSetToPolyDataFilter->SetInput( descriptorAttrubuteSetVTK->GetOutput() );
  descriptorSetToPolyDataFilter->Update();

  // write the features
  vtkSmartPointer< vtkXMLPolyDataWriter >  descriptorWriter = vtkSmartPointer< vtkXMLPolyDataWriter >::New();

  vtkSmartPointer< vtkPolyData >  poly_data = descriptorSetToPolyDataFilter->GetOutput();

  descriptorWriter->SetInput( poly_data );
  descriptorWriter->SetFileName( argv[3] );
  descriptorWriter->Write();
  
  return 0;
}
