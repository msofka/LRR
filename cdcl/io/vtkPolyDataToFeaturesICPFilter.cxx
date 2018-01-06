/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkPolyDataToFeaturesICPFilter.cxx,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkPolyDataToFeaturesICPFilter.h"

#include "vtkDataSet.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkPointData.h"
#include "vtkFloatArray.h"

#include "vtkLocator.h"
#include "vtkExecutive.h"
#include "vtkInformation.h"
#include "vtkGarbageCollector.h"
#include "vtkObjectFactory.h"
#include "vtkSmartPointer.h"

#include <vnl/vnl_cross.h>

vtkCxxRevisionMacro(vtkPolyDataToFeaturesICPFilter, "$Revision: 1.1 $");
vtkStandardNewMacro(vtkPolyDataToFeaturesICPFilter);

vtkStandardNewMacro(vtkFeatureICPAttributeSet);
vtkCxxRevisionMacro(vtkFeatureICPAttributeSet, "$Revision: 1.1 $");

//----------------------------------------------------------------------------
vtkPolyDataToFeaturesICPFilter::vtkPolyDataToFeaturesICPFilter()
{
  this->NumberOfRequiredInputs = 1;
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
  this->NumberOfOutputs = 1;

  // Copied from vtkPolyDataReader constructor:
  AttributeSetType *output = AttributeSetType::New();
  this->SetOutput(output);
  // Releasing data for pipeline parallism.
  // Filters will know it is empty. 
  output->ReleaseData();
  output->Delete();

}

//----------------------------------------------------------------------------
vtkPolyDataToFeaturesICPFilter::~vtkPolyDataToFeaturesICPFilter()
{
}

//----------------------------------------------------------------------------
// Specify the input data or filter.
void vtkPolyDataToFeaturesICPFilter::SetInput(vtkPolyData *  input)
{
  this->vtkProcessObject::SetNthInput(0, input);
}

//----------------------------------------------------------------------------
void vtkPolyDataToFeaturesICPFilter::SetOutput(AttributeSetType *output)
{
  this->GetExecutive()->SetOutputData(0, output);
}

//----------------------------------------------------------------------------
void vtkPolyDataToFeaturesICPFilter::Execute()
{
  this->GenerateOutput();
}

//----------------------------------------------------------------------------
// Specify the input data or filter.
vtkPolyData *vtkPolyDataToFeaturesICPFilter::GetInput()
{
  if (this->NumberOfInputs < 1)
    {
    return NULL;
    }
  
  return (vtkPolyData *)(this->Inputs[0]);
}

////----------------------------------------------------------------------------
//// Specify the input data or filter.
//void vtkPolyDataToFeaturesICPFilter::SetInput(vtkPolyData const &  input)
//{
//  this->vtkProcessObject::SetNthInput(0, input);
//}
//
////----------------------------------------------------------------------------
//// Specify the input data or filter.
//AttributeSetType *vtkPolyDataToFeaturesICPFilter::GetInput()
//{
//  if (this->NumberOfInputs < 1)
//    {
//    return NULL;
//    }
//  
//  return (vtkPolyData *)(this->Inputs[0]);
//}

//----------------------------------------------------------------------------
vtkPolyDataToFeaturesICPFilter::AttributeSetType *vtkPolyDataToFeaturesICPFilter::GetOutput()
{
  if (this->NumberOfOutputs < 1)
    {
    return NULL;
    }
  return (AttributeSetType *)(this->Outputs[0]);
}


void vtkPolyDataToFeaturesICPFilter::Update()
{
  this->GenerateOutput();
}

//----------------------------------------------------------------------------
// The conversion and poly data generation
void vtkPolyDataToFeaturesICPFilter::GenerateOutput()
{
  //AttributeSetType const &  input = this->GetInput();
  //if (input.size() == 0)
  //  {
  //  return;
  //  }
  vtkPolyData * poly_data = this->GetInput();
  if (!poly_data)
    {
    return;
    }
  // FIX -- TIME IS NOT CORRECT, SO UPDATE IS NOT DONE
  //int inputModified = (input->GetMTime() > this->GetMTime() ? 1 : 0);
  int inputModified = 1;


  //
  // If input to filter is modified, have to update
  //
  if ( inputModified )
    {
    // get points and covariances
    vtkSmartPointer< vtkPoints >     vtkpoints       = poly_data->GetPoints();
    //vtkSmartPointer< vtkDataArray >  covariances     = poly_data->GetPointData()->GetTensors( "covariances" );
    vtkSmartPointer< vtkDataArray >  errorProjectors = poly_data->GetPointData()->GetTensors( "errorProjectors" );
    vtkSmartPointer< vtkDataArray >  strengths       = poly_data->GetPointData()->GetScalars( "strengths" );
    vtkSmartPointer< vtkDataArray >  shapes          = poly_data->GetPointData()->GetScalars( "shapes" );
    vtkSmartPointer< vtkDataArray >  normals         = poly_data->GetPointData()->GetScalars( "normals" );
    vtkSmartPointer< vtkDataArray >  binormals       = poly_data->GetPointData()->GetScalars( "binormals" );


    // store the data
    AttributeSetType*  attributes = this->GetOutput();
    FeatureSetType &  features = attributes->m_FeatureSet;
    features.clear();
    const unsigned int dim = 3;
    unsigned int numCorners = 0;
    unsigned int numTubes = 0;
    unsigned int numSheets = 0;
    for( vtkIdType  i = 0; i < vtkpoints->GetNumberOfPoints(); ++i ) {
      vnl_vector_fixed< double, dim >  location;
      vtkpoints->GetPoint( i, location.data_block() );
  
      // ideally covariances should be stored as a 6-dimensional array but this is not the case for PolyData tensors
      // which are stored as a full 9-dimensional array
      //double *  covarianceArray = covariances->GetTuple( i );
      //vnl_matrix_fixed< double, dim, dim >  covariance( covarianceArray );

      double *  errorProjectorArray = errorProjectors->GetTuple( i );
      vnl_matrix_fixed< double, dim, dim >  errorProjector( errorProjectorArray );

      double strength = 0.0;
      if( strengths ) strengths->GetTuple( i, &strength );

      double shapeDouble = 0;
      if( shapes ) shapes->GetTuple( i, &shapeDouble );
      unsigned char shapeChar = (unsigned char) ( shapeDouble );
      
      FeatureType::shape_type  shape;
      vcl_vector< vnl_vector_fixed< double, dim > >  directions;
      switch( shapeChar ) {
        case 0: {
                shape = FeatureType::CORNER;
                double *  normalArray = normals->GetTuple( i );
                vnl_vector_fixed< double, dim >  normal( normalArray );
                double *  binormalArray = binormals->GetTuple( i );
                vnl_vector_fixed< double, dim >  binormal( binormalArray );
                vnl_vector_fixed< double, dim >  tangent = vnl_cross_3d( normal, binormal );
                directions.push_back( tangent );
                directions.push_back( binormal );
                directions.push_back( normal );
                ++numCorners;
      	        break;
                }
        case 1: {
                shape = FeatureType::TUBE;
                double *  normalArray = normals->GetTuple( i );
                vnl_vector_fixed< double, dim >  normal( normalArray );
                double *  binormalArray = binormals->GetTuple( i );
                vnl_vector_fixed< double, dim >  binormal( binormalArray );
                directions.push_back( binormal );
                directions.push_back( normal );
                ++numTubes;
      	        break;
                }
        case 2: {
                shape = FeatureType::SHEET;
                double *  normalArray = normals->GetTuple( i );
                vnl_vector_fixed< double, dim >  normal( normalArray );
                directions.push_back( normal );
                //std::cout << "n: " << normal << std::endl;
                ++numSheets;
      	        break;      	        
                }
      }
  
      #if HAS_DIRECTIONS
      FeaturePointer  feature = new FeatureType( location, strength, errorProjector, shape, directions );
      #else
      FeaturePointer  feature = new FeatureType( location, strength, errorProjector, shape );
      #endif
      features.push_back( feature );
    }

    std::cout << "Read: " << numCorners << " corners, " << numTubes << " tubes, " << numSheets << " sheets." << std::endl;

    }
}

//----------------------------------------------------------------------------
// Copy the update information across
void vtkPolyDataToFeaturesICPFilter::ComputeInputUpdateExtents(vtkDataObject *output)
{
  vtkDataObject *input = this->GetInput();

  if (input == NULL)
    {
    return;
    }
  
  this->vtkSource::ComputeInputUpdateExtents(output);
  input->RequestExactExtentOn();
}

//----------------------------------------------------------------------------
int
vtkPolyDataToFeaturesICPFilter
::FillInputPortInformation(int port, vtkInformation* info)
{
  if(!this->Superclass::FillInputPortInformation(port, info))
    {
    return 0;
    }
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
  return 1;
}

//----------------------------------------------------------------------------
int vtkPolyDataToFeaturesICPFilter
::FillOutputPortInformation(int port, vtkInformation* info)
{
  if(!this->Superclass::FillOutputPortInformation(port, info))
    {
    return 0;
    }
  info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkFeatureAttributeSet");
  return 1;
}

//----------------------------------------------------------------------------
void vtkPolyDataToFeaturesICPFilter::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}
