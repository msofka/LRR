/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkITKDescriptorPointSetToPolyDataFilter.cxx,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkITKDescriptorPointSetToPolyDataFilter.h"


// \brief Convert ITK PointSet to VTK PolyData format for saving

// \author: Michal Sofka, Rensselaer Polytechnic Institute (RPI)
//          sofka at cs dot rpi dot edu
// \date: 10/25/2007


#include "vtkDataSet.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkPointData.h"
#include "vtkFloatArray.h"

#include "vtkLocator.h"
#include "vtkInformation.h"
#include "vtkGarbageCollector.h"
#include "vtkObjectFactory.h"
#include "vtkSmartPointer.h"

vtkCxxRevisionMacro(vtkITKDescriptorPointSetToPolyDataFilter, "$Revision: 1.3 $");
vtkStandardNewMacro(vtkITKDescriptorPointSetToPolyDataFilter);

vtkStandardNewMacro(vtkDescriptorAttributeSet);
vtkCxxRevisionMacro(vtkDescriptorAttributeSet, "$Revision: 1.3 $");

//----------------------------------------------------------------------------
vtkITKDescriptorPointSetToPolyDataFilter::vtkITKDescriptorPointSetToPolyDataFilter()
{
  this->NumberOfRequiredInputs = 1;
  this->SetNumberOfInputPorts(1);
}

//----------------------------------------------------------------------------
vtkITKDescriptorPointSetToPolyDataFilter::~vtkITKDescriptorPointSetToPolyDataFilter()
{
}

//----------------------------------------------------------------------------
// Specify the input data or filter.
void vtkITKDescriptorPointSetToPolyDataFilter::SetInput(FeatureSetType *  input)
{
  this->vtkProcessObject::SetNthInput(0, input);
}

//----------------------------------------------------------------------------
void vtkITKDescriptorPointSetToPolyDataFilter::Execute()
{
  this->GenerateOutput();
}

//----------------------------------------------------------------------------
// Specify the input data or filter.
vtkDescriptorAttributeSet *vtkITKDescriptorPointSetToPolyDataFilter::GetInput()
{
  if (this->NumberOfInputs < 1)
    {
    return NULL;
    }
  
  return (vtkDescriptorAttributeSet *)(this->Inputs[0]);
}

////----------------------------------------------------------------------------
//// Specify the input data or filter.
//void vtkITKDescriptorPointSetToPolyDataFilter::SetInput(FeatureSetType const &  input)
//{
//  this->vtkProcessObject::SetNthInput(0, input);
//}
//
////----------------------------------------------------------------------------
//// Specify the input data or filter.
//FeatureSetType *vtkITKDescriptorPointSetToPolyDataFilter::GetInput()
//{
//  if (this->NumberOfInputs < 1)
//    {
//    return NULL;
//    }
//  
//  return (FeatureSetType *)(this->Inputs[0]);
//}

//----------------------------------------------------------------------------
vtkPolyData *vtkITKDescriptorPointSetToPolyDataFilter::GetOutput()
{
  if (this->NumberOfOutputs < 1)
    {
    return NULL;
    }
  
  return (vtkPolyData *)(this->Outputs[0]);
}

//----------------------------------------------------------------------------
vtkPolyData *vtkITKDescriptorPointSetToPolyDataFilter::GetOutput(int idx)
{
  return (vtkPolyData *) (this->Outputs[idx]);
}



//----------------------------------------------------------------------------
// The conversion and poly data generation
void vtkITKDescriptorPointSetToPolyDataFilter::GenerateOutput()
{
  //FeatureSetType const &  input = this->GetInput();
  //if (input.size() == 0)
  //  {
  //  return;
  //  }
  FeatureSetType * input = this->GetInput();
  if (!input)
    {
    return;
    }
  vtkPolyData *poly_data = this->GetOutput();
  // FIX -- TIME IS NOT CORRECT, SO UPDATE IS NOT DONE
  //int inputModified = (input->GetMTime() > this->GetMTime() ? 1 : 0);
  int inputModified = 1;


  //
  // If input to filter is modified, have to update
  //
  if ( inputModified )
    {


    PointSetType::PointDataContainerPointer   pointData = input->GetPointData();
    PointSetType::PointsContainerPointer      points = input->GetPoints();

    if( points->size() == 0 ) return;

    vtkSmartPointer< vtkPoints >      vtkpoints      = vtkSmartPointer< vtkPoints >::New();
    //vtkSmartPointer< vtkFloatArray >  covariances    = vtkSmartPointer< vtkFloatArray >::New();
    //vtkSmartPointer< vtkFloatArray >  strengths      = vtkSmartPointer< vtkFloatArray >::New();
    vtkSmartPointer< vtkFloatArray >  normals        = vtkSmartPointer< vtkFloatArray >::New();
    vtkSmartPointer< vtkFloatArray >  binormals      = vtkSmartPointer< vtkFloatArray >::New();
    vtkSmartPointer< vtkFloatArray >  descriptorData = vtkSmartPointer< vtkFloatArray >::New();

    descriptorData->SetName( "descriptors" );
    descriptorData->SetNumberOfComponents( (*pointData)[0].m_Descriptor.size() );

    normals->SetName( "normals" );
    normals->SetNumberOfComponents( (*pointData)[0].m_Direction.GetVectorDimension() );

    binormals->SetName( "binormals" );
    binormals->SetNumberOfComponents( (*pointData)[0].m_Bidirection.GetVectorDimension() );
    
    //covariances->SetName( "covariances" );
    
    // ideally covariances should be stored as a 6-dimensional array but this is not the case for PolyData tensors
    // which are stored as a full 9-dimensional array
    //unsigned int cov_dof = PointSetType::PixelType::SymmetricMatrixType::InternalDimension;

    //covariances->SetNumberOfComponents( dim*dim );

    //strengths->SetName( "strengths" );
    //strengths->SetNumberOfComponents( 1 );

    // write data points
    //
    // size type declared private FeatureSetType::PointsContainer::size_type -> use unsigned int
    for( unsigned int  i = 0; i < points->size(); ++i ) {

      PointSetType::PointType::ValueType*  point_array = (*points)[i].GetDataPointer();
      vtkpoints->InsertNextPoint( point_array );

      //vnl_matrix_fixed< float, dim, dim >  covariance;
      //for( unsigned int d1 = 0; d1 < dim; ++d1 ) {
      //  for( unsigned int d2 = 0; d2 < dim; ++d2 ) {
      //    covariance( d1, d2 ) = (*pointData)[i].m_Covariance( d1, d2 );
      //  }
      //}
      ////std::cout << covariance << std::endl << std::endl;
      //covariances->InsertNextTupleValue( covariance.data_block() );

      //strengths->InsertNextValue( (*pointData)[i].m_Strength );

      normals->InsertNextTupleValue( (*pointData)[i].m_Direction.GetDataPointer() );
      
      binormals->InsertNextTupleValue( (*pointData)[i].m_Bidirection.GetDataPointer() );
   
      descriptorData->InsertNextTupleValue( (*pointData)[i].m_Descriptor.data_block() );

    }

    poly_data->SetPoints( vtkpoints );
    //poly_data->GetPointData()->SetActiveScalars( "strengths" );
    //poly_data->GetPointData()->SetScalars( strengths );
    poly_data->GetPointData()->SetActiveVectors( "normals" );
    poly_data->GetPointData()->SetVectors( normals );
    poly_data->GetPointData()->SetActiveVectors( "binormals" );
    poly_data->GetPointData()->SetVectors( binormals );
    //poly_data->GetPointData()->SetTensors( covariances );
    poly_data->GetPointData()->SetActiveScalars( "descriptors" );
    poly_data->GetPointData()->SetScalars( descriptorData );
    }
}

//----------------------------------------------------------------------------
// Copy the update information across
void vtkITKDescriptorPointSetToPolyDataFilter::ComputeInputUpdateExtents(vtkDataObject *output)
{
  vtkDataObject *input = this->GetInput();

  if (input == NULL)
    {
    return;
    }
  
  this->vtkPolyDataSource::ComputeInputUpdateExtents(output);
  input->RequestExactExtentOn();
}

//----------------------------------------------------------------------------
int
vtkITKDescriptorPointSetToPolyDataFilter
::FillInputPortInformation(int port, vtkInformation* info)
{
  if(!this->Superclass::FillInputPortInformation(port, info))
    {
    return 0;
    }
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDescriptorAttributeSet");
  return 1;
}

//----------------------------------------------------------------------------
int vtkITKDescriptorPointSetToPolyDataFilter
::FillOutputPortInformation(int port, vtkInformation* info)
{
  if(!this->Superclass::FillOutputPortInformation(port, info))
    {
    return 0;
    }
  info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
  return 1;
}

//----------------------------------------------------------------------------
void vtkITKDescriptorPointSetToPolyDataFilter::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}
