/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkITKPointSetToPolyDataFilter.cxx,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkITKPointSetToPolyDataFilter.h"


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

vtkCxxRevisionMacro(vtkITKPointSetToPolyDataFilter, "$Revision: 1.5 $");
vtkStandardNewMacro(vtkITKPointSetToPolyDataFilter);

vtkStandardNewMacro(vtkPointAttributeSet);
vtkCxxRevisionMacro(vtkPointAttributeSet, "$Revision: 1.5 $");

//----------------------------------------------------------------------------
vtkITKPointSetToPolyDataFilter::vtkITKPointSetToPolyDataFilter()
{
  this->NumberOfRequiredInputs = 1;
  this->SetNumberOfInputPorts(1);
}

//----------------------------------------------------------------------------
vtkITKPointSetToPolyDataFilter::~vtkITKPointSetToPolyDataFilter()
{
}

//----------------------------------------------------------------------------
// Specify the input data or filter.
void vtkITKPointSetToPolyDataFilter::SetInput(FeatureSetType *  input)
{
  this->vtkProcessObject::SetNthInput(0, input);
}

//----------------------------------------------------------------------------
void vtkITKPointSetToPolyDataFilter::Execute()
{
  this->GenerateOutput();
}

//----------------------------------------------------------------------------
// Specify the input data or filter.
vtkPointAttributeSet *vtkITKPointSetToPolyDataFilter::GetInput()
{
  if (this->NumberOfInputs < 1)
    {
    return NULL;
    }
  
  return (vtkPointAttributeSet *)(this->Inputs[0]);
}

////----------------------------------------------------------------------------
//// Specify the input data or filter.
//void vtkITKPointSetToPolyDataFilter::SetInput(FeatureSetType const &  input)
//{
//  this->vtkProcessObject::SetNthInput(0, input);
//}
//
////----------------------------------------------------------------------------
//// Specify the input data or filter.
//FeatureSetType *vtkITKPointSetToPolyDataFilter::GetInput()
//{
//  if (this->NumberOfInputs < 1)
//    {
//    return NULL;
//    }
//  
//  return (FeatureSetType *)(this->Inputs[0]);
//}

//----------------------------------------------------------------------------
vtkPolyData *vtkITKPointSetToPolyDataFilter::GetOutput()
{
  if (this->NumberOfOutputs < 1)
    {
    return NULL;
    }
  
  return (vtkPolyData *)(this->Outputs[0]);
}

//----------------------------------------------------------------------------
vtkPolyData *vtkITKPointSetToPolyDataFilter::GetOutput(int idx)
{
  return (vtkPolyData *) (this->Outputs[idx]);
}



//----------------------------------------------------------------------------
// The conversion and poly data generation
void vtkITKPointSetToPolyDataFilter::GenerateOutput()
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
    PointSetType::PointsContainerPointer   points = input->GetPoints();

    vtkSmartPointer< vtkPoints >      vtkpoints       = vtkSmartPointer< vtkPoints >::New();
    vtkSmartPointer< vtkFloatArray >  errorProjectors = vtkSmartPointer< vtkFloatArray >::New();
    vtkSmartPointer< vtkFloatArray >  strengths       = vtkSmartPointer< vtkFloatArray >::New();
    vtkSmartPointer< vtkUnsignedCharArray >    shapes = vtkSmartPointer< vtkUnsignedCharArray >::New();
    vtkSmartPointer< vtkFloatArray >  normals     = vtkSmartPointer< vtkFloatArray >::New();
    vtkSmartPointer< vtkFloatArray >  binormals     = vtkSmartPointer< vtkFloatArray >::New();

    const unsigned int dim = 3;

    normals->SetName( "normals" );
    normals->SetNumberOfComponents( dim );

    binormals->SetName( "binormals" );
    binormals->SetNumberOfComponents( dim );

    errorProjectors->SetName( "errorProjectors" );

    #if COMPUTE_COVARIANCES
    vtkSmartPointer< vtkFloatArray >  covariances     = vtkSmartPointer< vtkFloatArray >::New();
    covariances->SetName( "covariances" );

    // ideally covariances should be stored as a 6-dimensional array but this is not the case for PolyData tensors
    // which are stored as a full 9-dimensional array
    //unsigned int cov_dof = PointSetType::PixelType::SymmetricMatrixType::InternalDimension;
    covariances->SetNumberOfComponents( dim*dim );
    #endif

    errorProjectors->SetNumberOfComponents( dim*dim );

    strengths->SetName( "strengths" );
    strengths->SetNumberOfComponents( 1 );

    shapes->SetName( "shapes" );
    shapes->SetNumberOfComponents( 1 );

    // write data points
    //
    // size type declared private FeatureSetType::PointsContainer::size_type -> use unsigned int
    for( unsigned int  i = 0; i < points->size(); ++i ) {

      PointSetType::PointType::ValueType*  point_array = (*points)[i].GetDataPointer();
      vtkpoints->InsertNextPoint( point_array );

      #if COMPUTE_COVARIANCE
      vnl_matrix_fixed< float, dim, dim >  covariance;
      for( unsigned int d1 = 0; d1 < dim; ++d1 ) {
        for( unsigned int d2 = 0; d2 < dim; ++d2 ) {
          covariance( d1, d2 ) = (*pointData)[i].m_Covariance( d1, d2 );
        }
      }
      covariances->InsertNextTupleValue( covariance.data_block() );
      #endif

      vnl_matrix_fixed< float, dim, dim >  errorProjector;
      for( unsigned int d1 = 0; d1 < dim; ++d1 ) {
        for( unsigned int d2 = 0; d2 < dim; ++d2 ) {
          errorProjector( d1, d2 ) = (*pointData)[i].m_ErrorProjector( d1, d2 );
        }

      }
      errorProjectors->InsertNextTupleValue( errorProjector.data_block() );


      strengths->InsertNextValue( (*pointData)[i].m_Strength );

      unsigned char shapeChar;
      if( (*pointData)[i].m_Shape == itk::PointAttribute::Corner ) {
        normals->InsertTupleValue( i, (*pointData)[i].m_Directions[0].GetVnlVector().data_block() );  // here could perhaps be .GetDataPointer() to avoid conversion through vnl_vector
        binormals->InsertTupleValue( i, (*pointData)[i].m_Directions[1].GetVnlVector().data_block() );
        //std::cout << "POINTSETTOPOLYDATA: " << (*pointData)[i].m_Directions[0].GetVnlVector() << std::endl << (*pointData)[i].m_Directions[1].GetVnlVector() << std::endl << (*pointData)[i].m_Directions[2].GetVnlVector() << std::endl;
        shapeChar = 0;
      }
      else if( (*pointData)[i].m_Shape == itk::PointAttribute::Tube ) {
        normals->InsertTupleValue( i, (*pointData)[i].m_Directions[0].GetVnlVector().data_block() );
        binormals->InsertTupleValue( i, (*pointData)[i].m_Directions[1].GetVnlVector().data_block() );
        shapeChar = 1;
      }
      else if( (*pointData)[i].m_Shape == itk::PointAttribute::Sheet ) {
        normals->InsertTupleValue( i, (*pointData)[i].m_Directions[0].GetVnlVector().data_block() );
        // poly data format is not smart enough: *each* point needs to have an attribute
        vnl_vector_fixed< float, 3 >  emptyVector( 0.0 );
        binormals->InsertTupleValue( i, emptyVector.data_block() );
        shapeChar = 2;
      }
      else std::cout << "Error: Unknown feature type." << std::endl;
      shapes->InsertNextValue( shapeChar );

    }

    poly_data->SetPoints( vtkpoints );
    poly_data->GetPointData()->SetActiveScalars( "strengths" );
    poly_data->GetPointData()->SetScalars( strengths );
    poly_data->GetPointData()->SetActiveScalars( "shapes" );
    poly_data->GetPointData()->SetScalars( shapes );
    poly_data->GetPointData()->SetActiveVectors( "normals" );
    poly_data->GetPointData()->SetVectors( normals );
    poly_data->GetPointData()->SetActiveVectors( "binormals" );
    poly_data->GetPointData()->SetVectors( binormals );
    #if COMPUTE_COVARIANCES
    poly_data->GetPointData()->SetActiveTensors( "covariances" );
    poly_data->GetPointData()->SetTensors( covariances );
    #endif
    poly_data->GetPointData()->SetActiveTensors( "errorProjectors" );
    poly_data->GetPointData()->SetTensors( errorProjectors );

    }
}

//----------------------------------------------------------------------------
// Copy the update information across
void vtkITKPointSetToPolyDataFilter::ComputeInputUpdateExtents(vtkDataObject *output)
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
vtkITKPointSetToPolyDataFilter
::FillInputPortInformation(int port, vtkInformation* info)
{
  if(!this->Superclass::FillInputPortInformation(port, info))
    {
    return 0;
    }
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointAttributeSet");
  return 1;
}

//----------------------------------------------------------------------------
int vtkITKPointSetToPolyDataFilter
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
void vtkITKPointSetToPolyDataFilter::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}
