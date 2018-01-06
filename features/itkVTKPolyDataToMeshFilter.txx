/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkVTKPolyDataToMeshFilter.txx,v $
  Language:  C++
  Date:      $Date: 2008/07/20 17:20:24 $
  Version:   $Revision: 1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkVTKPolyDataToMeshFilter_txx
#define __itkVTKPolyDataToMeshFilter_txx

#include "itkVTKPolyDataToMeshFilter.h"
#include <fstream>
#include <stdio.h>

#include "vtkDataArray.h"
#include "vtkFloatArray.h"
#include "vtkPointData.h"

namespace itk
{


/**
 *
 */
//template <class TOutputMesh>
//void 
//VTKPolyDataToMeshFilter<TOutputMesh>
//::SetInput(PolyDataConnectorType *input)
//{
//  this->ProcessObject::SetNthInput(0, input);
//}
//
//
///**
// *
// */
//template <class TOutputMesh>
//typename VTKPolyDataToMeshFilter<TOutputMesh>::PolyDataConnectorType *
//VTKPolyDataToMeshFilter<TOutputMesh>
//::GetInput()
//{
//  if (this->GetNumberOfInputs() < 1)
//    {
//    return 0;
//    }
//  
//  return static_cast<PolyDataConnectorType*>
//    (this->ProcessObject::GetInput(0));
//}
//
//  
///**
// *
// */
//template <class TOutputMesh>
//typename VTKPolyDataToMeshFilter<TOutputMesh>::PolyDataConnectorType *
//VTKPolyDataToMeshFilter<TOutputMesh>
//::GetInput(unsigned int idx)
//{
//  return static_cast<PolyDataConnectorType*>
//    (this->ProcessObject::GetInput(idx));
//}
//


// 
// Constructor
// 
template<class TOutputMesh>
VTKPolyDataToMeshFilter<TOutputMesh>
::VTKPolyDataToMeshFilter()
{
  //
  // Create the output
  //
  typename TOutputMesh::Pointer outputMesh = TOutputMesh::New();
  this->ProcessObject::SetNumberOfRequiredOutputs(1);
  this->ProcessObject::SetNthOutput(0, outputMesh.GetPointer());

  // Modify superclass default values, can be overridden by subclasses
  //this->SetNumberOfRequiredInputs(1);

  //PointsContainerPointer     points     = PointsContainer::New();
  //PointDataContainerPointer  pointData  = PointDataContainer::New();

  //outputMesh->SetPoints( points );
  //outputMesh->SetPointData( pointData );

}

template<class TOutputMesh>
void
VTKPolyDataToMeshFilter<TOutputMesh>
::GenerateData()
{
  typename OutputMeshType::Pointer outputMesh = this->GetOutput();

  //outputMesh->SetCellsAllocationMethod(
  //    OutputMeshType::CellsAllocatedDynamicallyCellByCell );

  //PointsContainerPointer      points    = outputMesh->GetPoints();
  //PointDataContainerPointer   pointData = outputMesh->GetPointData();

  PointsContainerPointer     points     = PointsContainer::New();
  PointDataContainerPointer  pointData  = PointDataContainer::New();

  outputMesh->SetPoints( points );
  outputMesh->SetPointData( pointData );



  //InputPolyDataType *  poly_data = this->GetInput()->m_PolyData;
  InputPolyDataType *  poly_data = this->m_Input;

  // get points and covariances
  vtkSmartPointer< vtkPoints >     vtkpoints       = poly_data->GetPoints();
  //vtkSmartPointer< vtkFloatArray >   covariances = vtkFloatArray::SafeDownCast( poly_data->GetPointData()->GetTensors( "covariances" ) );
  vtkSmartPointer< vtkFloatArray >   errorProjectors = vtkFloatArray::SafeDownCast( poly_data->GetPointData()->GetTensors( "errorProjectors" ) );

  vtkSmartPointer< vtkFloatArray >   strengths = vtkFloatArray::SafeDownCast( poly_data->GetPointData()->GetScalars( "strengths" ) );
  vtkSmartPointer< vtkUnsignedCharArray >   shapes = vtkUnsignedCharArray::SafeDownCast( poly_data->GetPointData()->GetScalars( "shapes" ) );
  vtkSmartPointer< vtkFloatArray >   normals = vtkFloatArray::SafeDownCast( poly_data->GetPointData()->GetScalars( "normals" ) );
  vtkSmartPointer< vtkFloatArray >   binormals = vtkFloatArray::SafeDownCast( poly_data->GetPointData()->GetScalars( "binormals" ) );


  const unsigned int dim = 3;
  unsigned int numCorners = 0;
  unsigned int numTubes = 0;
  unsigned int numSheets = 0;
  for( vtkIdType  i = 0; i < vtkpoints->GetNumberOfPoints(); ++i ) {
    double *  location = vtkpoints->GetPoint( i );//, location.data_block() );

    // ideally covariances should be stored as a 6-dimensional array but this is not the case for PolyData tensors
    // which are stored as a full 9-dimensional array
    //double *  covarianceArray = covariances->GetTuple( i );
    //vnl_matrix_fixed< double, dim, dim >  covariance( covarianceArray );

    vnl_matrix_fixed< float, dim, dim >  errorProjector;
    errorProjectors->GetTupleValue( i, errorProjector.data_block() );

    float strength = 0.0f;
    if( strengths ) strengths->GetTupleValue( i, &strength );

    unsigned char shapeChar = 0;
    if( shapes ) shapes->GetTupleValue( i, &shapeChar );
    
    PointAttribute::FeatureShape  shape;
    vcl_vector< vnl_vector_fixed< float, dim > >  directions;
    switch( shapeChar ) {
      case 0: {
              shape = PointAttribute::Corner;
              vnl_vector_fixed< float, dim >  normal;
              normals->GetTupleValue( i, normal.data_block() );
              vnl_vector_fixed< float, dim >  binormal;
              binormals->GetTupleValue( i, binormal.data_block() );
              vnl_vector_fixed< float, dim >  tangent = vnl_cross_3d( normal, binormal );
              directions.push_back( normal );
              directions.push_back( binormal );
              directions.push_back( tangent );
              ++numCorners;
    	        break;
              }
      case 1: {
              shape = PointAttribute::Tube;
              vnl_vector_fixed< float, dim >  normal;
              normals->GetTupleValue( i, normal.data_block() );
              vnl_vector_fixed< float, dim >  binormal;
              binormals->GetTupleValue( i, binormal.data_block() );
              directions.push_back( normal );
              directions.push_back( binormal );
              ++numTubes;
    	        break;
              }
      case 2: {
              shape = PointAttribute::Sheet;
              vnl_vector_fixed< float, dim >  normal;
              normals->GetTupleValue( i, normal.data_block() );
              directions.push_back( normal );
              //std::cout << "n: " << normal << std::endl;
              ++numSheets;
    	        break;      	        
              }
    }

    PointAttribute::SymmetricMatrixType  errorProjectorITK;
    for( unsigned int d1 = 0; d1 < dim; ++d1 ) {
      for( unsigned int d2 = 0; d2 < dim; ++d2 ) {
        errorProjectorITK( d1, d2 ) = errorProjector( d1, d2 );
      }
    }

    // VTK can't have points with type float, so we need to convert everything
    vcl_vector< itk::Vector< float, dim > >  directionsITK;
    for( typename std::vector< itk::Vector< float, dim > >::size_type  s = 0; s < directions.size(); ++s ) {
      itk::Vector< float, dim >  directionITK;
      directionITK[0] = directions[s][0];
      directionITK[1] = directions[s][1];
      directionITK[2] = directions[s][2];
      directionsITK.push_back( directionITK );
    }

    #if COMPUTE_COVARIANCES
    //vnl_matrix_fixed< float, dim, dim >  covariance( 0.0 );
    PointAttribute::SymmetricMatrixType  covariance( 0.0 );
    PointAttribute  attribute( strength, covariance, errorProjectorITK, shape, directionsITK );
    #else
    PointAttribute  attribute( strength, errorProjectorITK, shape, directionsITK );
    #endif
    pointData->push_back( attribute );
    itk::Point< float, dim >  itklocation;
    itklocation[0] = location[0];
    itklocation[1] = location[1];
    itklocation[2] = location[2];
    points->push_back( itklocation );

  }

  std::cout << "Read: " << numCorners << " corners, " << numTubes << " tubes, " << numSheets << " sheets." << std::endl;

}


template<class TOutputMesh>
void
VTKPolyDataToMeshFilter<TOutputMesh>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os,indent);

}


} //end of namespace itk


#endif
