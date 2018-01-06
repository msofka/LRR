/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkITKDescriptorPointSetToPolyDataFilter.h,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkITKDescriptorPointSetToPolyDataFilter
// .SECTION Description Convert ITK PointSet to VTK PolyData format for saving

// author: Michal Sofka, Rensselaer Polytechnic Institute (RPI)
//         sofka at cs dot rpi dot edu
// date: 10/25/2007

// .SECTION See Also
// vtkContourFilter vtkCutter vtkEdgePoints vtkExtractEdges
// vtkGeometryFilter vtkGlyph3D vtkHedgeHog vtkHyperStreamline
// vtkMaskPoints vtkOutlineFilter vtkStreamer vtkTensorGlyph
// vtkThresholdPoints vtkVectorTopology

#ifndef __vtkITKDescriptorPointSetToPolyDataFilter_h
#define __vtkITKDescriptorPointSetToPolyDataFilter_h

#include "vtkPolyDataSource.h"
#include "vtkDataObject.h"

#include "itkPointSet.h"

#include "itkDescriptorMeshFilter.h"

// Point attributes that are stored with each coordinate in the mesh
// (PointSet Pixel type in the PointSet, Mesh and related classes)
//namespace itk {
class DescriptorAttribute;
//}

// In order to be a child of vtkPolyDataSource, input must be of vtkObject type
// so it is not as easy to connect the ITK and VTK filters.
//
// Workaround here is to create a VTK wrapper around an ITK PointSet
// and mimick the pipeline with SetInput and GetOutput methods.
//
// It might be possible to hide this mechanism by including vtkDescriptorAttributeSet
// object in vtkITKDescriptorPointSetToPolyDataFilter
// and providing SetInput( PointSetPointer const &  PointSet ).

class vtkDescriptorAttributeSet : public vtkDataObject
{
public:
  static vtkDescriptorAttributeSet *New();// { return static_cast<vtkDescriptorAttributeSet*>(new itk::PointSet< itk::DescriptorAttribute >); };
  vtkTypeRevisionMacro(vtkDescriptorAttributeSet,vtkDataObject);

  typedef itk::PointSet< itk::DescriptorAttribute >  PointSetType;
  typedef PointSetType::Pointer  PointSetPointer;

  PointSetType::PointDataContainer *
  GetPointData(void) const
  {
    return m_PointSet->GetPointData();
  }

  PointSetType::PointsContainer *
  GetPoints(void) const
  {
    return m_PointSet->GetPoints();
  }

  void
  SetInput( PointSetPointer const &  PointSet )
  {
    m_PointSet = PointSet;
  }

  vtkDescriptorAttributeSet*
  GetOutput(void)
  {
    return this;
  }

  PointSetPointer  m_PointSet;

  virtual const char * 	GetClassName () { return "vtkDescriptorAttributeSet"; };
};


class VTK_FILTERING_EXPORT vtkITKDescriptorPointSetToPolyDataFilter : public vtkPolyDataSource
{
public:
  static vtkITKDescriptorPointSetToPolyDataFilter *New();
  vtkTypeRevisionMacro(vtkITKDescriptorPointSetToPolyDataFilter,vtkPolyDataSource);
  void PrintSelf(ostream& os, vtkIndent indent);

  // this class could be eventually templated over attributes
  typedef vtkDescriptorAttributeSet  FeatureSetType;
  typedef vtkDescriptorAttributeSet::PointSetType  PointSetType;
  typedef vtkDescriptorAttributeSet::PointSetPointer  PointSetPointer;
  //typedef itk::PointSet< itk::DescriptorAttribute >  FeatureSetType;
  // Description:
  // Set / get the input data or filter.
  //virtual void SetInput(FeatureSetType const &  input);
  //const FeatureSetType &  GetInput();
  virtual void SetInput(FeatureSetType*  input);
  FeatureSetType*  GetInput();

  vtkPolyData *GetOutput();
  vtkPolyData *GetOutput(int idx);
  void SetOutput(vtkPolyData *output);
  
  // Description:
  // Do not let images return more than requested.
  virtual void ComputeInputUpdateExtents( vtkDataObject *output );

protected:
  vtkITKDescriptorPointSetToPolyDataFilter();
  ~vtkITKDescriptorPointSetToPolyDataFilter();

  virtual int FillInputPortInformation(int, vtkInformation*);

  int FillOutputPortInformation(int, vtkInformation*);

  void Execute();
  void GenerateOutput();

private:
  vtkITKDescriptorPointSetToPolyDataFilter(const vtkITKDescriptorPointSetToPolyDataFilter&);  // Not implemented.
  void operator=(const vtkITKDescriptorPointSetToPolyDataFilter&);  // Not implemented.
};

#endif


