/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkPolyDataToFeaturesWithShapeFilter.h,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkPolyDataToFeaturesWithShapeFilter
// .SECTION Description Convert ITK PointSet to VTK PolyData format for saving

// author: Michal Sofka, Rensselaer Polytechnic Institute (RPI)
// date: 10/25/2007

// .SECTION See Also
// vtkContourFilter vtkCutter vtkEdgePoints vtkExtractEdges
// vtkGeometryFilter vtkGlyph3D vtkHedgeHog vtkHyperStreamline
// vtkMaskPoints vtkOutlineFilter vtkStreamer vtkTensorGlyph
// vtkThresholdPoints vtkVectorTopology

#ifndef __vtkPolyDataToFeaturesWithShapeFilter_h
#define __vtkPolyDataToFeaturesWithShapeFilter_h

#include "vtkPolyDataSource.h"
#include "vtkDataObject.h"
#include "vtkSmartPointer.h"

#include "itkPointSet.h"

//#include "itkFeatureImageFilter.h"

#include <cdcl/cdcl_feature_with_shape.h>

namespace itk {
class PointAttribute;
}

// vtkFeatureWithShapeAttributeSet is a wrapper around stl container (vector) to be used
// in the conversion filter vtkPolyDataToFeaturesWithShapeFilter below
class vtkFeatureWithShapeAttributeSet : public vtkDataObject
{
public:
  static vtkFeatureWithShapeAttributeSet *New();
  vtkTypeRevisionMacro(vtkFeatureWithShapeAttributeSet,vtkDataObject);

  vtkFeatureWithShapeAttributeSet() {}
  ~vtkFeatureWithShapeAttributeSet() {}

  typedef cdcl_feature_with_shape< 3 >  FeatureType;
  typedef FeatureType::sptr  FeaturePointer;
  typedef std::vector< FeaturePointer >  FeatureSetType;

  FeatureSetType &
  GetPoints(void)
  {
    return m_FeatureSet;
  }

  FeatureSetType  m_FeatureSet;

  virtual const char * 	GetClassName () { return "vtkFeatureWithShapeAttributeSet"; };
};


class VTK_FILTERING_EXPORT vtkPolyDataToFeaturesWithShapeFilter : public vtkSource
{
public:
  static vtkPolyDataToFeaturesWithShapeFilter *New();
  vtkTypeRevisionMacro(vtkPolyDataToFeaturesWithShapeFilter,vtkSource);
  void PrintSelf(ostream& os, vtkIndent indent);

  // this class could be eventually templated over attributes
  typedef vtkFeatureWithShapeAttributeSet  AttributeSetType;
  typedef vtkFeatureWithShapeAttributeSet::FeatureSetType  FeatureSetType;
  typedef vtkFeatureWithShapeAttributeSet::FeaturePointer  FeaturePointer;
  typedef vtkFeatureWithShapeAttributeSet::FeatureType  FeatureType;
  // Description:
  // Set / get the input data or filter.
  virtual void SetInput(vtkPolyData*  input);
  vtkPolyData*  GetInput();

  AttributeSetType *GetOutput();
  AttributeSetType *GetOutput(int idx);
  void SetOutput(AttributeSetType *output);
  
  // Description:
  // Do not let images return more than requested.
  virtual void ComputeInputUpdateExtents( vtkDataObject *output );

  void Update();

protected:
  vtkPolyDataToFeaturesWithShapeFilter();
  ~vtkPolyDataToFeaturesWithShapeFilter();

  virtual int FillInputPortInformation(int, vtkInformation*);

  int FillOutputPortInformation(int, vtkInformation*);

  void Execute();
  void GenerateOutput();

  vtkSmartPointer< AttributeSetType >  m_Attributes;

private:
  vtkPolyDataToFeaturesWithShapeFilter(const vtkPolyDataToFeaturesWithShapeFilter&);  // Not implemented.
  void operator=(const vtkPolyDataToFeaturesWithShapeFilter&);  // Not implemented.
};

#endif


