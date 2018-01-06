/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkMeshSpatialFilterClean.h,v $
  Language:  C++
  Date:      $Date: 2008/04/18 02:09:02 $
  Version:   $Revision: 1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkMeshSpatialFilterClean_h
#define __itkMeshSpatialFilterClean_h

#include "itkMeshToMeshFilter.h"
//#include "itkFeatureImageFilter.h"
#include "itkVectorContainer.h"

namespace itk
{

/** \class MeshSpatialFilterClean
 * \brief 
 *
 * Filter the points using a spatial filter to get good coverage but using only strongest features.
 * During filtering, the following recursive criteria is applied:  a point is retained if no point stronger
 * than it and within min_distance_between_points has also been retained.
 *
 * This filter gets the input mesh, filters the points, and returns sparser set in the output mesh.
 * 
 * \ingroup MeshFilters
 */
template <class TInputMesh, class TOutputMesh>
class ITK_EXPORT MeshSpatialFilterClean : 
    public MeshToMeshFilter<TInputMesh,TOutputMesh>
{
public:
  /** Standard class typedefs. */
  typedef MeshSpatialFilterClean  Self;
  typedef MeshToMeshFilter<TInputMesh,TOutputMesh> Superclass;
  typedef SmartPointer<Self>  Pointer;
  typedef SmartPointer<const Self>  ConstPointer;
  
  typedef TInputMesh InputMeshType;
  typedef TOutputMesh OutputMeshType;
  typedef typename InputMeshType::Pointer InputMeshPointer;
  typedef typename InputMeshType::ConstPointer InputMeshConstPointer;
  typedef typename OutputMeshType::Pointer OutputMeshPointer;
  typedef typename OutputMeshType::PointType        PointType;

  typedef typename TInputMesh::PointsContainer  InputPointsContainer;
  typedef typename TOutputMesh::PointsContainer OutputPointsContainer;

  typedef typename TInputMesh::PointsContainerPointer   InputPointsContainerPointer;
  typedef typename TInputMesh::PointsContainerConstPointer   InputPointsContainerConstPointer;
  typedef typename TOutputMesh::PointsContainerPointer  OutputPointsContainerPointer;

  typedef typename TInputMesh::PointDataContainer  InputPointDataContainer;
  typedef typename TOutputMesh::PointDataContainer OutputPointDataContainer;

  typedef typename TInputMesh::PointDataContainerPointer   InputPointDataContainerPointer;
  typedef typename TInputMesh::PointDataContainerConstPointer   InputPointDataContainerConstPointer;
  typedef typename TOutputMesh::PointDataContainerPointer  OutputPointDataContainerPointer;

  typedef typename itk::Point< float, TInputMesh::PointDimension >  PhysicalPointType;

  /** Type for representing coordinates. */
  typedef typename TInputMesh::CoordRepType  CoordRepType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);
  
  /** Run-time type information (and related methods). */
  itkTypeMacro(MeshSpatialFilterClean,MeshToMeshFilter);

  /** Set minimum distance between points. */
  itkSetMacro(MinDistanceBetweenPoints, float); 

  /** Set minimum distance between points. */
  itkGetMacro(MinDistanceBetweenPoints, float);

protected:
  MeshSpatialFilterClean();
  ~MeshSpatialFilterClean() {};
  void PrintSelf(std::ostream& os, Indent indent) const;
  
  /** Generate Requested Data */
  virtual void GenerateData( void );

  float m_MinDistanceBetweenPoints;

private:
  MeshSpatialFilterClean(const MeshSpatialFilterClean&); //purposely not implemented
  void operator=(const MeshSpatialFilterClean&); //purposely not implemented


  // Function object for sorting vector of keypoint indices
  class StrengthAtIndexIsGreater {
  public:
    StrengthAtIndexIsGreater( InputPointDataContainerConstPointer const &  pointData )
      : m_PointData( pointData ) {};

    bool operator()( int i, int j )
      {
        return (*m_PointData)[i].m_Strength > (*m_PointData)[j].m_Strength;
      }
  private:
    InputPointDataContainerConstPointer const &  m_PointData;
  };

  
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMeshSpatialFilterClean.txx"
#endif

#endif
