/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkVTKPolyDataToMeshFilter.h,v $
  Language:  C++
  Date:      $Date: 2008/07/20 17:20:24 $
  Version:   $Revision: 1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkVTKPolyDataToMeshFilter_h
#define __itkVTKPolyDataToMeshFilter_h

#include "itkMesh.h"
#include "itkMeshSource.h"
#include "itkTriangleCell.h"
#include "itkMapContainer.h"

namespace itk
{


class VTKPolyDataConnector : public Object {
public:

  typedef VTKPolyDataConnector         Self;
  typedef DataObject                   Superclass;
  typedef SmartPointer<Self>           Pointer;

  typedef vtkPolyData  PolyDataType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(VTKPolyDataConnector, DataObject);

  PolyDataType *  m_PolyData;

};

/** \class VTKPolyDataToMeshFilter
 * \brief
 * Reads a vtkPolyData file and create an itkMesh.
 *
 * Caveat: itkVTKPolyDataToMeshFilter can only read triangle meshes.
 *         Use vtkTriangleFilter to convert your mesh to a triangle mesh.
 */
template <class TOutputMesh>
class VTKPolyDataToMeshFilter : public MeshSource<TOutputMesh>
{
public:
  /** Standard "Self" typedef. */
  typedef VTKPolyDataToMeshFilter         Self;
  typedef MeshSource<TOutputMesh>   Superclass;
  typedef SmartPointer<Self>        Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(VTKPolyDataToMeshFilter, MeshSource);

  /** Hold on to the type information specified by the template parameters. */
  typedef TOutputMesh                          OutputMeshType;
  typedef typename OutputMeshType::MeshTraits  MeshTraits;
  typedef typename OutputMeshType::PointType   PointType;
  typedef typename MeshTraits::PixelType       PixelType;

  /** Some convenient typedefs. */
  typedef typename OutputMeshType::Pointer         OutputMeshPointer;
  //typedef typename OutputMeshType::CellTraits      CellTraits;
  //typedef typename OutputMeshType::CellIdentifier  CellIdentifier;
  //typedef typename OutputMeshType::CellType        CellType;
  //typedef typename OutputMeshType::CellAutoPointer CellAutoPointer;
  //typedef typename OutputMeshType::PointIdentifier PointIdentifier;
  //typedef typename CellTraits::PointIdIterator     PointIdIterator;

  typedef typename OutputMeshType::PointsContainerPointer
    PointsContainerPointer;
  
  typedef typename OutputMeshType::PointsContainer
    PointsContainer;

  typedef typename OutputMeshType::PointDataContainer PointDataContainer;
  typedef typename PointDataContainer::Pointer      PointDataContainerPointer;


  /** Define the triangular cell types which form the surface  */
  //typedef TriangleCell<CellType>      TriangleCellType;

  //typedef typename TriangleCellType::SelfAutoPointer
  //  TriangleCellAutoPointer;

  //typedef std::pair<unsigned long,unsigned long>     IndexPairType;
  //typedef MapContainer<IndexPairType, unsigned long> PointMapType;
  //typedef typename PointType::VectorType             VectorType;

  /** Some convenient typedefs. */
  typedef VTKPolyDataConnector PolyDataConnectorType;
  typedef PolyDataConnectorType::PolyDataType  InputPolyDataType;
  //typedef typename InputMeshType::Pointer InputMeshPointer;
  
  /** Set the mesh input of this process object.  */
  void SetInput(InputPolyDataType *  input) { m_Input = input; };

  ///** Set the mesh input of this process object.  */
  //void SetInput(PolyDataConnectorType *  input);

  //void SetInput(InputPolyDataType *  input) {
  //  typename PolyDataConnectorType::Pointer  connector = PolyDataConnectorType::New();
  //  connector->m_PolyData = input;
  //  this->SetInput( connector );
  //}

  ///** Get the mesh input of this process object.  */
  //PolyDataConnectorType * GetInput(void);
  //PolyDataConnectorType * GetInput(unsigned int idx);

protected:
  VTKPolyDataToMeshFilter();
  ~VTKPolyDataToMeshFilter() {}
  void PrintSelf(std::ostream& os, Indent indent) const;

  /** Reads the file */
  void GenerateData();

  InputPolyDataType *  m_Input;

private:
  VTKPolyDataToMeshFilter(const Self&); // purposely not implemented
  void operator=(const Self&); // purposely not implemented
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkVTKPolyDataToMeshFilter.txx"
#endif

#endif //_itkVTKPolyDataToMeshFilter_h
