/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkKdTreeForThreading.h,v $
  Language:  C++
  Date:      $Date: 2009-03-04 15:23:51 $
  Version:   $Revision: 1.28 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkKdTreeForThreading_h
#define __itkKdTreeForThreading_h

#include <queue>
#include <vector>

#include "itkMacro.h"
#include "itkPoint.h"
#include "itkSize.h"
#include "itkObject.h"
#include "itkNumericTraits.h"
#include "itkArray.h"

#include "itkSample.h"
#include "itkSubsample.h"

#include "itkEuclideanDistance.h"
#include "itkKdTree.h"  // for nodes

namespace itk {
namespace Statistics {


// Michal Sofka: This code is modified so that search can be run from a multithreaded pipeline.


template < class TSample >
class ITK_EXPORT KdTreeForThreading : public Object
{
public:
  /** Standard class typedefs */
  typedef KdTreeForThreading                   Self;
  typedef Object                   Superclass;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Run-time type information (and related methods) */
  itkTypeMacro(KdTreeForThreading, Object);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** typedef alias for the source data container */
  typedef TSample                                 SampleType;
  typedef typename TSample::MeasurementVectorType MeasurementVectorType;
  typedef typename TSample::MeasurementType       MeasurementType;
  typedef typename TSample::InstanceIdentifier    InstanceIdentifier;
  typedef typename TSample::FrequencyType         FrequencyType;

  typedef unsigned int                    MeasurementVectorSizeType;

  /** Get Macro to get the length of a measurement vector in the KdTreeForThreading.
   * The length is obtained from the input sample. */
  itkGetConstMacro( MeasurementVectorSize, MeasurementVectorSizeType );

  /** DistanceMetric type for the distance calculation and comparison */
  typedef EuclideanDistance< MeasurementVectorType > DistanceMetricType;

  /** Node type of the KdTreeForThreading */
  typedef KdTreeNode< TSample > KdTreeNodeType;

  /** Neighbor type. The first element of the std::pair is the instance
   * identifier and the second one is the distance between the measurement
   * vector identified by the first element and the query point. */
  typedef std::pair< InstanceIdentifier, double > NeighborType;

  typedef std::vector< InstanceIdentifier > InstanceIdentifierVectorType;

  /** Sets the number of measurement vectors that can be stored in a
   * terminal node */
  void SetBucketSize(unsigned int size);

  /** Sets the input sample that provides the measurement vectors to the k-d
   * tree */
  void SetSample(const TSample* sample);

  /** Returns the pointer to the input sample */
  const TSample* GetSample() const
    { return m_Sample; }

  unsigned long Size() const
    { return m_Sample->Size(); }

  /** Returns the pointer to the empty terminal node. A KdTreeForThreading object
   * has a single empty terminal node in memory. when the split process
   * has to create an empty terminal node, the single instance is reused
   * for this case */
  KdTreeNodeType* GetEmptyTerminalNode()
    { return m_EmptyTerminalNode; }

  /** Sets the root node of the KdTree that is a result of
   * KdTreeGenerator or WeightedCentroidKdTreeGenerator. */
  void SetRoot(KdTreeNodeType* root)
    { m_Root = root; }

  /** Returns the pointer to the root node. */
  KdTreeNodeType* GetRoot()
    { return m_Root; }

  /** Returns the measurement vector identified by the instance
   * identifier that is an identifier defiend for the input sample */
  const MeasurementVectorType & GetMeasurementVector(InstanceIdentifier id) const
    { return m_Sample->GetMeasurementVector(id); }

  /** Returns the frequency of the measurement vector identified by
   * the instance identifier */
  FrequencyType GetFrequency(InstanceIdentifier id) const
    { return m_Sample->GetFrequency( id ); }

  /** Get the pointer to the distance metric. */
  DistanceMetricType* GetDistanceMetric()
    { return m_DistanceMetric.GetPointer(); }

  /** Searches the neighbors fallen into a hypersphere */
  void Search(const MeasurementVectorType &query,
              double radius,
              InstanceIdentifierVectorType& result) const;

  /** Returns true if the intermediate k-nearest neighbors exist within
   * the the bounding box defined by the lowerBound and the
   * upperBound. Otherwise returns false. Returns false if the ball
   * defined by the distance between the query point and the farthest
   * neighbor touch the surface of the bounding box. */
  bool BallWithinBounds(const MeasurementVectorType &query,
                        MeasurementVectorType &lowerBound,
                        MeasurementVectorType &upperBound,
                        double radius) const;

  /** Returns true if the ball defined by the distance between the query
   * point and the farthest neighbor overlaps with the bounding box
   * defined by the lower and the upper bounds. */
  bool BoundsOverlapBall(const MeasurementVectorType &query,
                         MeasurementVectorType &lowerBound,
                         MeasurementVectorType &upperBound,
                         double radius) const;

  /** Deletes the node recursively */
  void DeleteNode(KdTreeNodeType *node);

  /** Prints out the tree information */
  void PrintTree( std::ostream & os ) const;

  /** Prints out the tree information */
  void PrintTree(KdTreeNodeType *node, unsigned int level,
                 unsigned int activeDimension,
                 std::ostream & os = std::cout ) const;

  /** Draw out the tree information to a ostream using 
   * the format of the Graphviz dot tool. */
  void PlotTree( std::ostream & os ) const;

  /** Prints out the tree information */
  void PlotTree(KdTreeNodeType *node, std::ostream & os = std::cout ) const;


  typedef typename TSample::Iterator      Iterator;
  typedef typename TSample::ConstIterator ConstIterator;

  Iterator Begin()
    {
    typename TSample::ConstIterator iter = m_Sample->Begin();
    return iter;
    }

  Iterator End()
    {
    Iterator iter = m_Sample->End();
    return iter;
    }

  ConstIterator Begin() const
    {
    typename TSample::ConstIterator iter = m_Sample->Begin();
    return iter;
    }

  ConstIterator End() const
    {
    ConstIterator iter = m_Sample->End();
    return iter;
    }

protected:
  /** Constructor */
  KdTreeForThreading();

  /** Destructor: deletes the root node and the empty terminal node. */
  virtual ~KdTreeForThreading();

  void PrintSelf(std::ostream& os, Indent indent) const;

  /** search loop */
  int SearchLoop(const KdTreeNodeType* node, const MeasurementVectorType &query,
                 MeasurementVectorType &lowerBound,
                 MeasurementVectorType &upperBound,
                 double const &  searchRadius,
                 InstanceIdentifierVectorType &  neighbors ) const;
private:
  KdTreeForThreading(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  /** Pointer to the input sample */
  const TSample* m_Sample;

  /** Number of measurement vectors can be stored in a terminal node. */
  int m_BucketSize;

  /** Pointer to the root node */
  KdTreeNodeType* m_Root;

  /** Pointer to the empty terminal node */
  KdTreeNodeType* m_EmptyTerminalNode;

  /** Distance metric smart pointer */
  typename DistanceMetricType::Pointer m_DistanceMetric;

  /** Measurement vector size */
  MeasurementVectorSizeType m_MeasurementVectorSize;
}; // end of class

} // end of namespace Statistics
} // end of namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkKdTreeForThreading.txx"
#endif

#endif
