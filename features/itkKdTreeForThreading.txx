/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkKdTreeForThreading.txx,v $
  Language:  C++
  Date:      $Date: 2009-03-04 15:23:51 $
  Version:   $Revision: 1.31 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkKdTreeForThreading_txx
#define __itkKdTreeForThreading_txx

#include "itkKdTreeForThreading.h"

namespace itk { 
namespace Statistics {


template< class TSample >
KdTreeForThreading< TSample >
::KdTreeForThreading()
{
  m_EmptyTerminalNode = 
    new KdTreeTerminalNode< TSample >();

  m_DistanceMetric = DistanceMetricType::New();
  m_Sample = 0;
  m_Root = 0;
  m_BucketSize = 16;
  m_MeasurementVectorSize = 0;
}

template< class TSample >
KdTreeForThreading< TSample >
::~KdTreeForThreading()
{
  if ( m_Root != 0 )
    {
    this->DeleteNode(m_Root);
    }
  
  delete m_EmptyTerminalNode;
}

template< class TSample >
void
KdTreeForThreading< TSample >
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);
  
  os << indent << "Input Sample: ";
  if ( m_Sample != 0 )
    {
    os << m_Sample << std::endl;
    }
  else
    {
    os << "not set." << std::endl;
    }

  os << indent << "Bucket Size: " << m_BucketSize << std::endl;
  os << indent << "Root Node: ";
  if ( m_Root != 0 )
    {
    os << m_Root << std::endl;
    }
  else
    {
    os << "not set." << std::endl;
    }
  os << indent << "MeasurementVectorSize: " << 
            m_MeasurementVectorSize << std::endl;
}

template< class TSample >
void
KdTreeForThreading< TSample >
::DeleteNode(KdTreeNodeType *node)
{
  if ( node->IsTerminal() )
    {
    // terminal node
    if (node == m_EmptyTerminalNode)
      {
      // empty node
      return;
      }
    delete node;
    return;
    }

  // non-terminal node
  if ( node->Left() != 0 )
    {
    this->DeleteNode( node->Left() );
    }

  if ( node->Right() != 0 )
    {
    this->DeleteNode( node->Right() );
    }

  delete node;
}

template< class TSample >
void
KdTreeForThreading< TSample >
::SetSample(const TSample* sample)
{
  m_Sample = sample;
  this->m_MeasurementVectorSize = m_Sample->GetMeasurementVectorSize();
  this->m_DistanceMetric->SetMeasurementVectorSize( 
                      this->m_MeasurementVectorSize );
  this->Modified();
}

template< class TSample >
void
KdTreeForThreading< TSample >
::SetBucketSize(unsigned int size)
{
  m_BucketSize = size;
}

template< class TSample >
void
KdTreeForThreading< TSample >
::Search(const MeasurementVectorType &query, double radius,
         InstanceIdentifierVectorType& result) const
{
  MeasurementVectorType lowerBound;
  MeasurementVectorType upperBound;
  MeasurementVectorTraits::SetLength( lowerBound, m_MeasurementVectorSize );
  MeasurementVectorTraits::SetLength( upperBound, m_MeasurementVectorSize );
  
  for (unsigned int d = 0; d < this->m_MeasurementVectorSize; d++)
    {
    lowerBound[d] = static_cast< MeasurementType >( -vcl_sqrt( -static_cast<double>( NumericTraits< MeasurementType >::NonpositiveMin() ) ) / 2.0 );
    upperBound[d] = static_cast< MeasurementType >(  vcl_sqrt(  static_cast<double>( NumericTraits< MeasurementType >::max() ) / 2.0 ) );
    }

  this->SearchLoop(m_Root, query, lowerBound, upperBound, radius, result);
}

template< class TSample >
inline int
KdTreeForThreading< TSample >
::SearchLoop(const KdTreeNodeType* node, 
             const MeasurementVectorType &query, 
             MeasurementVectorType &lowerBound,
             MeasurementVectorType &upperBound,
             double const &  searchRadius,
             InstanceIdentifierVectorType &  neighbors ) const
{
  unsigned int i;
  InstanceIdentifier tempId;
  double tempDistance;

  if ( node->IsTerminal() )
    {
    // terminal node
    if (node == m_EmptyTerminalNode)
      {
      // empty node
      return 0;
      }

    for (i = 0; i < node->Size(); i++)
      {
      tempId = node->GetInstanceIdentifier(i);
      tempDistance = 
        m_DistanceMetric->
        Evaluate(query, m_Sample->GetMeasurementVector(tempId));
      if (tempDistance < searchRadius )
        {
          neighbors.push_back(tempId);
        }
      }

    if ( this->BallWithinBounds(query, lowerBound, upperBound, 
                          searchRadius) )
      {
      return 1;
      }

    return 0;
    }


  unsigned int partitionDimension; 
  MeasurementType partitionValue;
  MeasurementType tempValue;
  node->GetParameters(partitionDimension, partitionValue);

  if (query[partitionDimension] <= partitionValue)
    {
    // search the closer child node
    tempValue = upperBound[partitionDimension];
    upperBound[partitionDimension] = partitionValue;
    if (SearchLoop(node->Left(), query, lowerBound, upperBound, searchRadius, neighbors))
      {
      return 1;
      }
    upperBound[partitionDimension] = tempValue;

    // search the other node, if necessary
    tempValue = lowerBound[partitionDimension];
    lowerBound[partitionDimension] = partitionValue;
    if ( this->BoundsOverlapBall(query, lowerBound, upperBound, 
                           searchRadius) )
      {
      SearchLoop(node->Right(), query, lowerBound, upperBound, searchRadius, neighbors);
      }
    lowerBound[partitionDimension] = tempValue;
    }
  else
    {
    // search the closer child node
    tempValue = lowerBound[partitionDimension];
    lowerBound[partitionDimension] = partitionValue;
    if (SearchLoop(node->Right(), query, lowerBound, upperBound, searchRadius, neighbors))
      {
      return 1;
      }
    lowerBound[partitionDimension] = tempValue;

    // search the other node, if necessary
    tempValue = upperBound[partitionDimension];
    upperBound[partitionDimension] = partitionValue;
    if ( this->BoundsOverlapBall(query, lowerBound, upperBound, 
                           searchRadius) )
      {
      SearchLoop(node->Left(), query, lowerBound, upperBound, searchRadius, neighbors);
      }
    upperBound[partitionDimension] = tempValue;
    }

  // stop or continue search
  if ( this->BallWithinBounds(query, lowerBound, upperBound, 
                        searchRadius) )
    {
    return 1;
    }  

  return 0;
}

template< class TSample >
inline bool
KdTreeForThreading< TSample >
::BallWithinBounds(const MeasurementVectorType &query, 
                   MeasurementVectorType &lowerBound,
                   MeasurementVectorType &upperBound,
                   double radius) const
{
  unsigned int dimension;
  for (dimension = 0; dimension < this->m_MeasurementVectorSize; dimension++)
    {
    if ((m_DistanceMetric->Evaluate(query[dimension] ,
                                    lowerBound[dimension]) <= 
         radius) ||
        (m_DistanceMetric->Evaluate(query[dimension] ,
                                    upperBound[dimension]) <= 
         radius))
      {
      return false;
      }
    }
  return true;
}

template< class TSample >
inline bool
KdTreeForThreading< TSample >
::BoundsOverlapBall(const MeasurementVectorType &query, 
                    MeasurementVectorType &lowerBound,
                    MeasurementVectorType &upperBound,
                    double radius) const
{
  double sum = NumericTraits< double >::Zero;
  double temp;
  unsigned int dimension;
  double squaredSearchRadius = radius * radius;
  for (dimension = 0; dimension < m_MeasurementVectorSize; dimension++)
    {

    if (query[dimension] <= lowerBound[dimension])
      {
      temp = m_DistanceMetric->Evaluate(query[dimension], 
                                        lowerBound[dimension]);
      sum += temp * temp;
      if (sum < squaredSearchRadius)
        {
        return true;
        }
      }
    else if (query[dimension] >= upperBound[dimension])
      {
      temp = m_DistanceMetric->Evaluate(query[dimension], 
                                        upperBound[dimension]);
      sum += temp * temp;
      if (sum < squaredSearchRadius)
        {
        return true;
        }
      }
    }
  return false;
}


template< class TSample >
void
KdTreeForThreading< TSample >
::PrintTree( std::ostream & os ) const
{
  const unsigned int topLevel = 0;
  const unsigned int activeDimension = 0;
  this->PrintTree( this->m_Root, topLevel, activeDimension, os );
}


template< class TSample >
void
KdTreeForThreading< TSample >
::PrintTree(KdTreeNodeType *node, unsigned int level, unsigned int activeDimension, std::ostream & os ) const
{
  level++;
  if ( node->IsTerminal() )
    {
    // terminal node
    if (node == m_EmptyTerminalNode)
      {
      // empty node
      os << "Empty node: level = " << level << std::endl;
      return;
      }

    os << "Terminal: level = " << level
       << " dim = " << activeDimension<< std::endl;
    os << "          ";
    for (unsigned int i = 0; i < node->Size(); i++)
      {
      os << "[" << node->GetInstanceIdentifier(i) << "] "
                << m_Sample->GetMeasurementVector(node->GetInstanceIdentifier(i)) << ", ";
      }
    os << std::endl;
    return;
    }
  
  unsigned int partitionDimension;
  MeasurementType partitionValue;

  node->GetParameters(partitionDimension, partitionValue);
  typename KdTreeNodeType::CentroidType centroid;
  node->GetWeightedCentroid(centroid);
  os << "Nonterminal: level = " << level << std::endl;
  os << "             dim = " << partitionDimension << std::endl;
  os << "             value = " << partitionValue << std::endl;
  os << "             weighted centroid = " << centroid;
  os << "             size = " << node->Size()<< std::endl;
  os << "             identifier = " << node->GetInstanceIdentifier(0);
  os << m_Sample->GetMeasurementVector(node->GetInstanceIdentifier(0)) << std::endl;
 
  this->PrintTree( node->Left(),  level, partitionDimension, os );
  this->PrintTree( node->Right(), level, partitionDimension, os );
}


template< class TSample >
void
KdTreeForThreading< TSample >
::PlotTree( std::ostream & os ) const
{
  // 
  // Graph header
  //
  os << "digraph G {" << std::endl;

  //
  // Recursively visit the tree and add entries for the nodes
  //
  this->PlotTree( this->m_Root, os );

  // 
  // Graph footer
  //
  os << "}" << std::endl;
}


template< class TSample >
void
KdTreeForThreading< TSample >
::PlotTree(KdTreeNodeType *node, std::ostream & os ) const
{
  unsigned int partitionDimension;
  MeasurementType partitionValue;

  node->GetParameters(partitionDimension, partitionValue);

  KdTreeNodeType * left  = node->Left();
  KdTreeNodeType * right = node->Right();

  char partitionDimensionCharSymbol = ('X'+partitionDimension);

  if( node->IsTerminal() )
    {
    // terminal node
    if( node != m_EmptyTerminalNode )
      {
      os << "\"" << node << "\" [label=\"";
      for( unsigned int i = 0; i < node->Size(); i++ )
        {
        os << this->GetMeasurementVector( node->GetInstanceIdentifier(i) );
        os << " ";
        }
      os << "\" ];" << std::endl;
      }
    }
  else
    {
    os << "\"" << node << "\" [label=\"";
    os << this->GetMeasurementVector( node->GetInstanceIdentifier(0) );
    os << " " << partitionDimensionCharSymbol << "=" << partitionValue;
    os << "\" ];" << std::endl;
    }


  if( left &&  ( left != m_EmptyTerminalNode ) )
    { 
    os << "\"" << node << "\" -> \"" << left << "\";" << std::endl;
    this->PlotTree( left, os );
    }

  if( right && ( right != m_EmptyTerminalNode ) )
    {
    os << "\"" << node << "\" -> \"" << right << "\";" << std::endl;
    this->PlotTree( right, os );
    }
}

} // end of namespace Statistics 
} // end of namespace itk

#endif
