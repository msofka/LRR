/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkMeshSpatialFilterClean.txx,v $
  Language:  C++
  Date:      $Date: 2008/08/13 21:23:55 $
  Version:   $Revision: 1.4 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkMeshSpatialFilterClean_txx
#define _itkMeshSpatialFilterClean_txx

#include "itkMeshSpatialFilterClean.h"
#include "itkExceptionObject.h"

#include "rsdl_bins.txx"


namespace itk {

  
/**
 *
 */
template <class TInputMesh, class TOutputMesh >
MeshSpatialFilterClean<TInputMesh,TOutputMesh>
::MeshSpatialFilterClean()
{
}


/**
 *
 */
template <class TInputMesh, class TOutputMesh >
void 
MeshSpatialFilterClean<TInputMesh,TOutputMesh>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);
  
  //os << indent << "Num Points: " << m_UnfilteredPoints.size() << std::endl;
  //os << indent << "Num Filtered Points: " << m_FilteredIndices.size() << std::endl;
}


/**
 * This method causes the filter to generate its output.
 */
template <class TInputMesh, class TOutputMesh >
void 
MeshSpatialFilterClean<TInputMesh,TOutputMesh>
::GenerateData(void) 
{
  
  InputMeshConstPointer    inputMesh      =  this->GetInput();
  OutputMeshPointer   outputMesh     =  this->GetOutput();

  InputPointsContainerConstPointer      inputPoints    = inputMesh->GetPoints();
  InputPointDataContainerConstPointer   inputPointData = inputMesh->GetPointData();

  OutputPointsContainerPointer      outputPoints    = OutputPointsContainer::New();
  OutputPointDataContainerPointer   outputPointData = OutputPointDataContainer::New();
  
  if( !inputMesh )
    {
    itkExceptionMacro(<<"Missing Input Mesh");
    }

  if( !outputMesh )
    {
    itkExceptionMacro(<<"Missing Output Mesh");
    }

  outputMesh->SetBufferedRegion( outputMesh->GetRequestedRegion() );



  // Keypoint filtering on the non-maxed points
  //
  unsigned int  max_number_of_points = 1000000;
  //float         m_MinDistanceBetweenPoints = std::min( spacing[0], std::min( spacing[1], spacing[2] ) );
  //float         m_MinDistanceBetweenPoints = std::sqrt( spacing[0]*spacing[0] + spacing[1]*spacing[1] + spacing[2]*spacing[2] );
  //float         m_MinDistanceBetweenPoints = 4.0;


  // vector of indices into the original vector of points (and point data)
  typedef std::vector< typename InputPointDataContainer::VectorContainerSizeType >  IndicesVectorType;
  IndicesVectorType  filteredIndices;
  filteredIndices.reserve( inputPointData->size() );
  for( typename InputPointDataContainer::VectorContainerSizeType  i = 0; i < inputPointData->size(); ++i ) {
    filteredIndices.push_back( i );
  }

  //  Sort the vector of subscripts into the original vector of
  //  features (m_FilteredIndices) according to increasing strength of the
  //  associated feature.
  StrengthAtIndexIsGreater  featureCompare( inputPointData );
  std::sort( filteredIndices.begin(), filteredIndices.end(), featureCompare );


  PointType const &  curr_point = (*inputPoints)[filteredIndices[0]];

  const unsigned int InputImageDimension = TInputMesh::PointDimension;//TInputMesh::TMeshTraits::PointDimension;
  vnl_vector_fixed<float,InputImageDimension>  min_point( curr_point.GetVnlVector() );
  vnl_vector_fixed<float,InputImageDimension>  max_point( curr_point.GetVnlVector() );
  for ( typename IndicesVectorType::size_type i=1; i< filteredIndices.size(); ++i )
  {
    PointType const &  loc = (*inputPoints)[filteredIndices[i]];

    for ( int k=0; k<InputImageDimension; ++k )
    {
      if ( loc[k] < min_point[k] ) 
        min_point[k] = loc[k];
      else if ( loc[k] > max_point[k] )
        max_point[k] = loc[k];
    }
  }

  //  Assign the bin sizes to be isotropic.  The default size is 2
  //  times the min point distance.  The subsequent for loop ensures
  //  that a reasonable size is set when the min_distance is too small.

  float bsize = m_MinDistanceBetweenPoints;
  const int max_bins_per_dim = 200;
  for ( unsigned int i=0; i<InputImageDimension; ++i )
  {
    float lower_bound_on_size = (max_point[i] - min_point[i]) / max_bins_per_dim;
    if ( bsize < lower_bound_on_size ) bsize = lower_bound_on_size;
  }
  vnl_vector_fixed<float,InputImageDimension>  bin_sizes( bsize );
  //  vcl_cout << "feature_filter:  min_point = " << min_point
  //   << ", max_point = " << max_point << ", bsize = " << bsize
  //   << ", m_MinDistanceBetweenPoints = " << m_MinDistanceBetweenPoints << vcl_endl;

  // The n-d version could be improved dramatically!  In particular,
  // the is_any_point_within_radius function.

  rsdl_bins< InputImageDimension, float, int >  subscript_bins( min_point, max_point, bin_sizes );





  // build a table of exponential function to save computation
  // Matlab plot: mindist=5.0; x = 0:0.1:20; for i=1:length(x), if x(i) > 3.0*mindist, y(i) = 0.0; else y(i) = 1.0 - 0.4*exp( -x(i)*x(i) / (mindist*mindist) ); end, end, plot( x, y );
  std::vector< double >  exponential;
  double resolution = 0.1;
  double sigma = m_MinDistanceBetweenPoints;
  for( double e = 0.0; e < 3.0*sigma; e += resolution ) {
    double scaled = e / sigma;
    // the larger the constant in front of the exponential, the more features we will get
    exponential.push_back( 1.0 - 0.5*std::exp( -scaled*scaled ) );
  }

  //  For each remaining point (in order of strength) only put it into
  //  the bins and the filtered result if there isn't already a point
  //  within the min distance of its physical location.

  typename vcl_vector< PointType >::size_type  subscript = 0;
  while ( subscript < filteredIndices.size() &&
    outputPointData->size() < max_number_of_points )
  {
    PointAttribute const &  currPointData = (*inputPointData)[ filteredIndices[subscript] ];
    PointType const &  currPoint = (*inputPoints)[ filteredIndices[subscript] ];
    PhysicalPointType const &  pointPhysical = currPoint;

    //if ( ! subscript_bins.is_any_point_within_radius( feature->location_, m_MinDistanceBetweenPoints ) )
    double strength_fraction = 1.0;

    vcl_vector< vnl_vector_fixed< float, InputImageDimension > >  nearest_neighbors;
    vcl_vector<int>  indices;
    subscript_bins.n_nearest( pointPhysical.GetVnlVector(), 1, nearest_neighbors, indices );


    // if the region in which we are searching for the nearest neighbor is not square/cube, then
    // once boundary along any dimension is reached during the search, the search is terminated
    // this can happen when there are not that many points in the bins yet, so include the point
    if( nearest_neighbors.size() == 0 ) {
      subscript_bins.add_point( pointPhysical.GetVnlVector(), filteredIndices[subscript] );
      outputPoints->push_back( (*inputPoints)[filteredIndices[subscript]] );
      outputPointData->push_back( (*inputPointData)[filteredIndices[subscript]] );
      ++subscript;
      continue;
    }


    double distance = ( pointPhysical.GetVnlVector() - nearest_neighbors[0] ).magnitude();

    // if the distance between current point and nearest neighbor already added is greater than 3*sigma (sigma = m_MinDistanceBetweenPoints)
    //   add this point
    // if the distance is smaller but greater than 0.5*sigma, decide based on the strength, points already added will always have higher strength
    //   so lower this by an exponential value and add the current point if it has strength higher than this lowered strength of an added point
    //
    // if there are too many points close to each other, increase the 1.0 minimum tolerated distance
    // if there are remaining too many points in the background, increase m_MinDistanceBetweenPoints
    // if there are points filtered along low contrast boundaries, decrease m_MinDistanceBetweenPoints
    //if( distance > 0.0 * sigma ) { // no spatial filtering (turn off the filter for debugging other components)
    if( distance > 2.0 * sigma ) {
      subscript_bins.add_point( pointPhysical.GetVnlVector(), filteredIndices[subscript] );
      outputPoints->push_back( (*inputPoints)[filteredIndices[subscript]] );
      outputPointData->push_back( (*inputPointData)[filteredIndices[subscript]] );
    }
    // if there is point closer than min_distance between points, compute fraction by which the strength
    // must be larger than already existing point
    else if( distance > 0.5 * sigma ) {
      unsigned int exponential_index = (unsigned int) (distance / resolution);
      strength_fraction = exponential[exponential_index];  // must be at least k time larger (k is a constant)

      float const &  nearestNeighborStrength = (*inputPointData)[indices[0]].m_Strength;
      float const &  currentPointStrength = currPointData.m_Strength;

      //std::cout << "cur: " << currentPointStrength << " nn: " << nearestNeighborStrength << " fr: " << strength_fraction << std::endl;
      if( strength_fraction * nearestNeighborStrength < currentPointStrength )
      {
        //vcl_cout << original_features[ indices[0] ]->strength_ << " " << original_features[ sorted_subscripts[subscript] ]->strength_ << " " << strength_fraction << vcl_endl;
        subscript_bins.add_point( pointPhysical.GetVnlVector(), filteredIndices[subscript] );
        outputPoints->push_back( (*inputPoints)[filteredIndices[subscript]] );
        outputPointData->push_back( (*inputPointData)[filteredIndices[subscript]] );
      }

    }

    //printf( "Filtered fraction %.6f features\r", double( subscript ) / double( m_FilteredIndices.size() ) );

    ++subscript;
  }
  std::cout << std::endl;  // printf only returned carriage


  outputMesh->SetPoints( outputPoints );

  outputMesh->SetPointData(  outputPointData );


  std::cout << "Points before filtering: " << inputPoints->size() << " filtered: " << outputPoints->size() << std::endl;


}

} // end namespace itk

#endif
