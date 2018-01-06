/*=========================================================================


author: Michal Sofka, Rensselaer Polytechnic Institute
  date: 10/10/2007
  

=========================================================================*/
#ifndef _itkFeatureImageFilter_txx
#define _itkFeatureImageFilter_txx
#include "itkFeatureImageFilter.h"

#include "itkInterpolateImageFunction.h"
#include "itkVectorLinearInterpolateImageFunction.h"
#include "itkConstNeighborhoodIterator.h"

#include "itkImageRegionIterator.h"
#include "itkNeighborhoodAlgorithm.h"
#include "itkZeroFluxNeumannBoundaryCondition.h"
#include "itkConstantBoundaryCondition.h"

#include "itkGaussianDerivativeImageFunction.h"
#include "itkDerivativeImageFilter.h"
#include "itkGaussianDerivativeImageFunction.h"
#include "itkRecursiveGaussianImageFilter.h"
#include "itkMultiResolutionPyramidImageFilter.h"

#include "itkKdTree.h"
#include "itkKdTreeGenerator.h"
#include "itkListSample.h"



#include <vnl/algo/vnl_determinant.h>
#include <vnl/vnl_trace.h>
#include <vnl/vnl_cross.h>

// for debug
#include "itkNumericSeriesFileNames.h"
#include "itkImageSeriesWriter.h"
#include "itkRescaleIntensityImageFilter.h"

#include "itkTimeProbe.h"

#define DEBUG 0

#if DEBUG
#define PRINTF(...)  printf(__VA_ARGS__)
#else
#define PRINTF(...)  
#endif


namespace itk
{

template <class TInputImage, class TOutputMesh>
FeatureImageFilter< TInputImage, TOutputMesh>
::FeatureImageFilter()
{

  // Modify superclass default values, can be overridden by subclasses
  this->SetNumberOfRequiredInputs(1);

  PointDataContainerPointer  pointData  = PointDataContainer::New();
  OutputMeshPointer mesh = this->GetOutput();
  mesh->SetPointData( pointData.GetPointer() );

  m_Sigma = 1.00;

  m_OnlyCorners = false;

  m_GradientImage = 0;
  m_ScoreImage = 0;
  m_FilteredPointsImage = 0;
  m_NonMaxedPoints = 0;
}


template <class TInputImage, class TOutputMesh>
void
FeatureImageFilter< TInputImage, TOutputMesh>
::GenerateOutputInformation()
{
}


template <class TInputImage, class TOutputMesh>
void
FeatureImageFilter< TInputImage, TOutputMesh>
::SetInput( const InputImageType * inputImage )
{

  // This const_cast is needed due to the lack of
  // const-correctness in the ProcessObject.
  this->SetNthInput( 0,
            const_cast< InputImageType * >( inputImage ) );

}


namespace {

template <class T, class Ran>
void rrel_util_median_and_scale( Ran first, Ran last,
                                 T& median, T& scale,
                                 int dof )
{
  long count = long(last-first); // VC6 & 7 has broken iterator_traits
  assert( count > 0 );
  assert( count > dof );

  //// ignore these values in median computation
  //T  notComputed = -1.0;
  //long notComputedCount = 0;

  //for ( Ran i=first; i!=last; ++i ) {
  //  if( *i == notComputed ) ++notComputedCount;
  //}

  Ran loc = first + count/2;// + notComputedCount;
  std::nth_element( first, loc, last );
  median = *loc;
  for ( Ran i=first; i!=last; ++i ) {
    *i = vnl_math_abs(*i-median);
  }
  ++loc;
  std::nth_element( first, loc, last );
  scale = T(1.4826 * (1 + 5.0/(count-dof)) * *loc);
}

template <class T, class InpIter>
void rrel_util_median_and_scale_copy( InpIter first, InpIter last,
                                      T& median, T& scale,
                                      int dof )
{
  // FIXME: scratch should be std::vector<
  // std::iterator_traits<InpIter>::value_type >, but this is not
  // supported under all compilers. In particular, VC++ doesn't
  // support it for vector iterators.
  //
  std::vector<T> scratch;
  for ( ; first != last; ++first )
    scratch.push_back( *first );
  rrel_util_median_and_scale( scratch.begin(), scratch.end(), median, scale, dof );
}


namespace {

// Function object for sorting vector of keypoint indices
template <class InternalRecord>
class FeatureAtIndexIsGreater {
public:
  FeatureAtIndexIsGreater( std::vector< InternalRecord > const &  points )
    : m_Points( points ) {};

  bool operator()( int i, int j )
    {
      return m_Points[i].m_Strength > m_Points[j].m_Strength;
    }
private:
  std::vector< InternalRecord > const &  m_Points;
};

}

}


//----------------------------------------------------------------------------
template <class TInputImage, class TOutputMesh>
int 
FeatureImageFilter< TInputImage, TOutputMesh>
::SplitRequestedRegion( int i, int num, typename InternalImageType::RegionType &  splitRegion, bool shrink )
{
  // Get the output pointer
  InternalImageType * outputPtr = this->m_ScoreImage;
  typename InternalImageType::RegionType const &  requestedRegionConst
    = outputPtr->GetRequestedRegion();

  typename InternalImageType::SizeType  requestedRegionConstSize
    = requestedRegionConst.GetSize();

  typename InternalImageType::IndexType  requestedRegionConstIndex
    = requestedRegionConst.GetIndex();

  typename InternalImageType::RegionType  requestedRegion;
  if( shrink ) {
    // Shrink region by 2 in each dimension.
    requestedRegionConstSize[0] -= 4;
    requestedRegionConstSize[1] -= 4;
    requestedRegionConstSize[2] -= 4;
    requestedRegionConstIndex[0] += 2;
    requestedRegionConstIndex[1] += 2;
    requestedRegionConstIndex[2] += 2;
    
    requestedRegion.SetSize( requestedRegionConstSize );
    requestedRegion.SetIndex( requestedRegionConstIndex );
  }
  else {
    requestedRegion = requestedRegionConst;
  }


  const typename InternalImageType::SizeType &  requestedRegionSize
    = requestedRegion.GetSize();

  int splitAxis;
  typename InternalImageType::IndexType splitIndex;
  typename InternalImageType::SizeType splitSize;

  // Initialize the splitRegion to the output requested region
  splitRegion = requestedRegion;
  splitIndex = splitRegion.GetIndex();
  splitSize = splitRegion.GetSize();

  // split on the outermost dimension available
  splitAxis = outputPtr->GetImageDimension() - 1;
  while (requestedRegionSize[splitAxis] == 1)
    {
    --splitAxis;
    if (splitAxis < 0)
      { // cannot split
      itkDebugMacro("  Cannot Split");
      return 1;
      }
    }

  // determine the actual number of pieces that will be generated
  typename InternalImageType::SizeType::SizeValueType range = requestedRegionSize[splitAxis];
  int valuesPerThread = (int)::ceil(range/(double)num);
  int maxThreadIdUsed = (int)::ceil(range/(double)valuesPerThread) - 1;

  // Split the region
  if (i < maxThreadIdUsed)
    {
    splitIndex[splitAxis] += i*valuesPerThread;
    splitSize[splitAxis] = valuesPerThread;
    }
  if (i == maxThreadIdUsed)
    {
    splitIndex[splitAxis] += i*valuesPerThread;
    // last thread needs to process the "rest" dimension being split
    splitSize[splitAxis] = splitSize[splitAxis] - i*valuesPerThread;
    }
  
  // set the split region ivars
  splitRegion.SetIndex( splitIndex );
  splitRegion.SetSize( splitSize );

  itkDebugMacro("  Split Piece: " << splitRegion );
  //std::cout << "  Split Piece: " << splitRegion << std::endl;

  return maxThreadIdUsed + 1;
}


template <class TInputImage, class TOutputMesh>
int 
FeatureImageFilter< TInputImage, TOutputMesh>
::SplitRequestedList( int i, int num,
                     typename std::list< InternalRecord >::iterator &  startIterator,
                     typename std::list< InternalRecord >::iterator &  endIterator )
{
  // determine the actual number of pieces that will be generated
  typename InternalImageType::SizeType::SizeValueType range = m_UnfilteredPoints.size();
  int valuesPerThread = (int)::ceil(range/(double)num);
  int maxThreadIdUsed = (int)::ceil(range/(double)valuesPerThread) - 1;

  unsigned int startIndex;
  unsigned int endIndex;
  // Split the region
  if (i < maxThreadIdUsed)
    {
    startIndex = i * valuesPerThread;
    endIndex = (i + 1) * valuesPerThread;
    }
  if (i == maxThreadIdUsed)
    {
    startIndex = i * valuesPerThread;
    // last thread needs to process the "rest" of the vector
    //endIndex = range - 1;
    endIndex = range;

    }
  startIterator = m_UnfilteredPoints.begin();
  endIterator = m_UnfilteredPoints.begin();
  for( unsigned int ind = 0; ind < startIndex; ++ind ) ++startIterator;
  for( unsigned int ind = 0; ind < endIndex; ++ind ) ++endIterator;

  itkDebugMacro("  Split Piece: " << startIndex << " " << endIndex );
  //std::cout << "  Split Piece: " << startIndex << " " << endIndex << std::endl;

  return maxThreadIdUsed + 1;
}


// The execute method created by the subclass.
template <class TInputImage, class TOutputMesh>
void 
FeatureImageFilter< TInputImage, TOutputMesh>
::ComputeScoreImageThreaded( typename InternalImageType::RegionType const &  outputRegion,
                             int  threadId )
{
  typedef typename GradientImageFilterType::OutputImageType  GradientImageType;

  typename InputImageType::SizeType  imageSize  = outputRegion.GetSize();
  typename InputImageType::IndexType  imageIndex = outputRegion.GetIndex();
  InputSizeType  radiusOne = { 1, 1, 1 };

  //ConstantBoundaryCondition<InputImageType> cbc;
  //cbc.SetConstant( NumericTraits<InputPixelType>::NonpositiveMin() );

  // Find the data-set boundary "faces"
  typedef typename NeighborhoodAlgorithm::ImageBoundaryFacesCalculator< GradientImageType >  FaceListCalculatorType;
  FaceListCalculatorType  bC;
  typename FaceListCalculatorType::FaceListType  faceList = bC( m_GradientImage, outputRegion, radiusOne );

  typename FaceListCalculatorType::FaceListType::iterator  fit;

  // Compute the score image
  //
  // Process each of the boundary faces.  These are N-d regions which border
  // the edge of the buffer.
  typename FaceListCalculatorType::FaceListType::iterator  fitEnd = faceList.begin();
  for (fit=faceList.begin(); fit != faceList.end(); ++fit)
    {
    ConstNeighborhoodIterator< GradientImageType >  bit = ConstNeighborhoodIterator< GradientImageType >( radiusOne, m_GradientImage, *fit );
    unsigned int neighborhoodSize = bit.Size();
    double oneOverNeighborhoodSize = 1.0 / neighborhoodSize;
    //bit.OverrideBoundaryCondition(&cbc);

    // compute score image
    for( bit.GoToBegin(); ! bit.IsAtEnd(); ++bit )
      {
      typename GradientImageType::IndexType  curr_index = bit.GetIndex();

      // Compute outer product of the gradients and score value
      //
      typename GradientImageType::PixelType  gradient2;
      gradient2.Fill( 0.0 );
      for (unsigned int i = 0; i < neighborhoodSize; ++i)
        {
        typename GradientImageType::PixelType  gradient = bit.GetPixel( i );

        // Diagonal of the outer product.
        gradient2[0] += gradient[0] * gradient[0];
        gradient2[1] += gradient[1] * gradient[1];
        gradient2[2] += gradient[2] * gradient[2];
        }

      // Trace of the outer product.
      double curr_score = oneOverNeighborhoodSize * ( gradient2[0] + gradient2[1] + gradient2[2] );

      m_ScoreImage->SetPixel( curr_index, curr_score );

      #if DEBUG
      double curr_one_index = (curr_index[0]-imageIndex[0]) + (curr_index[1]-imageIndex[1])*imageSize[0] + (curr_index[2]-imageIndex[2])*imageSize[0]*imageSize[1];
      double curr_one_size  = imageSize[0] * imageSize[1] * imageSize[2];
      PRINTF( "Computed %.4f of scores\r", curr_one_index / curr_one_size );
      //std::cout << curr_index << imageSize << curr_one_index << "  " << curr_one_size << "  " << curr_one_index / curr_one_size << std::endl;
      #endif
      }
    }

}


// Callback routine used by the threading library. This routine just calls
// the ThreadedGenerateData method after setting the correct region for this
// thread. 
template <class TInputImage, class TOutputMesh>
ITK_THREAD_RETURN_TYPE 
FeatureImageFilter< TInputImage, TOutputMesh>
::ThreaderCallbackScore( void *arg )
{
  ThreadStruct *str;
  int total, threadId, threadCount;

  threadId = ((MultiThreader::ThreadInfoStruct *)(arg))->ThreadID;
  threadCount = ((MultiThreader::ThreadInfoStruct *)(arg))->NumberOfThreads;

  str = (ThreadStruct *)(((MultiThreader::ThreadInfoStruct *)(arg))->UserData);

  // execute the actual method with appropriate output region
  // first find out how many pieces extent can be split into.
  typename TInputImage::RegionType splitRegion;
  total = str->Filter->SplitRequestedRegion(threadId, threadCount,
                                            splitRegion);

  if (threadId < total)
    {
    str->Filter->ComputeScoreImageThreaded( splitRegion, threadId );
    }
  // else
  //   {
  //   otherwise don't use this thread. Sometimes the threads dont
  //   break up very well and it is just as efficient to leave a 
  //   few threads idle.
  //   }
  
  return ITK_THREAD_RETURN_VALUE;
}


template <class TInputImage, class TOutputMesh>
ITK_THREAD_RETURN_TYPE 
FeatureImageFilter< TInputImage, TOutputMesh>
::ThreaderCallbackSubvoxel( void *arg )
{
  ThreadStruct *str;
  int total, threadId, threadCount;

  threadId = ((MultiThreader::ThreadInfoStruct *)(arg))->ThreadID;
  threadCount = ((MultiThreader::ThreadInfoStruct *)(arg))->NumberOfThreads;

  str = (ThreadStruct *)(((MultiThreader::ThreadInfoStruct *)(arg))->UserData);

  // execute the actual method with appropriate output region
  // first find out how many pieces extent can be split into.
  std::list< InternalRecord >::iterator  startIterator;
  std::list< InternalRecord >::iterator  endIterator;
  total = str->Filter->SplitRequestedList(threadId, threadCount,
                                          startIterator, endIterator);

  if (threadId < total)
    {
    str->Filter->SubvoxelAccuracyThreaded( startIterator, endIterator, threadId );
    }
  // else
  //   {
  //   otherwise don't use this thread. Sometimes the threads dont
  //   break up very well and it is just as efficient to leave a 
  //   few threads idle.
  //   }
  
  return ITK_THREAD_RETURN_VALUE;
}


// Callback routine used by the threading library. This routine just calls
// the ThreadedGenerateData method after setting the correct region for this
// thread. 
template <class TInputImage, class TOutputMesh>
ITK_THREAD_RETURN_TYPE 
FeatureImageFilter< TInputImage, TOutputMesh>
::ThreaderCallbackNonMax( void *arg )
{
  ThreadStruct *str;
  int total, threadId, threadCount;

  threadId = ((MultiThreader::ThreadInfoStruct *)(arg))->ThreadID;
  threadCount = ((MultiThreader::ThreadInfoStruct *)(arg))->NumberOfThreads;

  str = (ThreadStruct *)(((MultiThreader::ThreadInfoStruct *)(arg))->UserData);

  // execute the actual method with appropriate output region
  // first find out how many pieces extent can be split into.
  typename TInputImage::RegionType splitRegion;
  total = str->Filter->SplitRequestedRegion(threadId, threadCount,
                                            splitRegion, true );

  if (threadId < total)
    {
    str->Filter->NonMaximumSuppressionThreaded( splitRegion, threadId );
    }
  // else
  //   {
  //   otherwise don't use this thread. Sometimes the threads dont
  //   break up very well and it is just as efficient to leave a 
  //   few threads idle.
  //   }
  
  return ITK_THREAD_RETURN_VALUE;
}


template <class TInputImage, class TOutputMesh>
void
FeatureImageFilter< TInputImage, TOutputMesh>
::LocalContrastFiltering( typename InternalImageType::Pointer &  scoreImage )
{
  ConstantBoundaryCondition<InputImageType> cbc;
  cbc.SetConstant( NumericTraits<InputPixelType>::NonpositiveMin() );

  typename FilteredPointsIndicatorImageType::PixelType const  scoreTooLow = 255;

  typename InputImageType::SizeType  imageSize  = scoreImage->GetRequestedRegion().GetSize();
  typename InputImageType::IndexType  imageIndex = scoreImage->GetRequestedRegion().GetIndex();

  for( unsigned int d = 0; d < 3; ++d ) {
    // if any of the image dimensions too small, don't do filtering
    if( imageSize[d] / 2 < m_FilteringRadius[d] ) {
      std::cout << "Warning: image very small, skipping contrast filtering." << std::endl;
      return;
    }
  }


  // Local contrast filtering
  //
  typedef typename NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<InternalImageType>  FaceListCalculatorType;
  FaceListCalculatorType  facesCalculatorFiltering;
  typename FaceListCalculatorType::FaceListType  faceListFiltering = facesCalculatorFiltering(scoreImage, scoreImage->GetRequestedRegion(), m_FilteringRadius);
  typename FaceListCalculatorType::FaceListType::iterator  facesCalculatorFilteringEnd;

  unsigned int numFilteredOut = 0;
  // HACK: Since we are directly using the neighborhood and accessing the image, we can only use the central face (not boundary faces)
  // because the boundary condition would not get used.
  // Best would be to write this in terms of a new neighborhood iterator which steps in unit step of physical coordinates.
  // UPDATE: This hack has been avoided below, but more permanent solution is needed
  typename FaceListCalculatorType::FaceListType::iterator facesFilteringIt;
  typename FaceListCalculatorType::FaceListType::iterator facesFilteringEnd;
  facesFilteringEnd = faceListFiltering.begin();
  ++facesFilteringEnd;
  // if this happens (currently implemented), then boundary pixels do not have scores but still get used in neighborhood
  // filtering below (and for mean computation); therefore mean uses -1.0 (=notComputed) value
  std::cout << "Warning: nbrhood for contrast filt. same as for score img. computation." << std::endl;
  for (facesFilteringIt=faceListFiltering.begin(); facesFilteringIt != facesFilteringEnd; ++facesFilteringIt)
  //for (facesFilteringIt=faceListFiltering.begin(); facesFilteringIt != faceListFiltering.end(); ++facesFilteringIt)
    {
    //ConstNeighborhoodIterator< InternalImageType > bit = ConstNeighborhoodIterator<InternalImageType>(m_FilteringRadius,
    //                                                                                                  m_ScoreImage, m_ScoreImage->GetRequestedRegion());

    //NeighborhoodIterator< FilteredPointsIndicatorImageType > iit = NeighborhoodIterator<FilteredPointsIndicatorImageType>(m_FilteringRadius,
    //                                                                                                                      m_FilteredPointsImage, m_ScoreImage->GetRequestedRegion() );

    ConstNeighborhoodIterator< InternalImageType > bit = ConstNeighborhoodIterator<InternalImageType>(m_FilteringRadius,
                                                                                                      scoreImage, *facesFilteringIt);
    NeighborhoodIterator< FilteredPointsIndicatorImageType > iit = NeighborhoodIterator<FilteredPointsIndicatorImageType>(m_FilteringRadius,
                                                                                                                          m_FilteredPointsImage, *facesFilteringIt);
    typename ConstNeighborhoodIterator<InputImageType>::OffsetType  halfOffset;
    halfOffset[0] = m_FilteringRadius[0];
    halfOffset[1] = m_FilteringRadius[1];
    halfOffset[2] = m_FilteringRadius[2];
    unsigned int neighborhoodSize = bit.Size();
    //bit.OverrideBoundaryCondition(&cbc);
    bit.GoToBegin();
    iit.GoToBegin();

    typename InternalImageType::IndexType  start_index = bit.GetIndex();    

    // Filter score image
    bool within_region = false;
    do
      {
      typename ConstNeighborhoodIterator<InternalImageType>::NeighborhoodType const &  neighborhood = bit.GetNeighborhood();

      // FIXING THE HACK ABOVE: remove "notComputed" values from score neighborhood container so that they are not included in median computation
      std::vector< typename InternalImageType::PixelType >  neighborhoodWithoutNotComputed;
      const float scoreNotComputed = -1.0f;
      for( unsigned long  g = 0; g < neighborhood.Size(); ++g ) {
        if( neighborhood[g] != scoreNotComputed ) neighborhoodWithoutNotComputed.push_back( neighborhood[g] );
      }

      float median, median_scale;
      rrel_util_median_and_scale_copy( neighborhoodWithoutNotComputed.begin(), neighborhoodWithoutNotComputed.end(), median, median_scale, 1 );
      //rrel_util_median_and_scale_copy( neighborhood.Begin(), neighborhood.End(), median, median_scale, 1 );

      // Filter features in this neighborhood
      //
      const float ratioOfStdDev = -0.5f;
      const float threshold = median + ratioOfStdDev * median_scale;
      //std::cout << threshold << ": ";
      for (unsigned int i = 0; i < neighborhoodSize; ++i)
        {
        //std::cout << bit.GetPixel(i) << " ";
        if( bit.GetPixel(i) < threshold && iit.GetPixel(i) != scoreTooLow )
          {
          iit.SetPixel( i, scoreTooLow );
          //InputIndexType  index = bit.GetIndex(i);
          //m_FilteredPointsImage->SetPixel( index, scoreTooLow );
          ++numFilteredOut;
          }
        }
      //std::cout << std::endl;

      typename InternalImageType::IndexType  curr_index = bit.GetIndex();
      #if DEBUG
      double curr_one_index = (curr_index[0]-imageIndex[0]) + (curr_index[1]-imageIndex[1])*imageSize[0] + (curr_index[2]-imageIndex[2])*imageSize[0]*imageSize[1];
      double curr_one_size  = imageSize[0] * imageSize[1] * imageSize[2];
      PRINTF( "Filtered %.4f of scores\r", curr_one_index / curr_one_size );
      #endif
      //std::cout << curr_index << imageSize << curr_one_index << "  " << curr_one_size << "  " << curr_one_index / curr_one_size << std::endl;

      // step size is radius
      //bit += halfOffset;
      //if( !bit.InBounds() ) break;

      // manually stepping by half offset
      within_region = false;
      for( unsigned int d = 0; d < 3; ++d ) {
        curr_index[d] += halfOffset[d];
        bit.SetLocation( curr_index );
        iit.SetLocation( curr_index );
        if( bit.InBounds() ) {
          within_region = true;
          break;
        }
        else {
          curr_index[d] = start_index[d];
          bit.SetLocation( curr_index );
          iit.SetLocation( curr_index );
        }
      }

      } while( within_region );
    }
  std::cout << "Scores filtered out in contrast filtering: " << numFilteredOut << " in image size: " << imageSize[0]*imageSize[1]*imageSize[2] << std::endl;
  std::cout << std::endl;  // printf only returned carriage

}


template <class TInputImage, class TOutputMesh>
void
FeatureImageFilter< TInputImage, TOutputMesh>
::NonMaximumSuppressionThreaded( typename InternalImageType::RegionType const &  outputRegion,
                                 int  threadId )
{
  typename InputImageType::SizeType  imageSize  = outputRegion.GetSize();
  typename InputImageType::IndexType  imageIndex = outputRegion.GetIndex();

  typename FilteredPointsIndicatorImageType::PixelType const  scoreTooLow = 255;

  typedef typename GradientImageFilterType::OutputImageType  GradientImageType;

  // Interpolator for linearly intepolating the score 
  typedef itk::LinearInterpolateImageFunction< InternalImageType > InterpolatorType;
  typename InterpolatorType::Pointer  interp = InterpolatorType::New();
  interp->SetInputImage( m_ScoreImage );

  // Non-max suppression
  //
  unsigned int num_wout_maxima = 0;

  InputSizeType  radiusOne = { 1, 1, 1 };
  ConstNeighborhoodIterator< GradientImageType >  git = ConstNeighborhoodIterator< GradientImageType >( radiusOne, m_GradientImage, outputRegion );
  unsigned int neighborhoodSize = git.Size();
  double oneOverNeighborhoodSize = 1.0 / neighborhoodSize;

  ImageRegionConstIteratorWithIndex< InternalImageType > bit = ImageRegionConstIteratorWithIndex< InternalImageType >( m_ScoreImage, outputRegion );
  for( bit.GoToBegin(), git.GoToBegin(); ! bit.IsAtEnd(); ++bit, ++git )
  //for( bit.GoToBegin(); ! bit.IsAtEnd(); ++bit )
    {

    typename InternalImageType::IndexType  curr_index = bit.GetIndex();

    double curr_score = m_ScoreImage->GetPixel( curr_index );

    // Decide about a feature type by searching for a maximum along eigenvector directions
    //
    // Very small threshold on the score value (if the outer product matrix is just noise,
    // the code below fails because we can't get meaningful directions).
    if ( curr_score < 0.0001 ) {
      continue;
    };

    if( m_FilteredPointsImage->GetPixel( curr_index ) == scoreTooLow ) {
      continue;
    }

    InternalRecord  resultRecord;


    // Compute outer product of the gradients.
    //
    vnl_matrix_fixed< float, 3, 3 >  outerProductVnl( 0.0 );
    for (unsigned int i = 0; i < neighborhoodSize; ++i)
      {
      //typename GradientImageType::PixelType  gradient = gradientImage->GetPixel( index );
      typename GradientImageType::PixelType  gradient = git.GetPixel( i );
      vnl_vector_fixed< float, InputImageDimension >  gradientVector( gradient.GetDataPointer() );

      outerProductVnl += outer_product( gradientVector, gradientVector );
      }
    outerProductVnl *= oneOverNeighborhoodSize;
    SymmetricMatrixType  outerProduct;
    outerProduct( 0, 0 ) = outerProductVnl( 0, 0 );
    outerProduct( 0, 1 ) = outerProductVnl( 0, 1 );
    outerProduct( 0, 2 ) = outerProductVnl( 0, 2 );
    outerProduct( 1, 1 ) = outerProductVnl( 1, 1 );
    outerProduct( 1, 2 ) = outerProductVnl( 1, 2 );
    outerProduct( 2, 2 ) = outerProductVnl( 2, 2 );


    // ordering is ascending by default
    SymmetricMatrixType::EigenValuesArrayType  eigenValues;
    SymmetricMatrixType::EigenVectorsMatrixType  eigenVectors;
    outerProduct.ComputeEigenAnalysis( eigenValues, eigenVectors );


    // Find number of maxima along each dimension
    unsigned int num_maxima = 0;
    for( unsigned int d = 0; d < InputImageDimension; ++d ) // for each eigenvector
      {
      // Step in each direction and find out if the points have larger strength.
      // Do this with continuous index rather then physical point since the interpolation using physical points converts
      // the points to physical indices.
      //
      itk::Vector< float, InputImageDimension >   step_direction;
      step_direction[0] = eigenVectors( 0, d );
      step_direction[1] = eigenVectors( 1, d );
      step_direction[2] = eigenVectors( 2, d );

      // Test previous step.
      //
      InterpolatorType::ContinuousIndexType  prev_index = curr_index;
      prev_index[0] -= step_direction[0];
      prev_index[1] -= step_direction[1];
      prev_index[2] -= step_direction[2];

      // get neighbors using linear interpolation
      double prev_score = interp->EvaluateAtContinuousIndex( prev_index );

      if( curr_score > prev_score ) {

        // Test next step.
        //
        InterpolatorType::ContinuousIndexType  post_index = curr_index;
        post_index[0] += step_direction[0];
        post_index[1] += step_direction[1];
        post_index[2] += step_direction[2];

        double post_score = interp->EvaluateAtContinuousIndex( post_index );

        //std::cout << prev_score << " " << curr_score << " " << post_score << std::endl;

        // compare only immediate neighbors
        if( curr_score > post_score )
          {
          ++num_maxima;
          resultRecord.m_Directions.push_back( step_direction );
          }
        }

      } // end for each eigenvector


    //std::cout << "num_maxima: " << num_maxima << std::endl;
    // determine a feature type (for non-max suppression and sub-pixel accuracy) based on number of maxima along each dimension
    switch( num_maxima )
      {
      case 0: ++num_wout_maxima;
              continue; // no extremum -> no feature point
              break;
      case 1: {
              resultRecord.m_PointType = InternalRecord::Sheet;
              break;
              }
      case 2: {
              resultRecord.m_PointType = InternalRecord::Tube;
              break;
              }
      case 3: {
              resultRecord.m_PointType = InternalRecord::Corner;
              //resultRecord.m_Strength = curr_score;
              }
      }

    assert( num_maxima > 0 );

    m_ScoreImage->TransformIndexToPhysicalPoint( curr_index, resultRecord.m_Location );

    resultRecord.m_LocationIndex = curr_index;
    #if COMPUTE_COVARIANCES
    resultRecord.m_OuterProduct = outerProduct;
    #endif
    // push_back slows down multithreading
    m_NonMaxedPoints[threadId].push_back( resultRecord );

    #if DEBUG
    // progress output
    double curr_one_index = (curr_index[0]-imageIndex[0]) + (curr_index[1]-imageIndex[1])*imageSize[0] + (curr_index[2]-imageIndex[2])*imageSize[0]*imageSize[1];
    double curr_one_size  = imageSize[0] * imageSize[1] * imageSize[2];
    PRINTF( "Detected %.2d features, scores w/out max %.2d, fraction done %.2f\r", nonMaxedPoints.size(), num_wout_maxima, curr_one_index / curr_one_size );
    #endif


    } // end while bit!=IsAtEnd()
  //std::cout << std::endl;  // printf only returned carriage

  //std::cout << "Number of points after non-max suppression: " << m_NonMaxedPoints[threadId].size() << std::endl;
}


template <class TInputImage, class TOutputMesh>
void
FeatureImageFilter< TInputImage, TOutputMesh>
::ComputeDirections( std::list< InternalRecord > &  unfilteredPoints )
{
  typedef typename GradientImageFilterType::OutputImageType  GradientImageType;
  typedef itk::VectorLinearInterpolateImageFunction< GradientImageType > GradientInterpolatorType;
  typename GradientInterpolatorType::Pointer  interpGrad = GradientInterpolatorType::New();
  interpGrad->SetInputImage( m_GradientImage );

  typename GradientImageType::SpacingType  spacing = m_GradientImage->GetSpacing();
  double minSpacing = std::min( spacing[0], std::min( spacing[1], spacing[2] ) );
  typename GradientImageType::SizeType  radiusOne = { 1, 1, 1 };

  itk::ConstNeighborhoodIterator< GradientImageType >  bit = itk::ConstNeighborhoodIterator<GradientImageType>( radiusOne, m_GradientImage, m_GradientImage->GetRequestedRegion() );
  unsigned int neighborhoodSize = bit.Size();
  double oneOverNeighborhoodSize = 1.0 / neighborhoodSize;

  bit.GoToBegin();
  // stepping by -1, 0, 1 in physical coordinates in the local neighborhood
  // this is the same as stepping by minSpacing / spacing in the local neighborhood in terms of continous index
  // this precomputation saves time and using continuous index saves conversion in the interpolator
  std::vector< typename GradientInterpolatorType::ContinuousIndexType >  offsets;
  // note: it does not matter that we are using index at the beginning as the location is relative
  typename GradientImageType::IndexType  curr_index = bit.GetIndex();
  for (unsigned int i = 0; i < neighborhoodSize; ++i) {
    typename GradientImageType::IndexType  indexLocal = bit.GetIndex(i);
    typename GradientInterpolatorType::ContinuousIndexType  continuousIndex;
    continuousIndex[0] = ( indexLocal[0] - curr_index[0] ) * minSpacing / spacing[0];
    continuousIndex[1] = ( indexLocal[1] - curr_index[1] ) * minSpacing / spacing[1];
    continuousIndex[2] = ( indexLocal[2] - curr_index[2] ) * minSpacing / spacing[2];
    offsets.push_back( continuousIndex );
  }


  for( std::list< InternalRecord >::iterator  it = unfilteredPoints.begin(); it != unfilteredPoints.end(); ++it ) {

		itk::ContinuousIndex< float, 3 >  curr_cont_index;
    m_GradientImage->TransformPhysicalPointToContinuousIndex( it->m_Location, curr_cont_index );

    // Improve the directions based on averaging over neighborhood centered at the point physical location.
    //
    vnl_vector_fixed< double, 3 >  directionLocal( 0.0 );
    for (unsigned int n = 0; n < neighborhoodSize; ++n)
      {
        // stepping by -1, 0, 1 in physical coordinates in the local neighborhood
        typename GradientInterpolatorType::ContinuousIndexType  currIndexNeigh;
        currIndexNeigh[0] = curr_cont_index[0] + offsets[n][0];
        currIndexNeigh[1] = curr_cont_index[1] + offsets[n][1];
        currIndexNeigh[2] = curr_cont_index[2] + offsets[n][2];

        if( interpGrad->IsInsideBuffer( currIndexNeigh ) ) { // we could step outside the image
          typename GradientInterpolatorType::OutputType  gradient = interpGrad->EvaluateAtContinuousIndex( currIndexNeigh );
          vnl_vector_fixed< double, 3 >  gradientVector( gradient.GetDataPointer() );
          //gradientVector = gradientVector.normalize();
          directionLocal += gradientVector;
        }


      }

    directionLocal *= oneOverNeighborhoodSize;
    directionLocal.normalize();

    if( it->m_PointType == InternalRecord::Sheet ) { // only one direction
      it->m_Directions[0][0] = directionLocal[0];
      it->m_Directions[0][1] = directionLocal[1];
      it->m_Directions[0][2] = directionLocal[2];
    }
    else if( it->m_PointType == InternalRecord::Tube ) {
      // In non-max suppression, we saved the eigenvectors in order.
      // Retrieve it so that we can compute the other directions below.
      vnl_vector_fixed< double, 3 >  evector; // corresponding to the smallest outer product matrix eigenvalue (along the tube)
      evector[0] = it->m_Directions[0][0];
      evector[1] = it->m_Directions[0][1];
      evector[2] = it->m_Directions[0][2];

      // store the first direction
      it->m_Directions[0][0] = directionLocal[0];
      it->m_Directions[0][1] = directionLocal[1];
      it->m_Directions[0][2] = directionLocal[2];

      // limiting to half space for keypoint binormal (ellipsoid is mirrored and from the keypoint point of view
      // the two halves have the same orientation) -> flip the binormal
      if( evector[1] < 0.0 ) evector *= -1.0;

      // gram-schmidt orthogonalization; evector (binormal) projected onto keypoint_direction and this is then subtracted from evector
      // --> resulting in the vector orthogonal to keypoint_direction
      vnl_vector_fixed< double, 3 >  orthogonal = evector - dot_product( directionLocal, evector ) * directionLocal;

      orthogonal.normalize();
      it->m_Directions[1][0] = orthogonal[0];
      it->m_Directions[1][1] = orthogonal[1];
      it->m_Directions[1][2] = orthogonal[2];
      
    }
    else { // corner, it->m_PointType == InternalRecord::Corner

      // In non-max suppression, we saved the eigenvectors in order.
      // Retrieve it so that we can compute the other directions below.
      vnl_vector_fixed< double, 3 >  evector; // corresponding to the largest outer product matrix eigenvalue
      evector[0] = it->m_Directions[2][0];
      evector[1] = it->m_Directions[2][1];
      evector[2] = it->m_Directions[2][2];

      // store the first direction
      it->m_Directions[0][0] = directionLocal[0];
      it->m_Directions[0][1] = directionLocal[1];
      it->m_Directions[0][2] = directionLocal[2];

      // limiting to half space for keypoint binormal (ellipsoid is mirrored and from the keypoint point of view
      // the two halves have the same orientation) -> flip the binormal
      if( evector[1] < 0.0 ) evector *= -1.0;

      // gram-schmidt orthogonalization; evector (binormal) projected onto keypoint_direction and this is then subtracted from evector
      // --> resulting in the vector orthogonal to keypoint_direction
      vnl_vector_fixed< double, 3 >  orthogonal = evector - dot_product( directionLocal, evector ) * directionLocal;

      orthogonal.normalize();
      it->m_Directions[1][0] = orthogonal[0];
      it->m_Directions[1][1] = orthogonal[1];
      it->m_Directions[1][2] = orthogonal[2];

      vnl_vector_fixed< double, 3 >  tangent = vnl_cross_3d( directionLocal, orthogonal );
      tangent.normalize();
      it->m_Directions[2][0] = tangent[0];
      it->m_Directions[2][1] = tangent[1];
      it->m_Directions[2][2] = tangent[2];
    }

  }
}

template <class TInputImage, class TOutputMesh>
void
FeatureImageFilter< TInputImage, TOutputMesh>
::SubvoxelAccuracyThreaded( typename std::list< InternalRecord >::iterator &  startIterator,
                            typename std::list< InternalRecord >::iterator &  endIterator,
                            int threadId )
{
  // Subvoxel accuracy
  //

  // Precompute matrix A for least squares, where X = ( A^T A )^-1 A^T Y
  // What is needed here is to step in the cube neighborhood -1, 0, 1
  //
  ConstNeighborhoodIterator< InternalImageType >  sit;
  InputSizeType  radiusOne = { 1, 1, 1 };
  sit = ConstNeighborhoodIterator< InternalImageType >( radiusOne,
                                                        m_ScoreImage, m_ScoreImage->GetRequestedRegion() );

  sit.GoToBegin();
  InputIndexType  pointIndex = sit.GetIndex();

  typedef typename PhysicalPointType::CoordRepType  CoordRepType;
  itk::Matrix< float, 10, 10 >  A;
  A.Fill( 0.0f );
  for (unsigned int i = 0; i < sit.Size(); ++i)
    {
      InputIndexType  XYZindex = sit.GetIndex(i);

      CoordRepType Xi = static_cast<CoordRepType>( XYZindex[0] ) -  static_cast<CoordRepType>( pointIndex[0] );
      CoordRepType Yi = static_cast<CoordRepType>( XYZindex[1] ) -  static_cast<CoordRepType>( pointIndex[1] );
      CoordRepType Zi = static_cast<CoordRepType>( XYZindex[2] ) -  static_cast<CoordRepType>( pointIndex[2] );

      itk::Vector< float, 10 >  Y;
      Y[0] = 1.0;
      Y[1] = Zi;
      Y[2] = Yi;
      Y[3] = Xi;
      Y[4] = Zi * Zi;
      Y[5] = Yi * Zi;
      Y[6] = Xi * Zi;
      Y[7] = Yi * Yi;
      Y[8] = Xi * Yi;
      Y[9] = Xi * Xi;

      for( unsigned int r = 0; r < 10; ++r )
        for( unsigned int c = 0; c < 10; ++c ) {
          A[r][c] += Y[9-c] * Y[r];
        }

    }

  // this casting needs to be here since the constructor of itkMatrix is explicit
	// without the casting, the constructor is not called
  typedef itk::Matrix< float, 10, 10 >  MatrixType10;
  MatrixType10  At = MatrixType10( A.GetTranspose() );
  MatrixType10  AtA = At * A;
  MatrixType10  AtAinv = MatrixType10( AtA.GetInverse() );
  MatrixType10  LS = AtAinv * At;

//  vnl_matrix< float > At = A.GetTranspose();
//  vnl_matrix< float > AtA = At * A.GetVnlMatrix();
//  vnl_matrix< float > LSvnl = vnl_matrix_inverse< float >( AtA )*At;
//  itk::Matrix< float, 10, 10 > const &  LS = LSvnl;

  //std::cout << "LS: " << std::endl << LS << std::endl;



  typename InterpolatorType::Pointer  interpolator = InterpolatorType::New();
  interpolator->SetInputImage( m_ScoreImage );

  // Subvoxel localization. Remove points that are not at maximum.
  unsigned int counter = 0;
  for( std::list< InternalRecord >::iterator  it = startIterator; it != endIterator; ++it )
    {
    bool localMax = SubvoxelLocalizationLS( *it, interpolator, LS );
    if( localMax ) m_NonMaxedPoints[threadId].push_back( *it );
    }


  //std::cout << "Points that are a local maximum after fit: " << m_NonMaxedPoints[threadId].size() << std::endl;

}


template <class TInputImage, class TOutputMesh>
void
FeatureImageFilter< TInputImage, TOutputMesh>
::ComputeCovariances( std::list< InternalRecord > const &  unfilteredPoints )
{
  OutputMeshPointer           mesh      = this->GetOutput();
  
  PointsContainerPointer      points    = mesh->GetPoints();
  PointDataContainerPointer   pointData = mesh->GetPointData();

  // gamma-normalization to normalize strengths accross scales
  double sheetNormalization = m_Sigma;
  double tubeNormalization = m_Sigma * m_Sigma;
  double cornerNormalization = std::pow( m_Sigma, 1.5 );

  unsigned int numTubes = 0;
  unsigned int numCorners = 0;
  unsigned int numSheets = 0;
  for( std::list< InternalRecord >::const_iterator  it = unfilteredPoints.begin(); it != unfilteredPoints.end(); ++it )
    {

    #if COMPUTE_COVARIANCES
    // Compute outer product of the gradients and score value
    //
    vnl_matrix_fixed< double, 3, 3 >  outerProductVnl( 0.0 );
    for( unsigned int k = 0; k < 3; ++k )
      for( unsigned int l = 0; l < 3; ++l )
        outerProductVnl( k, l ) = it->m_OuterProduct( k, l );

    vnl_svd< double > svd( outerProductVnl );

    //vnl_vector_fixed< double, 3 >  location = (*points)[i].GetVnlVector();
    //double refinedStrength = 0.0;
    //FeaturePointerType  feature = new FeatureType( location, refinedStrength, svd.pinverse() );
    vnl_matrix_fixed< double, 3, 3 >  inverse = svd.pinverse();
    PointAttribute::SymmetricMatrixType  covariance;
    covariance( 0, 0 ) = inverse( 0, 0 );
    covariance( 0, 1 ) = inverse( 0, 1 );
    covariance( 0, 2 ) = inverse( 0, 2 );
    covariance( 1, 1 ) = inverse( 1, 1 );
    covariance( 1, 2 ) = inverse( 1, 2 );
    covariance( 2, 2 ) = inverse( 2, 2 );
    #endif

    vnl_matrix_fixed< float, 3, 3 > errorProjectorVnl;
    PointAttribute::FeatureShape  shape = PointAttribute::Corner;
    double strength = 0.0;
    if( it->m_PointType == InternalRecord::Sheet ) {
      errorProjectorVnl = outer_product( it->m_Directions[0].GetVnlVector(), it->m_Directions[0].GetVnlVector() );
      shape = PointAttribute::Sheet;
      strength = sheetNormalization * it->m_Strength;
      ++numSheets;
    }
    else if( it->m_PointType == InternalRecord::Tube ) {
      vnl_vector_fixed< float, 3 >  tangent = vnl_cross_3d( it->m_Directions[0].GetVnlVector(), it->m_Directions[1].GetVnlVector() );
      errorProjectorVnl.set_identity();
      errorProjectorVnl -= outer_product( tangent, tangent );
      shape = PointAttribute::Tube;
      strength = tubeNormalization * it->m_Strength;
      ++numTubes;
    }
    else if( it->m_PointType == InternalRecord::Corner ) {
      errorProjectorVnl.set_identity();
      shape = PointAttribute::Corner;
      strength = cornerNormalization * it->m_Strength;
      ++numCorners;
    }
    else std::cout << "Error: Unknown feature type." << std::endl;

    PointAttribute::SymmetricMatrixType  errorProjector;
    errorProjector( 0, 0 ) = errorProjectorVnl( 0, 0 );
    errorProjector( 0, 1 ) = errorProjectorVnl( 0, 1 );
    errorProjector( 0, 2 ) = errorProjectorVnl( 0, 2 );
    errorProjector( 1, 1 ) = errorProjectorVnl( 1, 1 );
    errorProjector( 1, 2 ) = errorProjectorVnl( 1, 2 );
    errorProjector( 2, 2 ) = errorProjectorVnl( 2, 2 );
 
    #if COMPUTE_COVARIANCES
    PointAttribute  attribute( strength, covariance, errorProjector, shape, it->m_Directions );
    #else
    PointAttribute  attribute( strength, errorProjector, shape, it->m_Directions );
    #endif 
    pointData->push_back( attribute );
    points->push_back( it->m_Location );

    PRINTF( "Computed %.4f of covariances\r", double(n) / double(filteredIndices.size()) );

    }
  std::cout << std::endl;  // printf only returned carriage
  std::cout << "Detected " << numSheets << " Sheets, " << numTubes << " Tubes, " << numCorners << " Corners." << std::endl;

}


template <class TInputImage, class TOutputMesh>
void
FeatureImageFilter< TInputImage, TOutputMesh>
::GenerateData()
{
  InputImageConstPointer      input     = this->GetInput(0);

  PointsContainerPointer      points    = PointsContainer::New();
  PointDataContainerPointer   pointData = PointDataContainer::New();

  OutputMeshPointer           mesh      = this->GetOutput();
  mesh->SetPoints( points );
  mesh->SetPointData( pointData );




  
  //typedef itk::MultiResolutionPyramidImageFilter< TInputImage, TInputImage >  MultiResolutionFilterType;
  typedef itk::MultiResolutionPyramidImageFilter< TInputImage, InternalImageType >  MultiResolutionFilterType;
  typename MultiResolutionFilterType::Pointer  multiResolutionFilter = MultiResolutionFilterType::New();

  unsigned int numberOfLevels = 3;
  multiResolutionFilter->SetInput( input );
  multiResolutionFilter->SetNumberOfLevels( numberOfLevels );
  multiResolutionFilter->Update();

  const unsigned int numberOfSigmas = 1;
  //double sigmas[numberOfSigmas] = {1.0, 1.26, 1.59, 2.0};
  //const unsigned int numberOfSigmas = 2;
  //double sigmas[numberOfSigmas] = {1.0, 2.0};
  double sigmas[numberOfSigmas] = {1.0};
  //const unsigned int numberOfSigmas = 3;
  //double sigmas[numberOfSigmas] = {1.0, 1.26, 1.59};
  //const unsigned int numberOfSigmas = 1;
  //double sigmas[numberOfSigmas] = {1.0};

  for( unsigned int level = 1; level < numberOfLevels; ++level ) {
  //for( unsigned int level = 1; level < 2; ++level ) {
    std::cout << "Processing level " << level << "..." << std::endl;

    typename MultiResolutionFilterType::OutputImagePointer  inputImage = multiResolutionFilter->GetOutput( level );


    for( unsigned int sigmaInd = 0; sigmaInd < numberOfSigmas; ++sigmaInd ) {
      m_Sigma = sigmas[sigmaInd];
      std::cout << "Processing sigma " << m_Sigma << "..." << std::endl;

      typename GradientImageFilterType::Pointer  gradientImageFilter = GradientImageFilterType::New();
      gradientImageFilter->SetInput( inputImage );
      gradientImageFilter->SetSigma( m_Sigma );
  //gradientImageFilter->SetNormalizeAcrossScale( true );
      gradientImageFilter->ReleaseDataFlagOn();
      gradientImageFilter->Update();

      m_GradientImage = gradientImageFilter->GetOutput();

      m_ScoreImage = InternalImageType::New();
      m_ScoreImage->SetRegions( m_GradientImage->GetRequestedRegion() );
      m_ScoreImage->CopyInformation( m_GradientImage );
      m_ScoreImage->Allocate();
      const float notComputed = -1.0f;
      m_ScoreImage->FillBuffer( notComputed );

      // Set up the multithreaded processing
      ThreadStruct str;
      str.Filter = this;

      unsigned int numThreads = this->GetNumberOfThreads();
      //this->GetMultiThreader()->SetNumberOfThreads( 2 );
      this->GetMultiThreader()->SetNumberOfThreads( numThreads );
      this->GetMultiThreader()->SetSingleMethod( this->ThreaderCallbackScore, &str );

      itk::TimeProbe timer0;
      timer0.Start();
      // multithread the execution
      this->GetMultiThreader()->SingleMethodExecute();
      timer0.Stop();
      std::cout << "Computing score image took " << timer0.GetMeanTime() << " seconds.\n";


      m_FilteredPointsImage = FilteredPointsIndicatorImageType::New();
      m_FilteredPointsImage->SetRegions( m_ScoreImage->GetRequestedRegion() );
      m_FilteredPointsImage->CopyInformation( m_ScoreImage );
      m_FilteredPointsImage->Allocate();
      typename FilteredPointsIndicatorImageType::PixelType const  filteredNotComputed = 0;
      m_FilteredPointsImage->FillBuffer( filteredNotComputed );


      itk::TimeProbe timer1;
      timer1.Start();
      this->LocalContrastFiltering( m_ScoreImage );
      timer1.Stop();
      std::cout << "Computing local contrast filtering took " << timer1.GetMeanTime() << " seconds.\n";

      m_NonMaxedPoints = new std::list< InternalRecord >[numThreads];

      itk::TimeProbe timer2;
      timer2.Start();
      // multithread the execution
      this->GetMultiThreader()->SetSingleMethod( this->ThreaderCallbackNonMax, &str );
      this->GetMultiThreader()->SingleMethodExecute();
      timer2.Stop();
      std::cout << "Nonmax suppression took " << timer2.GetMeanTime() << " seconds.\n";

      for( unsigned int i = 0; i < numThreads; ++i ) {
        m_UnfilteredPoints.splice( m_UnfilteredPoints.begin(), m_NonMaxedPoints[i] );
      }
      delete [] m_NonMaxedPoints;

      std::cout << "Number of points after non-max suppression: " << m_UnfilteredPoints.size() << std::endl;

      m_FilteredPointsImage = 0;  // Free memory.

      #if 0
      // Subvoxel accuracy is done before spatial filtering.
      // This is more work, since some points get discarded but subvoxel accuracy rejects points
      // that are not at a maximum (according to the polynomial fit).
      // Therefore by doing filtering first and then subvoxel accuracy, some points could get filtered
      // because of a nearby point that would get thrown away later for not being at a maximum.
      //
      itk::TimeProbe timer3;
      timer3.Start();
      // filtered indices will contain all points that are a local maximum according to the polynomial fit
      this->SubvoxelAccuracy( m_UnfilteredPoints );
      timer3.Stop();
      std::cout << "Subvoxel accuracy took " << timer3.GetMeanTime() << " seconds.\n";
      #else
      // Subvoxel accuracy is done before spatial filtering.
      // This is more work, since some points get discarded but subvoxel accuracy rejects points
      // that are not at a maximum (according to the polynomial fit).
      // Therefore by doing filtering first and then subvoxel accuracy, some points could get filtered
      // because of a nearby point that would get thrown away later for not being at a maximum.
      //
      itk::TimeProbe timer3;
      timer3.Start();
      // filtered indices will contain all points that are a local maximum according to the polynomial fit
      m_NonMaxedPoints = new std::list< InternalRecord >[numThreads];

      this->GetMultiThreader()->SetSingleMethod( this->ThreaderCallbackSubvoxel, &str );
      this->GetMultiThreader()->SingleMethodExecute();
      timer3.Stop();
      std::cout << "Subvoxel accuracy took " << timer3.GetMeanTime() << " seconds.\n";

      m_UnfilteredPoints.clear();
      for( unsigned int i = 0; i < numThreads; ++i ) {
        //std::copy( m_NonMaxedPoints[i].begin(), m_NonMaxedPoints[i].end(), std::back_inserter( unfilteredPoints ) );
        //for( std::list< InternalRecord >::iterator it = m_NonMaxedPoints[i].begin(); it != m_NonMaxedPoints[i].end(); ++it ) {
        //  unfilteredPoints.push_back( *it );
        //}
        m_UnfilteredPoints.splice( m_UnfilteredPoints.begin(), m_NonMaxedPoints[i] );
      }
      delete [] m_NonMaxedPoints;

      std::cout << "Points that are a local maximum after fit: " << m_UnfilteredPoints.size() << std::endl;

      #endif

      m_ScoreImage = 0;  // Free memory.

      itk::TimeProbe timer6;
      timer6.Start();
      this->ComputeDirections( m_UnfilteredPoints );
      timer6.Stop();
      std::cout << "Improving directions took " << timer6.GetMeanTime() << " seconds.\n";

      m_GradientImage = 0;  // Free memory.
      gradientImageFilter = 0;

      itk::TimeProbe timer4;
      timer4.Start();
      this->ComputeCovariances( m_UnfilteredPoints );
      timer4.Stop();
      std::cout << "Computing covariances took " << timer4.GetMeanTime() << " seconds.\n";

      // recommended hack on releasing memory is to swap in empty vector
      std::list< InternalRecord >  unfilteredEmpty;
      unfilteredEmpty.swap( m_UnfilteredPoints );

    }


  }


  // This indicates that the current BufferedRegion is equal to the
  // requested region. This action prevents useless rexecutions of
  // the pipeline.
  mesh->SetBufferedRegion( mesh->GetRequestedRegion() );

}


namespace {

//: Find the optimum place in quadratic sense
//  Assuming the points are on a parabola and
//  they are sampled on (x0, x1, x2) three points
//  The equation:  a*x^2 + b*x + c = y
//  Now we have three equations and three unknowns
//    a x0^2 + b x0 + c = s0
//    a x1^2 + b x1 + c = s1
//    a x2^2 + b x2 + c = s3
//  So,
//       (y0-y1)(x2-x1) - (y2-y1)(x0-x1)
//   a=  -------------------------------
//       (x0^2-x1^2)(x2-x1) - (x2^2-x1^2)(x0-x1)
//
//       (y2-y1)(x0^2-x1^2) - (y0-y1)(x2^2-x1^2)
//  a =  -------------------------------
//       (x0^2-x1^2)(x2-x1) - (x2^2-x1^2)(x0-x1)
//
//    c = s0 - ax0^2 - bx0
// The optimum pos:
//    xhat= -b / ( 2*a )
//    yhat = (4 a c - b^2 ) / (4*a^2)
static
bool  
find_quad_min( double& max_pos, double& max_strength, 
               double const& y0, double const& y1, double const& y2, 
               double const& x0, double const& x1, double const& x2,
               bool warn=true )
{
  // init 
  max_pos = x1;
  max_strength = y1;

  const double eps = 1e-12;
  const double x0sq = vnl_math_sqr(x0);
  const double x1sq = vnl_math_sqr(x1);
  const double x2sq = vnl_math_sqr(x2);
  
  const double div = (x0sq-x1sq)*(x2-x1) - (x2sq-x1sq)*(x0-x1);
  if( std::abs(div) <=eps ) {
    std::cerr << "Error(" << __FILE__ << ":" << __LINE__ << "): " << std::endl;
    return false;
  }
  const double a = ((y0-y1)*(x2-x1) - (y2-y1)*(x0-x1)) / div;
  const double b = ((y2-y1)*(x0sq-x1sq) - (y0-y1)*(x2sq-x1sq)) / div;
  const double c = y0 - a*x0sq - b*x0;
  
  
  if( std::abs(a) < eps ){
    std::cerr << "Error(" << __FILE__ << ":" << __LINE__ << "): " << std::endl;
    return false;
  }
  // double res1 = a*x1sq + b*x1 + c - y1;
  // double res2 = a*x2sq + b*x2 + c - y2;
  
  const double tmp_max_pos = -b / (2*a);
  const double tmp_max_strength = (4*a*c - b*b) / (4*a);
  
  if( tmp_max_pos < x0 || tmp_max_pos > x2 || tmp_max_strength < y1 )
    if( warn ) {
      std::cerr << "Error(" << __FILE__ << ":" << __LINE__ << "): "
      << " Cannot find maximum pos with points " << x0 << " " << x1 << " " << x2 << std::endl
               << " The scores are: " << y0 << "  " << y1 << "  " << y2 << std::endl;
    return false;
    }
    
  max_pos = tmp_max_pos;
  max_strength = tmp_max_strength;
  return true;
}


// Same as find_quad_min but ( x0, x1, x2 ) is actually ( x0, 0, x2 ).
// This saves a bit of computation.
static
bool  
find_quad_min0( double& max_pos, double& max_strength, 
                double const& y0, double const& y1, double const& y2, 
                double const& x0, double const& x1, double const& x2,
                bool warn=true )
{
  assert( x1 == 0.0 );

  // init 
  max_pos = 0.0;
  max_strength = y1;

  const double eps = 1e-12;
  const double x0sq = vnl_math_sqr(x0);
  const double x2sq = vnl_math_sqr(x2);
  
  const double div = x0sq*x2 - x2sq*x0;
  if( std::abs(div) <=eps ) {
    std::cerr << "Error(" << __FILE__ << ":" << __LINE__ << "): " << std::endl;
    return false;
  }
  const double a = ((y0-y1)*x2 - (y2-y1)*x0) / div;
  const double b = ((y2-y1)*x0sq - (y0-y1)*x2sq) / div;
  const double c = y0 - a*x0sq - b*x0;
  
  
  if( std::abs(a) < eps ){
    std::cerr << "Error(" << __FILE__ << ":" << __LINE__ << "): " << std::endl;
    return false;
  }
  // double res1 = a*x1sq + b*x1 + c - y1;
  // double res2 = a*x2sq + b*x2 + c - y2;
  
  const double tmp_max_pos = -b / (2*a);
  const double tmp_max_strength = (4*a*c - b*b) / (4*a);
  
  if( tmp_max_pos < x0 || tmp_max_pos > x2 || tmp_max_strength < y1 )
    if( warn ) {
      std::cerr << "Error(" << __FILE__ << ":" << __LINE__ << "): "
      << " Cannot find maximum pos with points " << x0 << " " << x1 << " " << x2 << std::endl
               << " The scores are: " << y0 << "  " << y1 << "  " << y2 << std::endl;
    return false;
    }
    
  max_pos = tmp_max_pos;
  max_strength = tmp_max_strength;
  return true;
}

//: Find the optimum place in quadratic sense
//  Assuming the points are on a parabola and
//  they are sampled on (-1, 0, 1) three points
//  The equation:  a*x^2 + b*x + c = s
//  Now we have three equations and three unknowns
//    a - b + c = s0
//            c = s1
//    a + b + c = s3
//  So,
//    a= (s0+s2)/2 - s1
//    b= (s2-s0)/2
//    c= s1
// The optimum pos:
//    x= -b / ( 2*a )
//

template<typename T>
void  
find_maximum_pos( double& max_pos, double& max_strength, 
                  T const& s0, T const& s1, T const& s2, bool warn=true )
{
  // init 
  max_pos = 0;
  max_strength = s1;

  const double b = (double(s2) - double(s0)) / 2;
  const double a2 = (double(s2) + double(s0)) - double(s1)*2;
  const double c = s1;
  
  if( a2 < -vnl_math::eps )
    max_pos = - b / (a2);
  
  // if it is not maximum, or outside of range [-1, 1]
  if( a2>=0 || max_pos > 1 || max_pos < -1 ) 
  { 
   max_pos = 0.0;
    //debug output
    if( warn )
      std::cerr << "Error(" << __FILE__ << ":" << __LINE__ << "): "
               << " Cannot find maximum pos within [-1,1].\n"
               << " The scores are: " << s0 << "  " << s1 << "  " << s2 << std::endl;
  } 
  else {
    
    max_strength = -b*b/(a2*2) + c;
  }
  
}

}


template <class TInputImage, class TOutputMesh>
bool
FeatureImageFilter< TInputImage, TOutputMesh>
::SubvoxelLocalizationLS( InternalRecord &  result,
                          typename InterpolatorType::Pointer const &  interpolator,
                          itk::Matrix< float, 10, 10 > const &  LS )
{
  const InternalImageType*  scoreImage = interpolator->GetInputImage();
  //const InterpolatorType::InputImageType*  scoreImage = interpolator->GetInputImage();

  typedef typename InternalImageType::IndexType  IndexType;

  IndexType const & pointIndex = result.m_LocationIndex;
  PhysicalPointType const &  pointPhysical = result.m_Location;

  typedef ConstNeighborhoodIterator< InternalImageType >  IteratorType;
  IteratorType  bit;
  typename IteratorType::RadiusType  radius = { 1, 1, 1 };
  bit = IteratorType( radius, scoreImage, scoreImage->GetRequestedRegion() );

  bit.SetLocation( pointIndex );

  typedef typename PhysicalPointType::CoordRepType  CoordRepType;
  itk::Vector< float, 10 >  Y;
  Y.Fill( 0.0f );
  for( unsigned int i = 0; i < bit.Size(); ++i )
    {
      InputIndexType  XYZindex = bit.GetIndex(i);

      CoordRepType Xi = static_cast<CoordRepType>( XYZindex[0] ) -  static_cast<CoordRepType>( pointIndex[0] );
      CoordRepType Yi = static_cast<CoordRepType>( XYZindex[1] ) -  static_cast<CoordRepType>( pointIndex[1] );
      CoordRepType Zi = static_cast<CoordRepType>( XYZindex[2] ) -  static_cast<CoordRepType>( pointIndex[2] );

#define UNIT_NORMAL_STEP 1

#if UNIT_NORMAL_STEP
      // abuse of the neighborhood iterator, treat the locations as if they were in the physical coordinates
      itk::Vector< float, 3 >  step;
      step[0] = Xi;
      step[1] = Yi;
      step[2] = Zi;
      PhysicalPointType  afterStep = pointPhysical + step;
      typename InternalImageType::PixelType  Ii = interpolator->Evaluate( afterStep );
#else
      typename InternalImageType::PixelType  Ii = bit.GetPixel(i);
#endif

      //if( (i % 3) == 0 ) std::cout << std::endl;
      //if( (i % 9) == 0 ) std::cout << std::endl;
      ////std::cout << "Xi: " << Xi << " Yi: " << Yi << " Zi: " << Zi << std::endl;
      //std::cout << Ii << " \t";

      Y[0] += Ii;
      Y[1] += Ii * Zi;
      Y[2] += Ii * Yi;
      Y[3] += Ii * Xi;
      Y[4] += Ii * Zi * Zi;
      Y[5] += Ii * Yi * Zi;
      Y[6] += Ii * Xi * Zi;
      Y[7] += Ii * Yi * Yi;
      Y[8] += Ii * Xi * Yi;
      Y[9] += Ii * Xi * Xi;
    }

  // parameters of the parabola
  itk::Vector< float, 10 >  X = LS * Y;

  itk::Matrix< float, 3, 3 >  A;
  A[0][0] = 2*X[0];
  A[0][1] =   X[1];
  A[0][2] =   X[3];
  
  A[1][1] = 2*X[2];
  A[1][2] =   X[4];
  
  A[2][2] = 2*X[5];
  
  A[1][0] = A[0][1];
  A[2][0] = A[0][2];
  A[2][1] = A[1][2];

  // The following rule can be rewritten in terms of determinants -> less expensive
  //
  // positive definite: all major subdeterminants are positive
  // negative definite: all major subdeterminants change signs, beginning with negative
  typedef itk::SymmetricSecondRankTensor< double, 3 >  SymmetricMatrixType;
  SymmetricMatrixType  AforEigen;
  for( unsigned int i = 0; i < 3; ++i )
    for( unsigned int j = 0; j < 3; ++j )
      AforEigen( i, j ) = A[i][j];
  SymmetricMatrixType::EigenValuesArrayType  eigenValues;
  SymmetricMatrixType::EigenVectorsMatrixType  eigenVectors;
  AforEigen.ComputeEigenAnalysis( eigenValues, eigenVectors );
  //for( unsigned int i = 0; i < 3; ++i )
  //  if( eigenValues[i] >= 0 ) {
  //    //std::cout << "Warning: Current point is not a maximum (matrix not negative definite)." << std::endl;
  //    return false;
  //  }
  if( eigenValues[0] >= 0 && eigenValues[1] >= 0 && eigenValues[2] >= 0 ) {
      //std::cout << "Warning: Current point is not a maximum (matrix not negative definite)." << std::endl;
      return false;
  }

  itk::Vector< double, 3 >  b;
  b[0] = - X[6];
  b[1] = - X[7];
  b[2] = - X[8];

  // solving for location, given parameters
  //itk::Matrix< float, 3, 3 >  Ainv = A.GetInverse();
  SymmetricMatrixType::EigenVectorsMatrixType  Dinv;
  Dinv.Fill( 0.0 );
  for( unsigned int i = 0; i < 3; ++i ) {
    if( eigenValues[i] > -1.0 ) eigenValues[i] = 0;   // this direction is not a maximum  (this should be > 0.0, but we don't care about small e-values)
    else Dinv[i][i] = 1.0 / eigenValues[i];
  }
	// this casting needs to be here since the explicit constructor of itkMatrix would not get called otherwise
  itk::Matrix< double, 3, 3 >  eigenVectorsT = itk::Matrix< double, 3, 3 >( eigenVectors.GetTranspose() );
  itk::Matrix< double, 3, 3 >  Ainv = eigenVectors * Dinv * eigenVectorsT;
  itk::Vector< double, 3 >  refinement = Ainv * b;

#if UNIT_NORMAL_STEP
  if( refinement.GetNorm() > 1.73 ) {  // sqrt( 3 )
    //std::cout << "Warning: refinement step too large: " << refinement.GetNorm() << std::endl;
    return false;
  }
  else
  {
    result.m_Location += refinement;
    //std::cout << refinement.GetNorm() << " " << (pointPhysical-refinedPoint).GetNorm() << std::endl;
  }
# else
  itk::ContinuousIndex< float, InputImageDimension >  refinedIndex;
  //if( refinement.GetNorm() > 1.73 ) {  // sqrt( 3 )
  //  return false;
  //}
  //else
  {
    refinedIndex[0] = pointIndex[0] + refinement[0];
    refinedIndex[1] = pointIndex[1] + refinement[1];
    refinedIndex[2] = pointIndex[2] + refinement[2];
  }

  //refinedIndex[0] = pointIndex[0];
  //refinedIndex[1] = pointIndex[1];
  //refinedIndex[2] = pointIndex[2];
  //if( eigenValues[0] < 0 ) refinedIndex[0] += refinement[0];
  //if( eigenValues[1] < 0 ) refinedIndex[1] += refinement[1];
  //if( eigenValues[2] < 0 ) refinedIndex[2] += refinement[2];
  scoreImage->TransformContinuousIndexToPhysicalPoint( refinedIndex, refinedPoint );
#endif

  float x = refinement[0];
  float y = refinement[1];
  float z = refinement[2];

  // compute the value at the interpolated location (=interpolated score)
  itk::Vector< float, 10 >  P;
  P[9] =   1.0;
  P[8] =     z;
  P[7] =     y;
  P[6] =     x;
  P[5] = z * z;
  P[4] = y * z;
  P[3] = x * z;
  P[2] = y * y;
  P[1] = x * y;
  P[0] = x * x;

  //float beforeInterpStrength = result.m_Strength;
  result.m_Strength = dot_product( P.GetVnlVector(), X.GetVnlVector() );
//if( result.m_Strength < 0.0 ) std::cout << "before: " << beforeInterpStrength << " after: " << result.m_Strength << std::endl;
if( result.m_Strength < 0.0 ) return false;// because of the hack during score image computation, some -1.0 (boundary) score values are included

  //std::cout << P.GetVnlVector() << std::endl << X.GetVnlVector() << std::endl;

  return true;
}


/**
 * Standard "PrintSelf" method
 */
template <class TInputImage, class TOutputMesh>
void
FeatureImageFilter< TInputImage, TOutputMesh>
::PrintSelf(
  std::ostream& os,
  Indent indent) const
{
  Superclass::PrintSelf( os, indent );

}

} // end namespace itk





#endif
