/*=========================================================================


author: Michal Sofka, Rensselaer Polytechnic Institute
  date: 10/10/2007
  

=========================================================================*/
#ifndef _itkDescriptorMeshFilter_txx
#define _itkDescriptorMeshFilter_txx
#include "itkDescriptorMeshFilter.h"

#include "itkVectorLinearInterpolateImageFunction.h"
#include "itkConstNeighborhoodIterator.h"

#include "itkImageRegionIterator.h"
#include "itkNeighborhoodAlgorithm.h"
#include "itkZeroFluxNeumannBoundaryCondition.h"
#include "itkConstantBoundaryCondition.h"

#include "itkGaussianDerivativeImageFunction.h"
#include "itkDerivativeImageFilter.h"
//#include "itkGradientImageFilter.h"
#include "itkGaussianDerivativeImageFunction.h"
#include "itkRecursiveGaussianImageFilter.h"
#include "itkMultiResolutionPyramidImageFilter.h"

#include "itkMeshSource.h"

#include "itkFeatureImageFilter.h"


#include <vnl/algo/vnl_determinant.h>
#include <vnl/vnl_trace.h>
#include "rsdl_bins.h"

// for debug
#include "itkNumericSeriesFileNames.h"
#include "itkImageSeriesWriter.h"
#include "itkRescaleIntensityImageFilter.h"

#include "itkTimeProbe.h"

#define DEBUG 0

#if DEBUG
#define PRINTF(...)  printf(__VA_ARGS__)
#else
#define PRINTF(...)  0
#endif

//#if DEBUG
//#define PRINTF1( text, arg1 ) printf( text, arg1 )
//#define PRINTF2( text, arg1, arg2 ) printf( text, arg1, arg2 )
//#define PRINTF3( text, arg1, arg2, arg3 ) printf( text, arg1, arg2, arg3 )
//#else
//#define PRINTF1( text, arg1 ) 0
//#define PRINTF2( text, arg1, arg2 ) 0
//#define PRINTF3( text, arg1, arg2, arg3 ) 0
//#endif

namespace itk
{


template <class TInputMesh1, class TInputMesh2, class TOutputMesh>
const double DescriptorMeshFilter< TInputMesh1, TInputMesh2, TOutputMesh>::m_2pi = 2.0 * vnl_math::pi; //6.28318530717959;


template <class TInputMesh1, class TInputMesh2, class TOutputMesh>
DescriptorMeshFilter< TInputMesh1, TInputMesh2, TOutputMesh>
::DescriptorMeshFilter()
{
  // Modify superclass default values, can be overridden by subclasses
  this->SetNumberOfRequiredInputs(2);

  //PointDataContainerPointer  pointData  = PointDataContainer::New();
  //OutputMeshPointer  mesh = this->GetOutput();
  //mesh->SetPointData( pointData.GetPointer() );

  m_Radius = 30.0f;

  // number of bins in each coordinate
  // flipping the angles w.r.t azimuth and elevation -> only half of bins needed
  m_RadiusNum = 5;
  m_OrientNum = 4;

  double logRadius = vcl_log( m_Radius );
  
  // bin sizes
  m_LogRadiusBinSize = logRadius / m_RadiusNum;  // on the logarithmic scale
  m_AzimuthBinSize   = vnl_math::pi / m_OrientNum;
  m_ElevationBinSize = vnl_math::pi / m_OrientNum;

  m_KeypointsWithDescriptors = 0;
  m_KeypointsWithDescriptorsData = 0;
  m_kdTreeFeatures = 0;
}


template <class TInputMesh1, class TInputMesh2, class TOutputMesh>
void
DescriptorMeshFilter< TInputMesh1, TInputMesh2, TOutputMesh>
::GenerateOutputInformation()
{
}


template <class TInputMesh1, class TInputMesh2, class TOutputMesh>
void
DescriptorMeshFilter< TInputMesh1, TInputMesh2, TOutputMesh>
::SetInputKeypoints( const TInputMesh1 * input )
{

  // This const_cast is needed due to the lack of
  // const-correctness in the ProcessObject.
  this->SetNthInput( 0,
            const_cast< TInputMesh1 * >( input ) );

}


template <class TInputMesh1, class TInputMesh2, class TOutputMesh>
void
DescriptorMeshFilter< TInputMesh1, TInputMesh2, TOutputMesh>
::SetInputFeatures( const TInputMesh2 * input )
{

  // This const_cast is needed due to the lack of
  // const-correctness in the ProcessObject.
  this->SetNthInput( 1,
            const_cast< TInputMesh2 * >( input ) );

}


template <class TInputMesh1, class TInputMesh2, class TOutputMesh>
int 
DescriptorMeshFilter< TInputMesh1, TInputMesh2, TOutputMesh>
::SplitRequestedList( int i, int num,
                      typename KeypointIteratorType &  startIterator,
                      typename KeypointIteratorType &  endIterator )
{
  typename KeypointMeshType::Pointer       keypoints    = dynamic_cast< TInputMesh1 * >( this->GetInput(0) );
  typename KeypointMeshType::PointsContainerPointer      keypointPoints    = keypoints->GetPoints();

  // determine the actual number of pieces that will be generated
  unsigned long  range = keypointPoints->Size();
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
  startIterator = keypointPoints->begin();
  endIterator = keypointPoints->begin();
  for( unsigned int ind = 0; ind < startIndex; ++ind ) ++startIterator;
  for( unsigned int ind = 0; ind < endIndex; ++ind ) ++endIterator;

  itkDebugMacro("  Split Piece: " << startIndex << " " << endIndex );
  //std::cout << "  Split Piece: " << startIndex << " " << endIndex << std::endl;

  return maxThreadIdUsed + 1;
}


template <class TInputMesh1, class TInputMesh2, class TOutputMesh>
void
DescriptorMeshFilter< TInputMesh1, TInputMesh2, TOutputMesh>
::AddToBin( itk::Vector< float, 3 > const &  norm_vector,
            double const &  keypoint_azimuth,
            double const &  keypoint_elevation,
            unsigned int ind,
            std::vector< std::vector< std::vector< vnl_vector_fixed< float, 3 > > > > &  tempDescriptor,
            itk::Vector< float, 3 > const &  keypointNormal )
{
  typename FeatureMeshType::Pointer        features     = dynamic_cast< TInputMesh2 * >( this->GetInput(1) );
  typename FeatureMeshType::PointDataContainerPointer      featureData    = features->GetPointData();


  /*               (NOT DRAWN TO SCALE)

                       \-----|-----/
                      / \ 16 |  9 / \  overflow area
                     /   \___|___/   \  
                    /	15 /\8 |1 /\ 10 \
                    |   /7 \_|_/ 2\   |
                    ----|--|_C_|--|----    
                    |   \6 /5|4\ 3/   |
                    \	14 \/__|__\/ 11 /
                     \   /   |   \   /
                      \ / 13 | 12 \ /
                       /-----|-----\
                             _
                            | |
                           _| |_
                           \   /
                            \ /
                             V

             +---+---+---+---+---+---+---+---+
           --> c | c | c | c | c | c | c | c -->         (C = sum(c_i))
             +---+---+---+---+---+---+---+---+
           --> 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 -->   
             +---+---+---+---+---+---+---+---+
           --> 9 | 10| 11| 12| 13| 14| 15| 16-->
             +---+---+---+---+---+---+---+---+
                     overflow
             +---+---+---+---+---+---+---+---+
       */


  //timerATBprepare.Start();
  double radius = norm_vector.GetNorm();
  if (radius < 1)  // apparently feature points can be less than 1 pixel away from the keypoint.
    radius = 1;  // if radius is < 1, then the log is negative, and it indexes a negative bin location.
  double log_radius = vcl_log( radius );


  // angle of the current context vector w.r.t. keypoint angle
  double azimuth = vcl_atan2( norm_vector[1], norm_vector[0] ) + vnl_math::pi;// add pi so that it's in the range <0,2*pi>
  if (azimuth >= m_2pi) azimuth -= m_2pi; // because 6.2831853946023664 / pi = 2.00000003. atan2 floating point error.
  double azimuth_wrt_keypoint = azimuth - keypoint_azimuth;
  if( azimuth_wrt_keypoint < 0.0 ) azimuth_wrt_keypoint += m_2pi;
  if( azimuth_wrt_keypoint >= vnl_math::pi ) azimuth_wrt_keypoint -= vnl_math::pi;


  double elevation = vcl_atan2( norm_vector[2], norm_vector[0] ) + vnl_math::pi;// add pi so that it's in the range <0,2*pi>
  if (elevation >= m_2pi) elevation -= m_2pi; // because 6.2831853946023664 / pi = 2.00000003. atan2 floating point error.
  double elevation_wrt_keypoint = elevation - keypoint_elevation;
  if( elevation_wrt_keypoint < 0.0 ) elevation_wrt_keypoint += m_2pi;
  if( elevation_wrt_keypoint >= vnl_math::pi ) elevation_wrt_keypoint -= vnl_math::pi;

  // feature point normal
  vnl_vector_fixed< float, 3 >  context_normal = (*featureData)[ind].m_Directions[0].GetVnlVector();
  // limit the normal to be in the half space
  if( context_normal[2] < 0.0 )  context_normal *= -1.0;


  double azFrac = azimuth_wrt_keypoint / m_AzimuthBinSize;
  unsigned int azimuthBin = (unsigned int)( vcl_floor( azFrac ) ) % m_OrientNum;  // module handles the case when azFrac == m_OrientNum
  azFrac -= (azimuthBin-.5);

  double elFrac = elevation_wrt_keypoint / m_ElevationBinSize;
  unsigned int elevationBin = (unsigned int)( vcl_floor( elFrac ) ) % m_OrientNum; // module handles the case when elFrac == m_OrientNum
  elFrac -= (elevationBin-.5);

  double radFrac = log_radius/m_LogRadiusBinSize;
  unsigned int radiusBin = (unsigned int)( vcl_floor(radFrac) );
  radFrac -= (radiusBin-.5);

/* if angleFrac > 1, add (frac-1)*stuff to next bin, else add (1-frac)*stuff to prev bin

          Example:

            2       3
        /\_/^\_/\   |
        +-------+---:---+-------+    * a point at 3, if the bin size is 2,
        |       |   :   |       |      falls in the middle of bin 1, because 3/2 = 1.5
        |       |...*...|       |     
        |       |   :   |       |      The percentage contribution to the bin is 100%, because
        +-------+---:---+-------+      1.5-.5 = 1
    bin:    0       1       2        
                                       A point to the right of center in a bin will get a frac
                                       in the range (1 , 1.5). In this case the contribution is 2-frac.
         +-----|-----+                              
         |           |                 A point to the left of center in a bin will get a frac
         |           |                 in the range (.5, 1). In this case the contribution is frac.
         -     *     -
         |           |                 Whatever percentage is remaining gets contributed to the next bin over.
         |           |
         +-----|-----+
        /      |      \
      0.5      1     1.499... 
    frac values for locations in bins
  */


  // only two bins (along each dimension) will get a contribution, the other is zero
  //
  // previous (or the next) azimuth bin which gets contribution; the other bin gets zero contribution
  unsigned int prevNextAzBin;
  if( azFrac > 1.0 ) {
    azFrac = 2.0 - azFrac;
    // prevNextAzBin = ( azimuthBin + 1 ) % m_OrientNum;
    prevNextAzBin = (azimuthBin < m_OrientNum-1) ? azimuthBin+1 : 0;   // next_azimuth_bin
  }
  else {
    // prevNextAzBin = ( azimuthBin + m_OrientNum - 1 ) % m_OrientNum;
    prevNextAzBin = (azimuthBin > 0) ? azimuthBin-1 : m_OrientNum-1;   // prev_azimuth_bin 
  }
  double prevNextAzFrac = 1.0 - azFrac;


  // previous (or the next) elevation bin which gets contribution; the other bin gets zero contribution
  unsigned int prevNextElBin;
  if( elFrac > 1.0 ) {
    elFrac = 2.0 - elFrac;
    // prevNextElBin = ( elevationBin + 1 ) % m_OrientNum;
    prevNextElBin = (elevationBin < m_OrientNum-1) ? elevationBin+1 : 0;   // next_elevation_bin
  }
  else {
    // prevNextElBin = ( elevationBin + m_OrientNum - 1 ) % m_OrientNum;
    prevNextElBin = (elevationBin > 0) ? elevationBin-1 : m_OrientNum-1;   // prev_elevation_bin
  }
  double prevNextElFrac = 1.0 - elFrac;

  
  // previous (or the next) radius bin which gets contribution; the other bin gets zero contribution
  unsigned int prevNextRadBin;
  if( radFrac > 1.0 ) {
    radFrac = 2.0 - radFrac;
    prevNextRadBin = radiusBin+1;   // next_radius_bin
  }
  else {
    prevNextRadBin = (radiusBin > 0) ? radiusBin-1 : 0;   // prev_radius_bin
  }
  double prevNextRadFrac = 1.0 - radFrac;
  //timerATBprepare.Stop();

  //timerATBcompute.Start();
  // fill in the ( 2x2x2 = 8 ) bins
  tempDescriptor[radiusBin][azimuthBin][elevationBin]          += float( radFrac * azFrac * elFrac ) * context_normal;
  tempDescriptor[radiusBin][azimuthBin][prevNextElBin]         += float( radFrac * azFrac * prevNextElFrac ) * context_normal;
  tempDescriptor[radiusBin][prevNextAzBin][elevationBin]       += float( radFrac * prevNextAzFrac * elFrac ) * context_normal;
  tempDescriptor[radiusBin][prevNextAzBin][prevNextElBin]      += float( radFrac * prevNextAzFrac * prevNextElFrac ) * context_normal;
  tempDescriptor[prevNextRadBin][azimuthBin][elevationBin]     += float( prevNextRadFrac * azFrac * elFrac ) * context_normal;
  tempDescriptor[prevNextRadBin][azimuthBin][prevNextElBin]    += float( prevNextRadFrac * azFrac * prevNextElFrac ) * context_normal;
  tempDescriptor[prevNextRadBin][prevNextAzBin][elevationBin]  += float( prevNextRadFrac * prevNextAzFrac * elFrac ) * context_normal;
  tempDescriptor[prevNextRadBin][prevNextAzBin][prevNextElBin] += float( prevNextRadFrac * prevNextAzFrac * prevNextElFrac ) * context_normal;
  //timerATBcompute.Stop();


  //std::cout << prev_radFrac << "\t" <<  radFrac << "\t" <<  next_radFrac << "\t" <<  prev_azFrac << "\t" <<  azFrac << "\t" <<  next_azFrac << "\t" <<  prev_elFrac << "\t" <<  elFrac << "\t" <<  next_elFrac << std::endl;
  //std::cout << "r:" << radius << " log_r: " << log_radius << " rBin: " << radiusBin << vcl_endl;
  //std::cout << "r:" << radius << " log_r: " << log_radius << " a:" << angle << " " << " ind: " << descriptor_index << " rBin: " << radiusBin << " aBin: " << angleBin << vcl_endl;
  //std::cout << "rbin_size: " << m_LogRadiusBinSize << " abin_size: " << angleBin_size << " rBin: " << radiusBin << " aBin: " << angleBin << vcl_endl;
}


template <class TInputMesh1, class TInputMesh2, class TOutputMesh>
void
DescriptorMeshFilter< TInputMesh1, TInputMesh2, TOutputMesh>
::CopyDescriptor( std::vector< std::vector< std::vector< vnl_vector_fixed< float, 3 > > > > &  tempDescriptor,
                  vnl_vector< float > &  descriptor )
{
  vnl_vector_fixed< float, 3 >  centerBin( 0.0 );
  for (unsigned int j = 0; j < m_OrientNum; j++)
    for (unsigned int k = 0; k < m_OrientNum; k++)
      centerBin += tempDescriptor[0][j][k];

  // normalize the center bin by its size
  double centerBin_radius = vcl_exp( m_LogRadiusBinSize );
  double wedge_volume = 4.0 / 3.0 * 3.1415 * vcl_pow( centerBin_radius, 3.0 );
  double wedge_volume_cube_root = vcl_pow( wedge_volume, 1.0/3.0 );
  // according to frome:eccv04, a good way of normalizing is by a cube root of the bin volume
  centerBin *= 1.0 / wedge_volume_cube_root;

  // initialize descriptor vector
  descriptor.set_size( (1 + (m_RadiusNum-1) * m_OrientNum * m_OrientNum) * 3 );
  descriptor.fill( 0.0 );

  descriptor[0] = centerBin[0];
  descriptor[1] = centerBin[1];
  descriptor[2] = centerBin[2];

  unsigned int loc = 3;
  double constant = 4.0 / 3.0 * 3.1415 / ( m_OrientNum * m_OrientNum );
  for( unsigned int i = 1; i < m_RadiusNum; i++) {
    double radiusBin = i;
    for( unsigned int j = 0; j < m_OrientNum; j++) {
      for( unsigned int k = 0; k < m_OrientNum; k++) {
        double larger_radius = vcl_exp( (radiusBin+1)*m_LogRadiusBinSize );
        double smaller_radius = vcl_exp( radiusBin*m_LogRadiusBinSize );

        // up to a constant: pi / m_OrientNum / m_OrientNum
        double wedge_volume = constant * ( vcl_pow( larger_radius, 3.0 ) - vcl_pow( smaller_radius, 3.0 ) );

        // according to frome:eccv04, a good way of normalizing is by a cube root of the bin volume
        double wedge_volume_cube_root = vcl_pow( wedge_volume, 1.0/3.0 );

        descriptor[loc]   = tempDescriptor[i][j][k][0] / wedge_volume_cube_root;
        descriptor[loc+1] = tempDescriptor[i][j][k][1] / wedge_volume_cube_root;
        descriptor[loc+2] = tempDescriptor[i][j][k][2] / wedge_volume_cube_root;
        loc += 3;
      }
    }
  }


  // Normalize the descriptor vector.
  double  two_norm = descriptor.two_norm();
  if( two_norm > 0 ) descriptor /= two_norm;

}


template <class TInputMesh1, class TInputMesh2, class TOutputMesh>
void
DescriptorMeshFilter< TInputMesh1, TInputMesh2, TOutputMesh>
::GenerateData()
{
  typename KeypointMeshType::Pointer       keypoints    = dynamic_cast< TInputMesh1 * >( this->GetInput(0) );
  typename FeatureMeshType::Pointer        features     = dynamic_cast< TInputMesh2 * >( this->GetInput(1) );

  PointsContainerPointer      points    = keypoints->GetPoints();
  //PointDataContainerPointer   pointData = keypoints->GetPointData();

  // creating new points, because some keypoints are discarded
  // (if there is not enough feature points to compute the descriptor or if the orientation cannot be computed reliably)
  PointsContainerPointer      pointsNew = PointsContainer::New();
  PointDataContainerPointer   pointData = PointDataContainer::New();

  OutputMeshPointer           mesh      = this->GetOutput();
  mesh->SetPoints( pointsNew );
  mesh->SetPointData( pointData );


  // timers
  itk::TimeProbe  timerBuildingTree;
  itk::TimeProbe  timerGettingFeatures;
  itk::TimeProbe  timerPreparingVectors;
  itk::TimeProbe  timerAddingToBins;
  itk::TimeProbe  timerCopyingDescriptor;
  itk::TimeProbe  timerComputingDescriptors;


  // compute a feature kd-tree for fast extraction of points within the descriptor neighborhood

  // setup kd-tree for locations (these will be used to gather major gradient orientation)
  //
  timerComputingDescriptors.Start();
  timerBuildingTree.Start();

  FeatureSampleType::Pointer    sampleFeatures = FeatureSampleType::New();


  // feature locations
  //
  typename FeatureMeshType::PointsContainerPointer      featurePoints    = features->GetPoints();
  pointsNew->reserve( featurePoints->size() );
  pointData->reserve( featurePoints->size() );
  for( unsigned int n = 0 ; n < featurePoints->size(); ++n ) {
    typename FeatureMeshType::PointType  featurePoint;
    features->GetPoint( n, &featurePoint );

    LocationVectorType  mv_k;
    mv_k[0] = featurePoint[0];
    mv_k[1] = featurePoint[1];
    mv_k[2] = featurePoint[2];
    //mv.Set_vnl_vector( m_FixedKeypoints[n]->location_.as_vector() );
    sampleFeatures->PushBack( mv_k );
  }


  // feature locations
  typename FeatureTreeGeneratorType::Pointer  featureTreeGenerator = FeatureTreeGeneratorType::New();

  featureTreeGenerator->SetSample( sampleFeatures );
  featureTreeGenerator->SetBucketSize( 16 );
  featureTreeGenerator->Update();

  m_kdTreeFeatures = featureTreeGenerator->GetOutput();
  timerBuildingTree.Stop();


  unsigned int numThreads = this->GetNumberOfThreads();

  m_KeypointsWithDescriptors = new std::list< KeypointMeshType::PointType >[numThreads];
  m_KeypointsWithDescriptorsData = new std::list< DescriptorAttribute >[numThreads];
  

  // Set up the multithreaded processing
  ThreadStruct str;
  str.Filter = this;

  this->GetMultiThreader()->SetNumberOfThreads( numThreads );
  this->GetMultiThreader()->SetSingleMethod( this->ThreaderCallback, &str );

  itk::TimeProbe timer0;
  timer0.Start();
  // multithread the execution
  this->GetMultiThreader()->SingleMethodExecute();
  timer0.Stop();
  std::cout << "Computing descriptors threaded took " << timer0.GetMeanTime() << " seconds.\n";


  for( unsigned int i = 0; i < numThreads; ++i ) {
    //m_UnfilteredPoints.splice( m_UnfilteredPoints.begin(), m_NonMaxedPoints[i] );
    std::list< DescriptorAttribute >::iterator itData = m_KeypointsWithDescriptorsData[i].begin();
    for( std::list< KeypointMeshType::PointType >::iterator  it = m_KeypointsWithDescriptors[i].begin(); it != m_KeypointsWithDescriptors[i].end(); ++it, ++itData ) {
      pointsNew->push_back( *it );
      pointData->push_back( *itData );
    }
  }

  delete [] m_KeypointsWithDescriptors;
  delete [] m_KeypointsWithDescriptorsData;



  timerComputingDescriptors.Stop();




  // This indicates that the current BufferedRegion is equal to the
  // requested region. This action prevents useless rexecutions of
  // the pipeline.
  mesh->SetBufferedRegion( mesh->GetRequestedRegion() );

}

#define USE_TIMER 0
#if USE_TIMER
#define TIMER( arg ) arg;
#else
#define TIMER( arg ) 0
#endif

template <class TInputMesh1, class TInputMesh2, class TOutputMesh>
void
DescriptorMeshFilter< TInputMesh1, TInputMesh2, TOutputMesh>
::ComputeDescriptorsThreaded( KeypointIteratorType &  startIterator,
                              KeypointIteratorType &  endIterator,
                              int threadId )
{
  typename FeatureMeshType::Pointer        features     = dynamic_cast< TInputMesh2 * >( this->GetInput(1) );
  typename FeatureMeshType::PointDataContainerPointer      featureData    = features->GetPointData();

  typename FeatureMeshType::PointsContainerPointer      featurePoints    = features->GetPoints();

#if USE_TIMER
  // timers
  itk::TimeProbe  timerBuildingTree;
  itk::TimeProbe  timerGettingFeatures;
  itk::TimeProbe  timerPreparingVectors;
  itk::TimeProbe  timerAddingToBins;
  itk::TimeProbe  timerCopyingDescriptor;
  itk::TimeProbe  timerComputingDescriptors;
#endif

  // not enough points to compute the descriptor for this keypoint
  unsigned int notEnoughPoints = 0;

  // cannot compute orientation for this keypoint
  unsigned int badOrientation = 0;

  for( KeypointIteratorType  it = startIterator; it != endIterator; ++it ) {
    // gather features in a neighborhood around keypoint
    //
    typename KeypointMeshType::PointType const &  keypointPoint = *it;

    TIMER( timerGettingFeatures.Start() );
    LocationVectorType  locationVector;
    locationVector[0] = keypointPoint[0];
    locationVector[1] = keypointPoint[1];
    locationVector[2] = keypointPoint[2];
    typename FeatureTreeType::InstanceIdentifierVectorType  featureNeighbors;
    m_kdTreeFeatures->Search( locationVector, m_Radius, featureNeighbors );
    TIMER( timerGettingFeatures.Stop() );

    // not enough points to compute the descriptor
    if( featureNeighbors.size() < 100 ) {
      //std::cout << "Not enough points to compute the descriptor." << std::endl;
      ++notEnoughPoints;
      continue;
    }

    //FeatureAtIndexIsGreater  sortOnDistances( featuresAll, filtered_features[i] );
    //vcl_sort( featureNeighbors.begin(), featureNeighbors.end(), sortOnDistances );


    // compute context vectors
    // context vectors are vectors created from feature point locations within the current shape context radius
    // and a keypoint location for which the descriptor is computed
    //
    std::vector< itk::Vector< float, 3 > >  context_vectors;
    context_vectors.reserve( featureNeighbors.size() );
    //
    // at the same time orient the keypoint based on the average of feature orientations
    //
    TIMER( timerPreparingVectors.Start() );
    itk::Vector< float, 3 >  keypoint_direction;
    keypoint_direction.Fill( 0.0 );
    itk::Vector< float, 3 >  keypoint_bidirection;
    keypoint_bidirection.Fill( 0.0 );
    double sum_w = 0.0;
    for( std::vector< unsigned int >::size_type  n = 0; n < featureNeighbors.size(); ++n )
      {
      const unsigned int ind = featureNeighbors[n];
      typename FeatureMeshType::PointType const &  featurePoint = (*featurePoints)[ind];

      // compute context vector
      //FeatureMeshType::PointType  context_vector = keypointPoint - featurePoint;
      itk::Vector< float, 3 >  context_vector;
      context_vector[0] = keypointPoint[0] - featurePoint[0];
      context_vector[1] = keypointPoint[1] - featurePoint[1];
      context_vector[2] = keypointPoint[2] - featurePoint[2];

      context_vectors.push_back( context_vector );


      // aggregate feature orientation
      //
      PointAttribute const &  featurePointData = (*featureData)[ind];

      if( featurePointData.m_Shape == PointAttribute::Sheet ) continue;
 
      keypoint_direction += featurePointData.m_Strength * featurePointData.m_Directions[0];
      sum_w += featurePointData.m_Strength;

      keypoint_bidirection += featurePointData.m_Strength * featurePointData.m_Directions[1];

      }

    if( sum_w == 0.0 ) {
      //std::cout << "Features for the descriptor are all sheets -> can't compute keypoint orientation." << std::endl;
      ++badOrientation;
      continue;
    }

    keypoint_direction *= ( 1.0 / sum_w );
    keypoint_direction.Normalize();

    keypoint_bidirection *= ( 1.0 / sum_w );
    keypoint_bidirection.Normalize();




    itk::Vector< float, 3 >  keypointNormal = keypoint_direction;
    itk::Vector< float, 3 >  keypointBinormal = keypoint_bidirection;

    // limiting to half space for keypoint normal (ellipsoid is mirrored and from the keypoint point of view
    // the two halves have the same orientation) -> flip the normal
    if( keypointNormal[2] < 0.0 ) keypointNormal *= -1.0;

    // limiting to half space for keypoint binormal (ellipsoid is mirrored and from the keypoint point of view
    // the two halves have the same orientation) -> flip the binormal
    if( keypointBinormal[1] < 0.0 ) keypointBinormal *= -1.0;

    double keypoint_azimuth = vcl_atan2( keypointNormal[1], keypointNormal[0] );
    double keypoint_elevation = vcl_atan2( keypointNormal[2], keypointNormal[0] );
    keypoint_azimuth += vnl_math::pi; //add pi so that it's in the range <0,2*pi>
    keypoint_elevation += vnl_math::pi; //add pi so that it's in the range <0,2*pi>

    vnl_matrix_fixed< float, 3, 3 >  rotation;
    rotation.set_column( 0, keypointNormal.GetVnlVector() );
    rotation.set_column( 1, keypointBinormal.GetVnlVector() );
    // ensure that the third vector forms a right hand coordinate system

    // ARE WE SAVING ALL DIRECTIONS IN m_Diretions? (if so, the following is not necessary)
    vnl_vector_fixed< float, 3 >  third_vector = vnl_cross_3d( keypointNormal.GetVnlVector(), keypointBinormal.GetVnlVector() );
    //std::cout << third_vector << vcl_endl << third_from_cross << vcl_endl << dot_product( third_vector, third_from_cross ) << vcl_endl;
    rotation.set_column( 2, third_vector );


    std::vector< std::vector< std::vector< vnl_vector_fixed< float, 3 > > > >  tempDescriptor(m_RadiusNum+1); //+1 to catch overflow off the outside edge
    // force initialization to 0
    for (unsigned int z = 0; z < m_RadiusNum+1; z++) {
      tempDescriptor[z].resize( m_OrientNum );
      for (unsigned int v = 0; v < m_OrientNum; v++) {
        tempDescriptor[z][v].resize( m_OrientNum );
        for (unsigned int w = 0; w < m_OrientNum; w++) {
          tempDescriptor[z][v][w].fill( 0 );
        }
      }
    }

    TIMER( timerPreparingVectors.Stop() );


    TIMER( timerAddingToBins.Start() );
    //for( std::vector< vnl_vector_fixed< double, 3 > >::size_type  i = 0; i < 1; ++i ) {
    for( std::vector< vnl_vector_fixed< double, 3 > >::size_type  i = 0; i < context_vectors.size(); ++i ) {

      // assuming that featureNeighbors and context_vectors have the same size and the order has not been changed
      unsigned int ind = featureNeighbors[i];

      this->AddToBin( context_vectors[i], keypoint_azimuth, keypoint_elevation, ind, tempDescriptor, keypointNormal );
    }


    // Rotate the values in the bin w.r.t. the keypoint orientation.
    for (unsigned int z = 0; z < m_RadiusNum+1; z++) {
      tempDescriptor[z].resize( m_OrientNum );
      for (unsigned int v = 0; v < m_OrientNum; v++) {
        tempDescriptor[z][v].resize( m_OrientNum );
        for (unsigned int w = 0; w < m_OrientNum; w++) {
          tempDescriptor[z][v][w] = rotation * tempDescriptor[z][v][w];
        }
      }
    }
    TIMER( timerAddingToBins.Stop() );




    // copy the descriptor bin values into the descriptor vector; normalize the bins
    //
    TIMER( timerCopyingDescriptor.Start() );
    vnl_vector< float >  descriptor;
    this->CopyDescriptor( tempDescriptor, descriptor );
    TIMER( timerCopyingDescriptor.Stop() );


    DescriptorAttribute  attribute( keypoint_direction, keypoint_bidirection, descriptor );
    m_KeypointsWithDescriptors[threadId].push_back( keypointPoint );
    m_KeypointsWithDescriptorsData[threadId].push_back( attribute );

  }

  //std::cout << "Not enough features to compute the descriptor for " << notEnoughPoints << " keypoints." << std::endl;
  //std::cout << "Cannot compute orientation for " << badOrientation << " keypoints." << std::endl;


#if USE_TIMER
  std::cout << "Building tree took " << timerBuildingTree.GetMeanTime() << " seconds.\n";
  std::cout << "Getting features took " << timerGettingFeatures.GetTotal() << " seconds.\n";
  std::cout << "Preparing vectors took " << timerPreparingVectors.GetTotal() << " seconds.\n";
  std::cout << "Adding to bins took " << timerAddingToBins.GetTotal() << " seconds.\n";
  std::cout << "Copying descriptor took " << timerCopyingDescriptor.GetTotal() << " seconds.\n";
  std::cout << "Computing all descriptors took " << timerComputingDescriptors.GetMeanTime() << " seconds.\n";
#endif

  //std::cout << "Computing ATBprepare took " << timerATBprepare.GetTotal() << " seconds.\n";
  //std::cout << "Computing ATBcompute took " << timerATBcompute.GetTotal() << " seconds.\n";

}


template <class TInputMesh1, class TInputMesh2, class TOutputMesh>
ITK_THREAD_RETURN_TYPE 
DescriptorMeshFilter< TInputMesh1, TInputMesh2, TOutputMesh>
::ThreaderCallback( void *arg )
{
  ThreadStruct *str;
  int total, threadId, threadCount;

  threadId = ((MultiThreader::ThreadInfoStruct *)(arg))->ThreadID;
  threadCount = ((MultiThreader::ThreadInfoStruct *)(arg))->NumberOfThreads;

  str = (ThreadStruct *)(((MultiThreader::ThreadInfoStruct *)(arg))->UserData);

  // execute the actual method with appropriate output region
  // first find out how many pieces extent can be split into.
  PointsContainer::iterator  startIterator;
  PointsContainer::iterator  endIterator;
  total = str->Filter->SplitRequestedList(threadId, threadCount,
                                          startIterator, endIterator);

  if (threadId < total)
    {
    str->Filter->ComputeDescriptorsThreaded( startIterator, endIterator, threadId );
    }
  // else
  //   {
  //   otherwise don't use this thread. Sometimes the threads dont
  //   break up very well and it is just as efficient to leave a 
  //   few threads idle.
  //   }
  
  return ITK_THREAD_RETURN_VALUE;
}



/**
 * Standard "PrintSelf" method
 */
template <class TInputMesh1, class TInputMesh2, class TOutputMesh>
void
DescriptorMeshFilter< TInputMesh1, TInputMesh2, TOutputMesh>
::PrintSelf(
  std::ostream& os,
  Indent indent) const
{
  Superclass::PrintSelf( os, indent );

}

} // end namespace itk





#endif
