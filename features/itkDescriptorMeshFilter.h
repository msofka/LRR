/*=========================================================================


author: Michal Sofka, Rensselaer Polytechnic Institute (RPI)
  date: 10/10/2007
  
acknowledgements: Inspired by LocalMaximumImageFilter ( Bryn A. Lloyd, Simon K. Warfield ),
by ImageToMeshFilter (ITK), and by RPI's feature extraction algorithms.

=========================================================================*/

#ifndef __itkDescriptorMeshFilter_h
#define __itkDescriptorMeshFilter_h

#include "itkMeshSource.h"


#include "itkGradientRecursiveGaussianImageFilter.h"
#include "itkSymmetricSecondRankTensor.h"

#include "itkKdTreeForThreading.h"
#include "itkKdTreeForThreadingGenerator.h"
#include "itkListSample.h"

#include "itkTimeProbe.h"
namespace itk
{
/* \class DescriptorMeshFilter
 \brief  Filter for computing Shape-Context-like keypoint descriptors.
         The first mesh are the keypoints and the second mesh features.
         Features are used to compute the descriptor for each keypoint.

 \sa Image

 \ingroup MeshFilters
 */



// Attributes (point data) stored with each point location (point) in the mesh
class DescriptorAttribute {
public:
  typedef DescriptorAttribute  Self;
  typedef SmartPointer< Self >  Pointer;

  itk::Vector< float, 3 >  m_Direction;
  itk::Vector< float, 3 >  m_Bidirection;
  vnl_vector< float >  m_Descriptor;

  DescriptorAttribute()
    : m_Direction( 0.0f ),
      m_Bidirection( 0.0f )
  { }

  DescriptorAttribute( itk::Vector< float, 3 > const &  direction,
                       itk::Vector< float, 3 > const &  bidirection,
                       vnl_vector< float > const &  descriptor )
    : m_Direction( direction ), m_Bidirection( bidirection ), m_Descriptor( descriptor )
  { }
};


template <class TInputMesh1, class TInputMesh2, class TOutputMesh>
class ITK_EXPORT DescriptorMeshFilter :
    public MeshSource<TOutputMesh>
{
public:

  /** Run-time type information (and related methods). */
  itkTypeMacro(DescriptorMeshFilter, ImageToMeshFilter);


  /** Some typedefs associated with the input meshes. */
  typedef TInputMesh1                               KeypointMeshType;
  typedef TInputMesh2                               FeatureMeshType;

  /** Some typedefs associated with the output mesh. */
  typedef TOutputMesh OutputMeshType;
  typedef typename OutputMeshType::PointType        PointType;
  typedef typename OutputMeshType::Pointer          OutputMeshPointer;
  typedef typename OutputMeshType::ConstPointer     OutputMeshConstPointer;
  typedef typename OutputMeshType::PointsContainer  PointsContainer;
  typedef typename OutputMeshType::PointIdentifier  PointIdentifier;
  typedef typename PointsContainer::Pointer         PointsContainerPointer;
  typedef typename PointsContainer::Iterator        PointsContainerIterator;
  typedef typename OutputMeshType::PointDataContainer PointDataContainer;
  typedef typename PointDataContainer::Pointer      PointDataContainerPointer;
  typedef typename PointDataContainer::Iterator     PointDataContainerIterator;

  typedef typename KeypointMeshType::PointsContainer::iterator  KeypointIteratorType;
  typedef itk::Vector< float, 3 >                                 LocationVectorType;

  typedef itk::Statistics::ListSample< LocationVectorType >       FeatureSampleType;
  typedef itk::Statistics::KdTreeForThreadingGenerator< FeatureSampleType >   FeatureTreeGeneratorType;
  typedef FeatureTreeGeneratorType::KdTreeType                    FeatureTreeType;

  /** Standard class typedefs. */
  typedef DescriptorMeshFilter Self;
  typedef MeshSource< OutputMeshType> Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** accept the input mesh */
  void SetInputFeatures( const TInputMesh2 *input );

  /** accept the input mesh */
  void SetInputKeypoints( const TInputMesh1 *input );

  /** Some typedefs associated with the output mesh. */
  void GenerateOutputInformation(void);

protected:
  DescriptorMeshFilter();
  virtual ~DescriptorMeshFilter() {}
  void PrintSelf(std::ostream& os, Indent indent) const;

  void GenerateData();

private:
  DescriptorMeshFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented


  // Add a value to the descriptor bin.
  __forceinline void AddToBin( itk::Vector< float, 3 > const &  norm_vector,
                 double const &  keypoint_azimuth,
                 double const &  keypoint_elevation,
                 unsigned int ind,
                 std::vector< std::vector< std::vector< vnl_vector_fixed< float, 3 > > > > &  tempDescriptor,
                 itk::Vector< float, 3 > const &  keypointNormal );

  // Copy the descriptor bin values into the descriptor vector; normalize the bins.
  void CopyDescriptor( std::vector< std::vector< std::vector< vnl_vector_fixed< float, 3 > > > > &  tempDescriptor,
                       vnl_vector< float > &  descriptor );

  struct ThreadStruct
    {
    Pointer Filter;
    };

  // Split list for multithreading.
  int SplitRequestedList( int i, int num,
                          KeypointIteratorType &  startIterator,
                          KeypointIteratorType &  endIterator );

  void ComputeDescriptorsThreaded( KeypointIteratorType &  startIterator,
                                   KeypointIteratorType &  endIterator,
                                   int threadId );

  static ITK_THREAD_RETURN_TYPE ThreaderCallback( void *arg );

  std::list< typename KeypointMeshType::PointType > *  m_KeypointsWithDescriptors;
  std::list< DescriptorAttribute > *  m_KeypointsWithDescriptorsData;

  // Radius of the descriptor region.
  float m_Radius;

  // number of bins in each coordinate
  // flipping the angles w.r.t azimuth and elevation -> only half of bins needed
  unsigned int m_RadiusNum;
  unsigned int m_OrientNum;

//itk::TimeProbe  timerATBprepare;
//itk::TimeProbe  timerATBcompute;

  // bin sizes
  double m_LogRadiusBinSize;  // on the logarithmic scale
  double m_AzimuthBinSize;
  double m_ElevationBinSize;

  typename FeatureTreeType::Pointer  m_kdTreeFeatures;

  // 2*pi
  static const double m_2pi;

};


} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkDescriptorMeshFilter.txx"
#endif

#endif
