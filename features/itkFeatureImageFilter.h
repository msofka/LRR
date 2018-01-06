/*=========================================================================


author: Michal Sofka, Rensselaer Polytechnic Institute (RPI)
  date: 10/10/2007
  
acknowledgements: Inspired by LocalMaximumImageFilter ( Bryn A. Lloyd, Simon K. Warfield ),
by ImageToMeshFilter (ITK), and by RPI's feature extraction algorithms.

=========================================================================*/

#ifndef __itkFeatureImageFilter_h
#define __itkFeatureImageFilter_h

#include "itkImageToMeshFilter.h"
#include "itkImage.h"

#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkImageRegionIteratorWithIndex.h"

#include "itkLinearInterpolateImageFunction.h"

#include "itkGradientRecursiveGaussianImageFilter.h"
#include "itkSymmetricSecondRankTensor.h"


namespace itk
{
/* \class FeatureImageFilter
 \brief  Extracting features as a measure of the outer products of the local gradients.
         This is followed by non-maximum suppression, local contrast filtering,
         spatial filtering and (optional) covariance computation.

 \sa Image

 \ingroup IntensityImageFilters
 */


#define COMPUTE_COVARIANCES 0


// Attributes (point data) stored with each point location (point) in the mesh
class PointAttribute {
public:
  typedef PointAttribute  Self;
  typedef SmartPointer< Self >  Pointer;
  typedef SymmetricSecondRankTensor< float, 3 >  SymmetricMatrixType;
  enum FeatureShape { Corner, Tube, Sheet };

  float m_Strength;
  #if COMPUTE_COVARIANCES       
  SymmetricMatrixType  m_Covariance;
  #endif
  SymmetricMatrixType  m_ErrorProjector;
  FeatureShape         m_Shape;
  // directions which define the feature (3 directions for a corner, 2 for a tube and 1 for a sheet)
  // the ordering is: tangent, binormal, normal (this is given by the eigen value ordering
  // if this ordering is changed, then reading / writing / converting filters need to be changed
  std::vector< itk::Vector< float, 3 > >  m_Directions;

  PointAttribute()
    : m_Strength( -1.0f )
  { }

  #if COMPUTE_COVARIANCES 
  PointAttribute( float const &  Strength, SymmetricMatrixType const &  Covariance, SymmetricMatrixType const &  ErrorProjector, FeatureShape Shape, std::vector< itk::Vector< float, 3 > > const &  Directions )
    : m_Strength( Strength ), m_Covariance( Covariance ), m_ErrorProjector( ErrorProjector ), m_Shape( Shape ), m_Directions( Directions )
  { }
  #else
  PointAttribute( float const &  Strength, SymmetricMatrixType const &  ErrorProjector, FeatureShape Shape, std::vector< itk::Vector< float, 3 > > const &  Directions )
    : m_Strength( Strength ), m_ErrorProjector( ErrorProjector ), m_Shape( Shape ), m_Directions( Directions )
  { }
  #endif
};


template <class TInputImage, class TOutputMesh>
class ITK_EXPORT FeatureImageFilter :
    public ImageToMeshFilter<TInputImage,TOutputMesh>
{
public:

  /** Run-time type information (and related methods). */
  itkTypeMacro(FeatureImageFilter, ImageToMeshFilter);


  /** Extract dimension from input and output image. */
  itkStaticConstMacro(InputImageDimension, unsigned int,
                      TInputImage::ImageDimension);

  /** Some typedefs associated with the input images. */
  typedef TInputImage                               InputImageType;
  typedef typename InputImageType::Pointer          InputImagePointer;
  typedef typename InputImageType::ConstPointer     InputImageConstPointer;
  typedef typename InputImageType::RegionType       InputImageRegionType;
  typedef ImageRegionConstIteratorWithIndex<InputImageType>
                                                    InputImageIterator;

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

  typedef Image< float, InputImageDimension >  InternalImageType;
  // Gradient vector, type float as elements rather than default double.
  typedef typename CovariantVector< typename NumericTraits< typename InternalImageType::PixelType >::FloatType, InputImageDimension >  GradientVectorType;
  typedef typename Image< GradientVectorType, InputImageDimension >  GradientImageType;

  typedef itk::GradientRecursiveGaussianImageFilter< InternalImageType, GradientImageType >  GradientImageFilterType;
  typedef SymmetricSecondRankTensor< float, 3 >  SymmetricMatrixType;

  /** Image typedef support. */
  typedef typename InputImageType::PixelType InputPixelType;
  typedef typename OutputMeshType::PixelType OutputPixelType;



  typedef typename InputImageType::SizeType InputSizeType;
  typedef typename InputImageType::IndexType InputIndexType;

  typedef itk::LinearInterpolateImageFunction< InternalImageType >  InterpolatorType;

  /** Standard class typedefs. */
  typedef FeatureImageFilter Self;
  typedef ImageToMeshFilter< InputImageType, OutputMeshType> Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Set the radius of the neighborhood used to filter the scores. */
  itkSetMacro(FilteringRadius, InputSizeType);
  /** Get the radius of the neighborhood used to filter the scores */
  itkGetConstReferenceMacro(FilteringRadius, InputSizeType);

  itkSetMacro(FilterAcrossScales, bool);
  itkSetMacro(OnlyCorners, bool);

  /** accept the input image */
  void SetInput( const InputImageType * inputImage );

  /** Some typedefs associated with the output mesh. */
  void GenerateOutputInformation(void);

protected:
  FeatureImageFilter();
  virtual ~FeatureImageFilter() {}
  void PrintSelf(std::ostream& os, Indent indent) const;

  void GenerateData();

private:
  FeatureImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  // Storage of intermediate feature computation results.
  class InternalRecord {
  public:
    enum FeatureType { Corner, Tube, Sheet };
    
    float m_Strength;
    FeatureType  m_PointType;
    // after subvoxel accuracy, this field will hold the refined location
    itk::Point< float, InputImageDimension >  m_Location;
    typename InputImageType::IndexType  m_LocationIndex;
    // directions which define the feature (3 directions for a corner, 2 for a tube and 1 for a sheet)
    // the ordering is: tangent, binormal, normal (this is given by the eigen value ordering
    // if this ordering is changed, then reading / writing / converting filters need to be changed
    std::vector< itk::Vector< float, InputImageDimension > >  m_Directions;

    #if COMPUTE_COVARIANCES
    // Outer product of the gradients.
    SymmetricMatrixType  m_OuterProduct;
    #endif
    
    InternalRecord()
    : m_Strength( -1.0f )
    {   }

  };

  typedef itk::Image< unsigned char, InputImageDimension >  FilteredPointsIndicatorImageType;

  typedef typename itk::Point< float, InputImageDimension >  PhysicalPointType;

  // Split requested region for multithreading (from ImageSource).
  int SplitRequestedRegion( int i, int num, typename InternalImageType::RegionType &  splitRegion, bool shrink = false );

  // Compute score image by evaluating a response measure at each pixel as a function of the outer product of the gradients.
  void ComputeScoreImageThreaded( const typename InternalImageType::RegionType &  outputRegion,
                                  int  threadId );

  // Filter the features based on contrast not to give rise large amounts of features in textured regions.
  // This is done by computing median response in a region and thresholding at median + std. dev.
  void  LocalContrastFiltering( typename InternalImageType::Pointer &  scoreImage );

  // Non maximum suppression of responses.
  void NonMaximumSuppressionThreaded( typename InternalImageType::RegionType const &  outputRegion,
                                      int  threadId );

  // Split list for multithreading.
  int SplitRequestedList( int i, int num,
                          typename std::list< InternalRecord >::iterator &  startIterator,
                          typename std::list< InternalRecord >::iterator &  endIterator );

  // Run subvoxel localization algorithm for all detected points.
  void SubvoxelAccuracyThreaded( typename std::list< InternalRecord >::iterator &  startIterator,
                                 typename std::list< InternalRecord >::iterator &  endIterator,
                                 int threadId );

  // Subvoxel localization for a feature using all points in the neighborhood and fitting a polynomial.
  bool SubvoxelLocalizationLS( InternalRecord &  result,
                               typename InterpolatorType::Pointer const &  interpolator,
                               itk::Matrix< float, 10, 10 > const &  LS );

  // Compute improved directions using average of (interpolated) local derivatives.
  void ComputeDirections( std::list< InternalRecord > &  unfilteredPoints );

  // Covariance computation for the filtered points.
  void ComputeCovariances( std::list< InternalRecord > const &  unfilteredPoints );

  static ITK_THREAD_RETURN_TYPE ThreaderCallbackScore( void *arg );
  static ITK_THREAD_RETURN_TYPE ThreaderCallbackNonMax( void *arg );
  static ITK_THREAD_RETURN_TYPE ThreaderCallbackSubvoxel( void *arg );

  struct ThreadStruct
    {
    Pointer Filter;
    };

  typename InternalImageType::Pointer  m_ScoreImage;
  typename GradientImageFilterType::OutputImagePointer  m_GradientImage;
  typename FilteredPointsIndicatorImageType::Pointer  m_FilteredPointsImage;
  std::list< InternalRecord > *  m_NonMaxedPoints;

  std::list< InternalRecord >  m_UnfilteredPoints;

  // Radius of the spatial filtering.
  InputSizeType  m_FilteringRadius;

  double m_Sigma;

  // Indicates whether to filter scores accross scales.
  bool m_FilterAcrossScales;

  // Detect only corners (not tubes or sheets).
  // Because of filtering, different set of corners is obtained than when all feature types are allowed.
  bool m_OnlyCorners;

};


} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkFeatureImageFilter.txx"
#endif

#endif
