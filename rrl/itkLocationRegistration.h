#ifndef __itkLocationRegistration_h
#define __itkLocationRegistration_h

#include <vul/vul_file.h>
#include <vul/vul_arg.h>

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <fstream>

#include "itkImageSeriesReader.h"
#include "itkImageSeriesWriter.h"
#include "itkNumericSeriesFileNames.h"
#include "itkRegularExpressionSeriesFileNames.h"

#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"

#include "itkRawImageIO.h"

#include "itkKdTree.h"
#include "itkKdTreeGenerator.h"
#include "itkListSample.h"

#include "itkImage.h"
#include "itkWarpImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkIdentityTransform.h"
#include "itkImageRegion.h"
#include "itkVectorLinearInterpolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"

#include "itkBSplineDeformableTransform.h"
#include "itkTransformFileReader.h"
#include "itkTransformFileWriter.h"
#include "itkTransformFactory.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkCenteredAffineTransform.h"
#include "itkResampleImageFilter.h"
#include "itkRegionOfInterestImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkCheckerBoardImageFilter.h"
#include "itkTileImageFilter.h"
#include "itkMaximumProjectionImageFilter.h"
#include "itkImageSeriesWriter.h"
#include "itkCastImageFilter.h"
#include "itkVTKImageIO.h"
#include "itkTimeProbe.h"

#include "itkDanielssonDistanceMapImageFilter.h"

#include "itkLightProcessObject.h"

#include "itkRescaleIntensityImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkMaximumProjectionImageFilter.h"
#include "itkResampleImageFilter.h"
#include "itkTileImageFilter.h"
#include "itkPNGImageIO.h"

#include <cdcl/cdcl_keypoint.h>
#include <cdcl/cdcl_estimation_ICP.h>
//#include <cdcl/cdcl_estimation_ICP_with_shape.h>
#include <cdcl/cdcl_estimation_ICP_matching_all.h>
//#include <cdcl/cdcl_estimation_symmetric_ICP.h>
#include <cdcl/cdcl_estimation_symmetric_ICP_matching_all.h>
#include <rrl/rrl_estimation_symmetric_ICP_matching_all.h>
#include <cdcl/cdcl_estimation_transfer.h>
#include <cdcl/cdcl_trans_affine.h>
#include <cdcl/cdcl_trans_rigid3d.h>
#include <cdcl/cdcl_trans.h>
#include <cdcl/cdcl_feature_with_shape.h>
#include <cdcl/cdcl_feature_ICP.h>
#include <cdcl/cdcl_match.h>
#include <cdcl/cdcl_utils.h>
#include <cdcl/cdcl_utils_VTK.h>

#include <cdcl/io/vtkPolyDataToFeaturesWithShapeFilter.h>
#include <cdcl/io/vtkPolyDataToFeaturesICPFilter.h>

#include "vtkPolyData.h"
#include "vtkPolyDataWriter.h"
#include "vtkSmartPointer.h"

#include <cdcl/io/itkImageSlicesWithGeometryFilter.h>

#include <vnl/vnl_inverse.h>
#include <vnl/vnl_cross.h>
#include <vnl/vnl_trace.h>
#include <vnl/vnl_vector_fixed.txx>


#include <libsvm/svm_interface.h>


//#include <rrl/rrl_cdcl_registration.h>

//:
// \file
// \brief  Register volume location using rrl_cdcl_registration object.
// \author Michal Sofka
// \date   Oct 2007

#define SYMMETRIC_MATCHING 1
#define WRITE_PANELS 0
#define LAST_ITER 1
#define DECISION_DEBUGGING 0


namespace itk {


class ITK_EXPORT LocationRegistration : public LightProcessObject
{
public:
  LocationRegistration();


  /** Standard class typedefs. */
  typedef LocationRegistration  Self;
  //typedef MeshToMeshFilter<TInputMesh,TOutputMesh> Superclass;

  typedef SmartPointer<Self>  Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);  

  const static unsigned int Dimension = 3;
  typedef signed short PixelType;
  //typedef float PixelType;

  typedef itk::RGBPixel< unsigned char > RGBPixelType;

  typedef itk::Image<RGBPixelType, 2>                           OutputImageType;

  typedef itk::Image< PixelType, Dimension >  ImageType;

  typedef itk::Image< unsigned int, Dimension >   SegmentedImageType;

  typedef itk::Image< unsigned int, 3 >  VoronoiMapImageType;

  #if CDC
  typedef cdcl_feature_with_shape< 3 >  feature_type;
  typdef vtkPolyDataToFeaturesWithShapeFilter  PolyDataToFeaturesFilterType;
  typdef vtkFeatureWithShapeAttributeSet  FeatureAttributeSet;
  #else
  typedef cdcl_feature_ICP< 3 >  feature_type;
  typedef vtkPolyDataToFeaturesICPFilter  PolyDataToFeaturesFilterType;
  typedef vtkFeatureICPAttributeSet  FeatureAttributeSet;
  #endif
  typedef feature_type::sptr  feature_sptr_type;

  itkSetMacro(FixedImageDir, std::string);
  itkSetMacro(MovingImageDir, std::string);
  itkSetMacro(TransformFile, std::string);
  itkSetMacro(Prefix, std::string);
  itkSetMacro(SegmentationSuffix, std::string);
  
  int Run();

  void GenerateSlices( vtkSmartPointer< vtkPolyData > const &  polyDataMovingMapped,
                           vtkSmartPointer< vtkPolyData > const &  polyDataFixed,
                           ImageType::Pointer const &  m_MovingMappedImage,
                           ImageType::Pointer const & m_FixedImage,
                           itk::ImageRegion< 3 > const &  m_ROI,
                           OutputImageType::Pointer &  m_SlicesImage );

  // Read moving and fixed images. These are used for visualization.
  void ReadImages();

  // Read oversegmented moving and fixed images.
  void ReadSegmentedImages();

  // Read features used for registration.
  void ReadFeatures();

  // Read voronoi maps for rapid matching, these are volumes holding indices to features vector.
  int ReadVoronoiMaps();

  // Read moving keypoints and descriptors and build a kd-tree for each.
  void ReadMovingKeypointDescriptors();

  void ReadFixedKeypointDescriptors();

  int SetupGroundTruthTransform();

  static const unsigned int dim = 3;

  //static const unsigned int dof = 6;
  //typedef cdcl_trans_rigid3d  trans_type;
  static const unsigned int dof = 12;
  typedef cdcl_trans_affine<dim>  trans_type;


  // Get CDC transform from ITK transform.
  cdcl_trans< dim, dof >::sptr  GetCDCTransform();

  // Update ITK transform given CDC transform.
  void UpdateITKTransform();

  void SetupFinalTransform();

  void ReadFinalTransform( std::string  TransformFile );

  // Using the oversegmented images, get features that will be used in estimation.
  void FeaturesInRegions( vcl_vector< feature_sptr_type > &  moving_inside,
                          vcl_vector< feature_sptr_type > &  fixed_inside,
                          vnl_vector_fixed< double, dim >  const &  moving_x0,
                          vnl_vector_fixed< double, dim >  const &  moving_x1,
                          vnl_vector_fixed< double, dim >  const &  fixed_x0,
                          vnl_vector_fixed< double, dim >  const &  fixed_x1 );

  // Get features inside the regions of interest, without using the oversegmented images. 
  void FeaturesInROIs( vcl_vector< feature_sptr_type > &  moving_inside,
                       vcl_vector< feature_sptr_type > &  fixed_inside,
                       vnl_vector_fixed< double, dim >  const &  moving_x0,
                       vnl_vector_fixed< double, dim >  const &  moving_x1,
                       vnl_vector_fixed< double, dim >  const &  fixed_x0,
                       vnl_vector_fixed< double, dim >  const &  fixed_x1 );

  // Create estimation object and initialize.
  bool SetupEstimation();

  void GetFeaturesInMovingROI();

  void WarpVolume();

  void WriteResults( cdcl_trans<dim,dof>::sptr const &  initial_transform,
                     double const &  initialRMS,
                     bool converged,
                     unsigned int k, // iteration
                     unsigned int candidateNum,
                     std::vector< double > const &  measurements,
                     std::string  prefix = "" );

  void ComputeDecisionMeasurements( std::vector< double > &  measurements );

private:

  const static unsigned int ImageDimension = 3;
	const static unsigned int SpaceDimension = ImageDimension;
	const static unsigned int SplineOrder = 3;
	typedef double CoordinateRepType;

  // BSpline
	typedef itk::BSplineDeformableTransform<
													CoordinateRepType,
													SpaceDimension,
													SplineOrder >     BSplineTransformType;


  // Deformation field
  typedef itk::Vector< float, ImageDimension >           VectorPixelType;
  typedef itk::Image<  VectorPixelType, ImageDimension > DeformationFieldType;


  void ComputeTransferErrorCovariance( itk::Image< signed short, dim >::Pointer const &  Image,
                                     itk::ImageRegion< dim > const &  m_ROImoving,
                                     vbl_smart_ptr< cdcl_trans<dim, dof> > const & transform,
                                     vnl_matrix_fixed< double, dim, dim > &  covarianceJ,
                                     double &  maxTrace,
                                     double &  maxEval );

  bool validTransform( vbl_smart_ptr< cdcl_trans< dim, dof > > const &  trans );


  void MeanAndScaleOfAffineApproximationError( DeformationFieldType::Pointer const &  deformationField,
                                             itk::ImageRegion< 3 > const &  m_ROI,
                                             vnl_vector_fixed< double, 3 > const &  fixed_keypoint_location,
                                             vnl_vector_fixed< double, 3 > const &  moving_keypoint_location,
                                             vnl_matrix_fixed< double, 3, 3 > const &  A,
                                             vnl_vector_fixed< double, 3 > const &  t,
                                             double &  medianDeformationComponentLSError,
                                             double &  scaleDeformationComponentLSError );

  void LocationRegistration::FitAffineTransformToDeformationField( DeformationFieldType::Pointer const &  deformationField,
                                           itk::ImageRegion< 3 > const &  m_ROI,
                                           vnl_vector_fixed< double, 3 > const &  fixed_keypoint_location,
                                           vnl_vector_fixed< double, 3 > const &  moving_keypoint_location,
                                           vnl_matrix_fixed< double, 3, 3 > &  A,
                                           vnl_vector_fixed< double, 3 > &  t );


  // Kd-tree for moving keypoint locations
  typedef itk::Vector< float, 3 >                                 KeypointVectorType;
  typedef itk::Statistics::ListSample< KeypointVectorType >       KeypointSampleType;
  typedef itk::Statistics::KdTreeGenerator< KeypointSampleType >  KeypointTreeGeneratorType;
  typedef KeypointTreeGeneratorType::KdTreeType                   KeypointTreeType;

  // Kd-tree for descriptor vectors
  //static const unsigned int descriptorSize = 771;
  //typedef itk::Vector< float, descriptorSize >                      DescriptorVectorType;
  typedef itk::Array< float >                      DescriptorVectorType;
  typedef itk::Statistics::ListSample< DescriptorVectorType >       DescriptorSampleType;
  typedef itk::Statistics::KdTreeGenerator< DescriptorSampleType >  DescriptorTreeGeneratorType;
  typedef DescriptorTreeGeneratorType::KdTreeType                   DescriptorTreeType;


  DescriptorTreeType::Pointer  kdTreeMovingDescriptors;
  KeypointTreeType::Pointer  kdTreeMovingKeypoints;


  BSplineTransformType::Pointer  bsplineTransform;


  DeformationFieldType::Pointer  deformationField;

  enum deformationTypeEnum { NONE, BSPLINE, FIELD };
  deformationTypeEnum  deformationType;



  ImageType::Pointer  m_FixedImage;
  ImageType::Pointer  m_MovingImage;

  SegmentedImageType::Pointer  m_FixedImageSegmented;
  SegmentedImageType::Pointer  m_MovingImageSegmented;

  // Mapped moving image (in the m_ROImoving).
  ImageType::Pointer  m_MovingMappedImage;

  vcl_vector< vcl_vector< feature_sptr_type > >  m_FeaturesFixedAll;
  vcl_vector< vcl_vector< feature_sptr_type > >  m_FeaturesMovingAll;
  vcl_vector< vcl_vector< feature_sptr_type > >  m_FeaturesMovingROI;

  vcl_vector< cdcl_keypoint< 3 >::sptr >               fixedKeypoints;
  vcl_vector< vnl_vector< float > >                    fixedDescriptors;

#if SYMMETRIC_MATCHING
  //typedef cdcl_estimation_symmetric_ICP< dim, dof >  estimation_symmetric_type;
  //typedef cdcl_estimation_symmetric_ICP_matching_all< dim, dof >  estimation_symmetric_type;
  typedef rrl_estimation_symmetric_ICP_matching_all< dim, dof >  estimation_symmetric_type;
  estimation_symmetric_type::sptr  estimation;
#else
  //typedef cdcl_estimation_ICP_matching_all< dim, dof >  estimation_type;
  //typedef cdcl_estimation_ICP_with_shape< dim, dof >  estimation_type;

  #define CDC 0
  //typedef feature_type  feature_type;
  typedef cdcl_estimation<dim, dof>  estimation_type;
  estimation_type::sptr  estimation;
#endif

  typedef itk::CenteredAffineTransform< double, 3 >  TransformType;
  TransformType::Pointer  m_FinalTransform;

  cdcl_trans<dim,dof>::sptr  transform;

  itk::Point< double, 3 >  queryPointMapped;
  itk::Point< double, 3 >  queryPoint;
  cdcl_keypoint<3>::sptr  moving_keypoint;
  cdcl_keypoint<3>::sptr  fixed_keypoint;

  // Current region of interest (points and image box) in the fixed coordinate system.
  itk::ImageRegion<3>  m_ROI;

  // Current region of interest m_ROI inverse mapped to the moving coordinate system.
  itk::ImageRegion<3>  m_ROImoving;

  vcl_vector< feature_sptr_type >  moving_inside;
  vcl_vector< feature_sptr_type >  fixed_inside;

  VoronoiMapImageType::Pointer  m_FixedVoronoiMap;
  VoronoiMapImageType::Pointer  m_MovingVoronoiMap;

  unsigned int m_ROISize[3];

  std::string  m_FixedImageDir;
  std::string  m_MovingImageDir;
  std::string  m_TransformFile;
  std::string  m_Prefix;
  std::string  m_SegmentationSuffix;

  std::string  locRegDir;
  std::string  matchesDirPath;


};

}


#endif
