/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkDemonsImageToImageMetric.h,v $
  Language:  C++
  Date:      $Date: $
  Version:   $Revision: $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkImageToImageObjectMetric_h
#define __itkImageToImageObjectMetric_h

#include "itkCovariantVector.h"
#include "itkCentralDifferenceImageFunction.h"
#include "itkGradientRecursiveGaussianImageFilter.h"
#include "itkObjectToObjectMetric.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkBSplineInterpolateImageFunction.h"
#include "itkSpatialObject.h"
#include "itkImageToData.h"

namespace itk
{
/** \class ImageToImageObjectMetric
 *
 * Templated over the fixed and moving image types, as well as an optional
 * VirtualImage type to define the virtual domain. The VirtualImage type
 * defaults to TFixedImage.
 *
 * Image gradient calculation is performed in one of three ways:
 * 1) The BSplineInterpolator is used if it has been set by user as the
 * interpolator.
 * 2) otherwise, GradientRecursiveGaussianImageFilter is used by default to
 * precompute the gradient for each image. This requires memory to hold the
 * fully computed gradient image.
 * 3) lastly, CentralDifferenceImageFunction is used if m_PrecomputeGradient
 * option has been set to false.
 *
 * Derived classes must implement pure virtual methods declared in
 * ObjectToObjectMetric.
 */
template<class TFixedImage,class TMovingImage,class TVirtualImage = TFixedImage>
class ITK_EXPORT ImageToImageObjectMetric :
public ObjectToObjectMetric<TFixedImage, TMovingImage>
{
public:

  /** Standard class typedefs. */
  typedef ImageToImageObjectMetric                          Self;
  typedef ObjectToObjectMetric<TFixedImage, TMovingImage>   Superclass;
  typedef SmartPointer<Self>                                Pointer;
  typedef SmartPointer<const Self>                          ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro(Self, Superclass);

  /** Type used internally for computations */
  typedef double  InternalComputationValueType;

  /** Type used for representing parameter values  */
  typedef typename Superclass::CoordinateRepresentationType
                                                  CoordinateRepresentationType;

  /**  Type of the parameters. */
  typedef typename Superclass::ParametersType       ParametersType;
  typedef typename Superclass::ParametersValueType  ParametersValueType;

  /** Derivative source type */
  typedef typename Superclass::DerivativeSourceType DerivativeSourceType;

  /** Image-accessor typedefs */
  typedef TFixedImage                             FixedImageType;
  typedef typename FixedImageType::PixelType      FixedImagePixelType;
  typedef typename FixedImageType::Pointer        FixedImagePointer;
  typedef typename FixedImageType::ConstPointer   FixedImageConstPointer;
  typedef typename FixedImageType::PointType      FixedImagePointType;
  typedef typename FixedImageType::IndexType      FixedImageIndexType;
  typedef TMovingImage                            MovingImageType;
  typedef typename MovingImageType::PixelType     MovingImagePixelType;
  typedef typename MovingImageType::Pointer       MovingImagePointer;
  typedef typename MovingImageType::ConstPointer  MovingImageConstPointer;
  typedef typename MovingImageType::PointType     MovingImagePointType;
  typedef typename MovingImageType::RegionType    MovingImageRegionType;
  typedef typename MovingImageType::IndexType     MovingImageIndexType;
  /** Types for the virtual domain */
  typedef TVirtualImage                             VirtualImageType;
  typedef typename VirtualImageType::PixelType      VirtualImagePixelType;
  typedef typename VirtualImageType::Pointer        VirtualImagePointer;
  typedef typename VirtualImageType::RegionType     VirtualRegionType;
  typedef typename VirtualRegionType::SizeType      VirtualSizeType;
  typedef typename VirtualImageType::SpacingType    VirtualSpacingType;
  typedef typename VirtualImageType::PointType      VirtualOriginType;
  typedef typename VirtualImageType::PointType      VirtualPointType;
  typedef typename VirtualImageType::DirectionType  VirtualDirectionType;
  typedef typename VirtualImageType::SizeType       VirtualRadiusType;
  typedef typename VirtualImageType::IndexType      VirtualIndexType;

  /* Image dimension accessors */
  itkStaticConstMacro(FixedImageDimension, unsigned int,
      ::itk::GetImageDimension<FixedImageType>::ImageDimension);
  itkStaticConstMacro(MovingImageDimension, unsigned int,
      ::itk::GetImageDimension<MovingImageType>::ImageDimension);
  itkStaticConstMacro(VirtualImageDimension, unsigned int,
      ::itk::GetImageDimension<VirtualImageType>::ImageDimension);

  /**  Type of the Transform Base classes */
  typedef Transform<CoordinateRepresentationType,
    itkGetStaticConstMacro( MovingImageDimension ),
    itkGetStaticConstMacro( VirtualImageDimension )>  MovingTransformType;

  typedef Transform<CoordinateRepresentationType,
    itkGetStaticConstMacro( FixedImageDimension ),
    itkGetStaticConstMacro( VirtualImageDimension )> FixedTransformType;

  typedef typename FixedTransformType::Pointer         FixedTransformPointer;
  typedef typename FixedTransformType::InputPointType  FixedInputPointType;
  typedef typename FixedTransformType::OutputPointType FixedOutputPointType;
  typedef typename FixedTransformType::ParametersType
                                                FixedTransformParametersType;
  typedef typename FixedTransformType::JacobianType
                                                FixedTransformJacobianType;

  typedef typename MovingTransformType::Pointer         MovingTransformPointer;
  typedef typename MovingTransformType::InputPointType  MovingInputPointType;
  typedef typename MovingTransformType::OutputPointType MovingOutputPointType;
  typedef typename MovingTransformType::ParametersType
                                                MovingTransformParametersType;
  typedef typename MovingTransformType::JacobianType
                                                MovingTransformJacobianType;

  /**  Type for the mask of the fixed image. Only pixels that are "inside"
       this mask will be considered for the computation of the metric */
  typedef SpatialObject< itkGetStaticConstMacro(FixedImageDimension) >
                                                       FixedImageMaskType;
  typedef typename FixedImageMaskType::Pointer         FixedImageMaskPointer;
  typedef typename FixedImageMaskType::ConstPointer
                                                  FixedImageMaskConstPointer;

  /**  Type for the mask of the moving image. Only pixels that are "inside"
       this mask will be considered for the computation of the metric */
  typedef SpatialObject< itkGetStaticConstMacro(MovingImageDimension) >
                                                        MovingImageMaskType;
  typedef typename MovingImageMaskType::Pointer         MovingImageMaskPointer;
  typedef typename MovingImageMaskType::ConstPointer
                                                   MovingImageMaskConstPointer;

  /**  Type of the Interpolator Base class */
  typedef InterpolateImageFunction< FixedImageType,
                                    CoordinateRepresentationType >
                                                      FixedInterpolatorType;
  typedef InterpolateImageFunction< MovingImageType,
                                    CoordinateRepresentationType >
                                                      MovingInterpolatorType;
  typedef typename FixedInterpolatorType::Pointer     FixedInterpolatorPointer;
  typedef typename MovingInterpolatorType::Pointer    MovingInterpolatorPointer;

  /** Type of the default linear interpolators. */
  typedef LinearInterpolateImageFunction< FixedImageType,
                                          CoordinateRepresentationType >
                                                  FixedLinearInterpolatorType;
  typedef LinearInterpolateImageFunction< MovingImageType,
                                          CoordinateRepresentationType >
                                                  MovingLinearInterpolatorType;

  /**
   * If a BSplineInterpolationFunction is used, this class obtain
   * image derivatives from the BSpline interpolator. Otherwise,
   * image derivatives are computed using central differencing
   * or Gaussian gradient fileter. See main documentation.
   */
  typedef BSplineInterpolateImageFunction< MovingImageType,
                                           CoordinateRepresentationType >
                                                FixedBSplineInterpolatorType;
  typedef BSplineInterpolateImageFunction< MovingImageType,
                                           CoordinateRepresentationType >
                                                MovingBSplineInterpolatorType;

  /** Image derivatives types */
  typedef   CovariantVector< CoordinateRepresentationType,
                             itkGetStaticConstMacro(FixedImageDimension) >
                                                      FixedImageDerivativesType;
  typedef   CovariantVector< CoordinateRepresentationType,
                             itkGetStaticConstMacro(MovingImageDimension) >
                                                      MovingImageDerivativesType;

  /** Gaussian filter types to compute the gradient of the images.
   * This is used by default to compute image gradients. See comments
   * in main documentation */
  typedef typename NumericTraits< FixedImagePixelType >::RealType
                                                    FixedRealType;
  typedef CovariantVector< FixedRealType,
                           itkGetStaticConstMacro(FixedImageDimension) >
                                                    FixedGradientPixelType;
  typedef Image< FixedGradientPixelType,
                 itkGetStaticConstMacro(FixedImageDimension) >
                                                    FixedGradientImageType;
  typedef GradientRecursiveGaussianImageFilter< FixedImageType,
                                                FixedGradientImageType >
                                                 FixedGradientImageFilterType;

  typedef typename NumericTraits< MovingImagePixelType >::RealType
                                                    MovingRealType;
  typedef CovariantVector< MovingRealType,
                           itkGetStaticConstMacro(MovingImageDimension) >
                                                    MovingGradientPixelType;
  typedef Image< MovingGradientPixelType,
                 itkGetStaticConstMacro(MovingImageDimension) >
                                                    MovingGradientImageType;
  typedef GradientRecursiveGaussianImageFilter< MovingImageType,
                                                MovingGradientImageType >
                                                 MovingGradientImageFilterType;

  /** Image gradient calculator types. */
  typedef CentralDifferenceImageFunction<
                                FixedImageType,
                                CoordinateRepresentationType>
                                               FixedGradientCalculatorType;
  typedef CentralDifferenceImageFunction<
                                MovingImageType,
                                CoordinateRepresentationType>
                                               MovingGradientCalculatorType;

  /**  Type of the measure. */
  typedef typename Superclass::MeasureType    MeasureType;

  /**  Type of the derivative. */
  typedef typename Superclass::DerivativeType DerivativeType;

  /* Set/get images */
  /** Connect the Fixed Image.  */
  itkSetConstObjectMacro(FixedImage, FixedImageType);
  /** Get the Fixed Image. */
  itkGetConstObjectMacro(FixedImage, FixedImageType);
  /** Connect the Moving Image.  */
  itkSetConstObjectMacro(MovingImage, MovingImageType);
  /** Get the Moving Image. */
  itkGetConstObjectMacro(MovingImage, MovingImageType);

  /** Define the virtual reference space. This space defines the resolution,
   *  at which the registration is performed as well as the physical coordinate
   *  system.  Useful for unbiased registration.  If the user does not set this
   *  explicitly then it is taken from the fixed image in \c Initialize method.
   *  This method will allocate \c m_VirtualDomainImage with the passed
   *  information, and set \c m_VirtualDomainRegion, the region over which
   *  metric evaluation is performed. The region can also be set separately
   *  using \c SetVirtualDomainRegion, which replaces the region information set
   *  with this method.
   *  To set all settings from an existing image, use \c SetVirtualDomainImage.
   */
  void CreateVirtualDomainImage( VirtualSpacingType &,
                                    VirtualOriginType &,
                                    VirtualDirectionType &,
                                    VirtualSizeType &,
                                    VirtualIndexType & );
  /** Alternate version that takes the region itself */
  void CreateVirtualDomainImage( VirtualSpacingType &,
                                    VirtualOriginType &,
                                    VirtualDirectionType &,
                                    VirtualRegionType & );

  const VirtualSpacingType GetVirtualDomainSpacing( void );
  const VirtualOriginType GetVirtualDomainOrigin( void );
  const VirtualDirectionType GetVirtualDomainDirection( void );
  /* It gets messy to have individual set accessors for the domain
   * settings. Have to track if each has been set before initializing.
  itkSetMacro(VirtualDomainSpacing, VirtualSpacingType);
  itkGetMacro(VirtualDomainSpacing, VirtualSpacingType);
  itkSetMacro(VirtualDomainOrigin, VirtualOriginType);
  itkGetMacro(VirtualDomainOrigin, VirtualOriginType);
  itkSetMacro(VirtualDomainDirection, VirtualDirectionType);
  itkGetMacro(VirtualDomainDirection, VirtualDirectionType);
  */
  /** Set the size of the virtual domain region over which to compute
   * the metric. Any subsequent call to SetVirtualRegion will override
   * these values. */
  //void SetVirtualDomainSize( VirtualSizeType& index );
  //const VirtualSizeType GetVirtualDomainSize( void );

  /** Set/Get the index of the virtual domain region over which to compute
   * the metric. Any subsequent call to SetVirtualRegion will override
   * these values. */
  //void SetVirtualDomainIndex( VirtualDomainIndex& index );
  //const VirtualIndexType GetVirtualDomainIndex( void );

  /** Set/Get the region of the virtual domain over which to compute
   * the metric. If not set by user either here or by
   * CreateVirtualDomainImage, or by SetVirtualDomainImage, it will be set
   * from m_FixedImage.GetRequestedRegion() in Intitialize.
   * The region passed here will be replaced by
   * any subsequent calls to CreateVirtualDomainImage. */
  // FIXME: handle const-ness properly
  void SetVirtualDomainRegion( VirtualRegionType & region )
    { SetVirtualDomainImage( static_cast<const VirtualRegionType& >(region) ); }
  void SetVirtualDomainRegion( const VirtualRegionType & region )
    {
    if( region != m_VirtualDomainRegion || ! m_VirtualDomainRegionHasBeenSet )
      {
      m_VirtualDomainRegionHasBeenSet = true;
      m_VirtualDomainRegion = region;
      this->Modified();
      }
    }
  itkGetMacro(VirtualDomainRegion, VirtualRegionType);

  /** Set all virtual domain setings at once via an image.
   * The image is expected to have already been allocated.
   * The image's requested region will be used for VirtualDomainRegion,
   * which can be then changed by calling SetVirtualDomainRegion. */
  virtual void SetVirtualDomainImage(VirtualImageType* image);
  /** Get the virtual domain image */
  itkGetConstObjectMacro(VirtualDomainImage, VirtualImageType);

  /** Connect the fixed transform. */
  itkSetObjectMacro(FixedImageTransform, FixedTransformType);
  /** Get a pointer to the fixed transform.  */
  itkGetConstObjectMacro(FixedImageTransform, FixedTransformType);
  /** Connect the moving transform. */
  itkSetObjectMacro(MovingImageTransform, MovingTransformType);
  /** Get a pointer to the moving transform.  */
  itkGetConstObjectMacro(MovingImageTransform, MovingTransformType);

  /** Connect the fixed interpolator. */
  itkSetObjectMacro(FixedInterpolator, FixedInterpolatorType);
  /** Get a pointer to the fixed interpolator.  */
  itkGetConstObjectMacro(FixedInterpolator, FixedInterpolatorType);

  /** Connect the Moving interpolator. */
  itkSetObjectMacro(MovingInterpolator, MovingInterpolatorType);
  /** Get a pointer to the Moving interpolator.  */
  itkGetConstObjectMacro(MovingInterpolator, MovingInterpolatorType);

  /** Set/Get the moving image mask. */
  itkSetObjectMacro(MovingImageMask, MovingImageMaskType);
  itkSetConstObjectMacro(MovingImageMask, MovingImageMaskType);
  itkGetConstObjectMacro(MovingImageMask, MovingImageMaskType);

  /** Set/Get the fixed image mask. */
  itkSetObjectMacro(FixedImageMask, FixedImageMaskType);
  itkSetConstObjectMacro(FixedImageMask, FixedImageMaskType);
  itkGetConstObjectMacro(FixedImageMask, FixedImageMaskType);

  /** Set/Get gradient computation via GradientRecursiveGaussianImageFilter. */
  itkSetMacro(PrecomputeGradient, bool);
  itkGetConstReferenceMacro(PrecomputeGradient, bool);
  itkBooleanMacro(PrecomputeGradient);

  /** Get number of threads to use.
   * \warning This value can change during Initialize, if the threader
   * determines that fewer threads would be more efficient. The value is
   * initialized in constructor to the default value returned by the
   * assigned threader object, and will never be changed in Initialize to
   * a larger value. Derived classes
   * should only use this value once Superclass::Initialize() has been called.
   */
  itkGetConstMacro( NumberOfThreads, ThreadIdType );

  /** Set number of threads to use.
   * \warning See discussion for GetNumberOfThreads. */
  void SetNumberOfThreads( ThreadIdType );

  /** Computes the gradients of the fixed and moving images, using the
   * GradientRecursiveGaussianImageFilter, assigning the output to
   * to m_[Fixed|Moving]GaussianGradientImage. */
  virtual void ComputeGaussianGradient(void);

  /** Get Fixed Gradient Image. */
  itkGetConstObjectMacro(FixedGaussianGradientImage, FixedGradientImageType);
  /** Get Moving Gradient Image. */
  itkGetConstObjectMacro(MovingGaussianGradientImage, MovingGradientImageType);

  /** Get number of valid points from most recent update */
  itkGetConstMacro( NumberOfValidPoints, SizeValueType );

  /** Get the measure value, as computed by most recent evaluation */
  itkGetConstMacro( Value, MeasureType );

  /** Return the number of parameters.
   * NOTE: this is required because we derive from CostFunction. But
   * it's not really relevant, and may be removed if we decide not
   * to derive from SingleValuedCostFunction. */
  /* We return the number of parameters
   * in moving image transform because (for now at least) that's the
   * transform we're optimizing. */
  virtual unsigned int GetNumberOfParameters() const
    { return this->m_MovingImageTransform->GetNumberOfParameters(); }

  /** Return the size of derivative, taking DerivativeSource
   * into account. */
  virtual SizeValueType GetLocalDerivativeSize() const;

  /* Initialize the metric before calling GetValue or GetDerivative.
   * Derived classes must call this Superclass version if they override
   * this to perform their own initialization.
   * \note This is meant to be called once for a particular metric setup.
   * That is, when used in registration, this method would be called once
   * before entering the registration loop, during which GetValue or
   * GetDerivative will be called repeatedly. */
  virtual void Initialize(void) throw ( itk::ExceptionObject );

protected:

  /** Transform a point from VirtualImage domain to FixedImage domain.
   * This function also checks if mapped point is within the mask if
   * one is set, and that is within the moving image buffer, in which
   * case \c pointIsValid will be true on return.
   * \c mappedFixedPoint and \c fixedImageValue are also returned, and
   * \c fixedImageGradient is returned if \c computeImageGradient is set,
   * and it is mapped into the virtual space.
   * \note It would be better for maintainence to have a single method
   * that could work for either fixed or moving domains. However setting
   * that up is complicated because dimensionality and pixel type may
   * be different between the two. */
  virtual void TransformAndEvaluateFixedPoint(const VirtualPointType & point,
                              FixedImagePointType & mappedFixedPoint,
                              bool & pointIsValid,
                              FixedImagePixelType & fixedImageValue,
                              bool computeImageGradient,
                              FixedImageDerivativesType & fixedGradient,
                              ThreadIdType threadID) const;

  /** Transform a point from VirtualImage domain to MovingImage domain,
   * as is done in \c TransformAndEvaluateMovingPoint. */
  virtual void TransformAndEvaluateMovingPoint(const VirtualPointType & point,
                              MovingImagePointType & mappedMovingPoint,
                              bool & pointIsValid,
                              MovingImagePixelType & movingImageValue,
                              bool computeImageGradient,
                              MovingImageDerivativesType & movingGradient,
                              ThreadIdType threadID) const;

  /** Compute image derivatives at a point. */
  virtual void ComputeFixedImageDerivatives(
                                    const FixedImagePointType & mappedPoint,
                                    FixedImageDerivativesType & gradient,
                                    ThreadIdType threadID) const;
  virtual void ComputeMovingImageDerivatives(
                                    const MovingImagePointType & mappedPoint,
                                    MovingImageDerivativesType & gradient,
                                    ThreadIdType threadID) const;

  /** Store derivative result from a single point calculation. */
  virtual void StoreDerivativeResult(  DerivativeType & derivative,
                                        const VirtualIndexType & virtualIndex,
                                        ThreadIdType threadID );

  /** Perform any initialization required before each evaluation of
   * value and derivative. This is distinct from Initialize, which
   * is called only once before a number of iterations, e.g. before
   * a registration loop.
   * Called from \c GetValueAndDerivativeMultiThreadedInitiate before
   * threading starts. */
  virtual void InitializeForIteration(void);

  /** Called from \c GetValueAndDerivativeMultiThreadedInitiate after
   * threading is complete, to count the total number of valid points
   * used during calculations, storing it in \c m_NumberOfValidPoints */
  virtual void CollectNumberOfValidPoints(void);

  FixedImageConstPointer  m_FixedImage;
  FixedTransformPointer   m_FixedImageTransform;
  MovingImageConstPointer m_MovingImage;
  MovingTransformPointer  m_MovingImageTransform;
  VirtualImagePointer     m_VirtualDomainImage;

  VirtualRegionType       m_VirtualDomainRegion;
  bool                    m_VirtualDomainRegionHasBeenSet;
  //VirtualSpacingType      m_VirtualDomainSpacing;
  //VirtualOriginType       m_VirtualDomainOrigin;
  //VirtualDirectionType    m_VirtualDomainDirection;

  /** Pointers to interpolators */
  FixedInterpolatorPointer                        m_FixedInterpolator;
  MovingInterpolatorPointer                       m_MovingInterpolator;
  /** Pointers to BSplineInterpolators */
  typename FixedBSplineInterpolatorType::Pointer  m_FixedBSplineInterpolator;
  typename MovingBSplineInterpolatorType::Pointer m_MovingBSplineInterpolator;
  /** Boolean to indicate if the fixed interpolator is BSpline. */
  bool                               m_FixedInterpolatorIsBSpline;
  /** Boolean to indicate if the moving interpolator is BSpline. */
  bool                               m_MovingInterpolatorIsBSpline;

  /** Flag to control use of precomputed Gaussian gradient filter or gradient
   * calculator for image gradient calculations. */
  bool                               m_PrecomputeGradient;
  /** Gradient images to store Gaussian gradient filter output. */
  typename FixedGradientImageType::Pointer    m_FixedGaussianGradientImage;
  typename MovingGradientImageType::Pointer   m_MovingGaussianGradientImage;

  /** Image gradient calculators */
  typename FixedGradientCalculatorType::Pointer   m_FixedGradientCalculator;
  typename MovingGradientCalculatorType::Pointer  m_MovingGradientCalculator;

  /** Derivative results holder. */
  DerivativeType                              m_DerivativeResult;
  /** Store the number of points used during most recent value and derivative
   * calculation. */
  SizeValueType                               m_NumberOfValidPoints;

  /** Masks */
  FixedImageMaskConstPointer                  m_FixedImageMask;
  MovingImageMaskConstPointer                 m_MovingImageMask;

  /** Metric value, stored after evaluating */
  MeasureType             m_Value;

  /** User worker method to calculate value and derivative
   * given the point, value and image derivative for both fixed and moving
   * spaces. The provided values have been calculated from \c virtualPoint,
   * which is provided in case it's needed.
   * Must be overriden by derived classes to do anything meaningful.
   * \c mappedMovingPoint and \c mappedFixedPoint will be valid points
   * that have passed bounds checking, and lie within any mask that may
   * be assigned.
   * Results are returned in \c metricValueResult and \c derivativesResult, and
   * will be managed by this base class.
   * \c threadID may be used as needed, for example to access any per-thread
   * data cached during pre-processing by the derived class.
   * \warning The derived class should use \c m_NumberOfThreads from this base
   * class only after ImageToImageObjectMetric:: Initialize has been called, to
   * assure that the same number of threads are used.
   * \warning  This is called from the threader, and thus must be thread-safe.
   */
  virtual inline bool GetValueAndDerivativeProcessPoint(
        const VirtualPointType &           itkNotUsed(virtualPoint),
        const FixedImagePointType &        itkNotUsed(mappedFixedPoint),
        const FixedImagePixelType &        itkNotUsed(fixedImageValue),
        const FixedImageDerivativesType &  itkNotUsed(FixedImageDerivatives),
        const MovingImagePointType &       itkNotUsed(mappedMovingPoint),
        const MovingImagePixelType &       itkNotUsed(movingImageValue),
        const MovingImageDerivativesType & itkNotUsed(movingImageDerivatives),
        MeasureType &                      itkNotUsed(metricValueResult),
        DerivativeType &                   itkNotUsed(derivativesResult),
        ThreadIdType                       itkNotUsed(threadID) )
    /* ImageToImageMetric had this method simply return false. But why not
     * make it pure virtual ? */
    {return false;}

  /*
   * Multi-threading variables and methods
   */

  /** Initiates multi-threading to calculate the current metric value
   * and derivatives.
   * Derived classes should call this in their GetValueAndDerivative method.
   * This will end up calling the derived class'
   * GetValueAndDerivativeProcessPoint for each valid point in
   * VirtualDomainRegion.
   * See ... */
  virtual void GetValueAndDerivativeMultiThreadedInitiate();

  /** Default post-processing after multi-threaded calculation of
   * value and derivative. Typically called by derived classes after
   * GetValueAndDerivativeMultiThreadedInitiate. Collects the results
   * from each thread and sums them.
   * Results are stored in m_Value and m_DerivativeResult on return.
   * Pass true for \c doAverage to use the number of valid points, \c
   * m_NumberOfValidPoints, to average the value sum, and to average derivative
   * sums for global transforms only (i.e. transforms without local support).
   * Derived classes need not call this if they require special handling.
   */
  virtual void GetValueAndDerivativeMultiThreadedPostProcess( bool doAverage );

  /** Type of the default threader used for GetValue and GetDerivative.
   * This splits an image region in per-thread sub-regions over the outermost
   * image dimension. */
  typedef ImageToData< VirtualImageDimension, Self >
                                             ValueAndDerivativeThreaderType;
  typedef typename ValueAndDerivativeThreaderType::InputObjectType
                                             ThreaderInputObjectType;

  /* Optinally set the threader type to use. This performs the splitting of the
   * virtual region over threads, and user may wish to provide a different
   * one that does a different split. The default is ImageToData. */
  itkSetObjectMacro(ValueAndDerivativeThreader,ValueAndDerivativeThreaderType);

  /** Threader used to evaluate value and deriviative. */
  typename ValueAndDerivativeThreaderType::Pointer
                                              m_ValueAndDerivativeThreader;


  /** Intermediary threaded metric value storage. */
  std::vector<InternalComputationValueType>   m_MeasurePerThread;
  std::vector< DerivativeType >               m_DerivativesPerThread;
  std::vector< SizeValueType >                m_NumberOfValidPointsPerThread;

  ImageToImageObjectMetric();
  virtual ~ImageToImageObjectMetric() {}

  void PrintSelf(std::ostream& os, Indent indent) const
    {
    Superclass::PrintSelf(os, indent);
    os << indent << "ImageToImageObjectMetric: TODO..." << std::endl;
    }


private:

  /** Multi-threader callback used to iterate over image region by thread,
   * and call the derived class' user worker method to calculate
   * value and derivative.
   * If a derived class needs to implement its own callback to replace this,
   * define a static method with a different name, and assign it to the
   * threader in the class' constructor by calling
   * \c m_ValueAndDerivativeThreader->SetThreadedGenerateData( mycallback ) */
  static void GetValueAndDerivativeMultiThreadedCallback(
                          const ThreaderInputObjectType& virtualImageSubRegion,
                          ThreadIdType threadId,
                          Self * dataHolder);

  /* The number of threads to use.
   * /warning See discussion in GetNumberOfThreads.
   */
  ThreadIdType            m_NumberOfThreads;

  //purposely not implemented
  ImageToImageObjectMetric(const Self &);
  //purposely not implemented
  void operator=(const Self &);

};
}//namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkImageToImageObjectMetric.txx"
#endif

#endif
