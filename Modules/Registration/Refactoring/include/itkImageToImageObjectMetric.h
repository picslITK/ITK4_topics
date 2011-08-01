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
#include "itkDeformationFieldTransform.h"
#include "itkCompositeTransform.h"
#include "itkIdentityTransform.h"
#include "itkResampleImageFilter.h"

namespace itk
{
/** \class ImageToImageObjectMetric
 *
 * Computes similarity between regions of two images, using two
 * user-supplied transforms.
 *
 * \note Currently support for using only \c GetValueAndDerivative is
 * implemented. \c GetValue
 * will follow after final implementation details are worked out.
 *
 * Templated over the fixed and moving image types, as well as an optional
 * VirtualImage type to define the virtual domain. The VirtualImage type
 * defaults to TFixedImage.
 *
 * This class uses a virtual reference space. This space defines the resolution,
 *  at which the registration is performed as well as the physical coordinate
 *  system.  Useful for unbiased registration.  If the user does not set this
 *  explicitly then it is taken from the fixed image in \c Initialize method.
 * The user can define a virtual domain by calling either
 * \c CreateVirtualDomainImage or \c SetVirtualDomainImage. The virtual
 * domain region can be changed via \c SetVirtualDomainRegion. See these
 * methods for details.
 *
 * Both transforms are initialized to an IdentityTransform.
 *
 * The PreWarpImages option can be set for computational efficiency. This
 * creates a warped version for each image at the begin of each iteration.
 * However the cost is more memory usage (2x VirtualDomain size).
 *
 * \warning When using PreWarpImages flag, it is assumed that the images are
 * already in the virtual domain, for example via an initial affine
 * transformation in the first step of a CompositeTransform.
 *
 * \note Use of PreWarpImages option is not yet supported when also using
 * image masks.
 *
 * Image gradient calculation is performed in one of these ways:
 * If \c PreWarpImages is enabled:
 *  CentralDifferenceImageFunction.
 *  Other options will be supported in the future.
 * If \c PreWarpImages is disabled:
 *  1) The BSplineInterpolator is used if it has been set by user as the
 *  interpolator.
 *  2) otherwise, GradientRecursiveGaussianImageFilter is used by default to
 *  precompute the gradient for each image. This requires memory to hold the
 *  fully computed gradient image.
 *  3) lastly, CentralDifferenceImageFunction is used if
 *  m_PrecomputeImageGradient option has been set to false.
 *
 * Image masks are supported using SetMovingImageMask or SetFixedImageMask.
 * Random sampling or user-supplied point lists are not supported except
 * via a user-supplied mask.
 *
 * Derived classes need to override at least:
 *  GetValueAndDerivativeProcessPoint
 *  Pure virtual methods declared in the base class.
 *
 * See ImageToImageObjectMetricTest for a clear example of what a
 * derived class must implement and do.
 *
 * This class also implements EstimateScales and ComputeMaximumVoxelShift.
 * These methods depends on the fixed image, moving image, virtual image and
 * the transform objects. Therefore, it is natural to put them here.
 *
 * \note: EstimateScales and ComputeMaximumVoxelShift might be put into a
 * separate class for two reasons: 1) Users might want to provide different
 * strategies; 2) There is more flexibility to manipulate transform objects
 * if the information about specific tranform classes is available.
 *
 * \ingroup ITKRegistrationRefactoring
 */
template<class TFixedImage,class TMovingImage,class TVirtualImage = TFixedImage>
class ITK_EXPORT ImageToImageObjectMetric :
  public ObjectToObjectMetric
{
public:

  /** Standard class typedefs. */
  typedef ImageToImageObjectMetric                          Self;
  typedef ObjectToObjectMetric                              Superclass;
  typedef SmartPointer<Self>                                Pointer;
  typedef SmartPointer<const Self>                          ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro(ImageToImageObjectMetric, ObjectToObjectMetric);

  /** Type used internally for computations */
  typedef typename Superclass::InternalComputationValueType
                                                  InternalComputationValueType;

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
    itkGetStaticConstMacro( VirtualImageDimension )> MovingTransformType;

  typedef Transform<CoordinateRepresentationType,
    itkGetStaticConstMacro( FixedImageDimension ),
    itkGetStaticConstMacro( VirtualImageDimension )> FixedTransformType;

  /** Deformation field typedef for convenience */
  typedef DeformationFieldTransform<CoordinateRepresentationType,
    itkGetStaticConstMacro( MovingImageDimension ) >
                                          MovingDeformationFieldTransformType;
  /** CompositeTransform typedef for convenience */
  typedef CompositeTransform<CoordinateRepresentationType,
    itkGetStaticConstMacro( MovingImageDimension ) >
                                          MovingCompositeTransformType;

  /** Identity transform typedef's for convenience */
  typedef IdentityTransform<CoordinateRepresentationType,
    itkGetStaticConstMacro( MovingImageDimension ) >
                                          MovingIdentityTransformType;
  typedef IdentityTransform<CoordinateRepresentationType,
    itkGetStaticConstMacro( MovingImageDimension ) >
                                          FixedIdentityTransformType;

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

  typedef   CovariantVector< CoordinateRepresentationType,
                             itkGetStaticConstMacro(VirtualImageDimension) >
                                                      VirtualImageDerivativesType;

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

  /** ResampleImageFilter types for image pre-warping */
  typedef ResampleImageFilter< MovingImageType,
                               VirtualImageType,
                               MovingRealType >
                                             MovingWarpResampleImageFilterType;
  typedef typename MovingWarpResampleImageFilterType::Pointer
                                          MovingWarpResampleImageFilterPointer;
  typedef ResampleImageFilter< FixedImageType,
                               VirtualImageType,
                               FixedRealType >
                                             FixedWarpResampleImageFilterType;
  typedef typename FixedWarpResampleImageFilterType::Pointer
                                          FixedWarpResampleImageFilterPointer;

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

  /** Set/Get the region of the virtual domain over which to compute
   * the metric. If not set by user either here or by
   * CreateVirtualDomainImage, or by SetVirtualDomainImage, it will be set
   * from m_FixedImage.GetBufferedRegion() in Intitialize.
   * The region passed here will be replaced by
   * any subsequent calls to CreateVirtualDomainImage. */
  void SetVirtualDomainRegion( VirtualRegionType & region )
  // FIXME: handle const-ness properly
  { SetVirtualDomainRegion( static_cast<const VirtualRegionType& >(region) ); }

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
   * The image's buffered region will be used for VirtualDomainRegion,
   * which can be then changed by calling SetVirtualDomainRegion. */
  virtual void SetVirtualDomainImage(VirtualImageType* image);
  /** Get the virtual domain image */
  itkGetConstObjectMacro(VirtualDomainImage, VirtualImageType);

  /** Connect the fixed transform. */
  itkSetObjectMacro(FixedTransform, FixedTransformType);
  /** Get a pointer to the fixed transform.  */
  itkGetConstObjectMacro(FixedTransform, FixedTransformType);
  /** Connect the moving transform. */
  itkSetObjectMacro(MovingTransform, MovingTransformType);
  /** Get a pointer to the moving transform.  */
  itkGetConstObjectMacro(MovingTransform, MovingTransformType);
  /** Connect the moving transform using a backwards-compatible name */
  void SetTransform( MovingTransformType* transform )
  { SetMovingTransform( transform ); }
  /** Get the moving transform using a backwards-compatible name */
  const MovingTransformType * GetTransform()
  { return GetMovingTransform(); }

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
  itkSetMacro(PrecomputeImageGradient, bool);
  itkGetConstReferenceMacro(PrecomputeImageGradient, bool);
  itkBooleanMacro(PrecomputeImageGradient);

  /** Set/Get pre-warping of images option. */
  itkSetMacro(PreWarpImages, bool);
  itkGetConstReferenceMacro(PreWarpImages, bool);
  itkBooleanMacro(PreWarpImages);

  /** Get pre-warmed images */
  itkGetConstObjectMacro( MovingWarpedImage, MovingImageType );
  itkGetConstObjectMacro( FixedWarpedImage, FixedImageType );

  /** Get number of threads to use.
   * \warning This value can change during Initialize, if the threader
   * determines that fewer threads would be more efficient. The value is
   * initialized in constructor to the default value returned by the
   * assigned threader object, and will never be changed in Initialize to
   * a larger value. Derived classes
   * should only use this value once Superclass::Initialize() has been called.
   */
  itkGetConstMacro( NumberOfThreads, ThreadIdType );

  /** Set number of threads to use. The default is the number of threads
   * available as reported by MultiThreader.
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

  /** Get the measure value in m_Value, as computed by most recent evaluation.
   * Need to differentiate it from GetValue method in base class. */
  MeasureType GetValueResult()
  { return m_Value; }

  /** Return the number of parameters. Convenience methods for derived
   * classes mainly.
   * We're always optimizing the moving image transform, so return
   * its number of parameters. */
  virtual unsigned int GetNumberOfParameters() const
    { return this->m_MovingTransform->GetNumberOfParameters(); }

  /** Get a const reference to the moving transform's parameters */
  virtual const ParametersType & GetParameters() const
  { return this->m_MovingTransform->GetParameters(); }

  /** Update the metric's transform parameters.
   * This call is passed through directly to the transform.
   * \c factor is a scalar multiplier for each value in update, and
   * defaults to 1.0 .
   * \c derivative must be the proper size, as retrieved
   * from GetNumberOfParameters. */
  virtual void UpdateTransformParameters( DerivativeType & derivative,
                                          ParametersValueType factor = 1.0);

  /** FIXME: documentation. See GetNumberOfParameters */
  virtual unsigned int GetNumberOfLocalParameters() const
    { return this->m_MovingTransform->GetNumberOfLocalParameters(); }

  /** FIXME: documentation. See GetNumberOfParameters */
  virtual bool HasLocalSupport() const
    { return this->m_MovingTransform->HasLocalSupport(); }

  /** Return the size of derivative result, taking DerivativeSource
   * into account. */
  /* Actually I think this is wrong. Local derivate is always relevant to
   * moving transform since we're always optimizing moving transform. And
   * DerivativeSource is only for source of image derivatives. */
  //virtual SizeValueType GetLocalDerivativeSize() const;

  /* Initialize the metric before calling GetValue or GetDerivative.
   * Derived classes must call this Superclass version if they override
   * this to perform their own initialization.
   * \note This is meant to be called once for a particular metric setup.
   * That is, when used in registration, this method would be called once
   * before entering the registration loop, during which GetValue or
   * GetDerivative will be called repeatedly. It must be called again if
   * metric settings are changed before beginning a new registration. */
  virtual void Initialize(void) throw ( itk::ExceptionObject );

protected:

  /** User worker method to calculate value and derivative
   * given the point, value and image derivative for both fixed and moving
   * spaces. The provided values have been calculated from \c virtualPoint,
   * which is provided in case it's needed.
   * Must be overriden by derived classes to do anything meaningful.
   * \c mappedMovingPoint and \c mappedFixedPoint will be valid points
   * that have passed bounds checking, and lie within any mask that may
   * be assigned.
   * Results are returned in \c metricValueReturn and \c localDerivativeReturn,
   * and will be processed by this base class.
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
        MeasureType &                      itkNotUsed(metricValueReturn),
        DerivativeType &                   itkNotUsed(localDerivativeReturn),
        ThreadIdType                       itkNotUsed(threadID) )
    /* ImageToImageMetric had this method simply return false. But why not
     * make it pure virtual ? */
    {return false;}

  /** Perform any initialization required before each evaluation of
   * value and derivative. This is distinct from Initialize, which
   * is called only once before a number of iterations, e.g. before
   * a registration loop.
   * Called from \c GetValueAndDerivativeMultiThreadedInitiate before
   * threading starts. */     //NOTE: make this private?
  virtual void InitializeForIteration(void);

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
  virtual void TransformAndEvaluateFixedPoint(
                              const VirtualIndexType & index,
                              const VirtualPointType & point,
                              FixedImagePointType & mappedFixedPoint,
                              bool & pointIsValid,
                              FixedImagePixelType & fixedImageValue,
                              bool computeImageGradient,
                              FixedImageDerivativesType & fixedGradient,
                              ThreadIdType threadID) const;

  /** Transform a point from VirtualImage domain to MovingImage domain,
   * as is done in \c TransformAndEvaluateMovingPoint. */
  virtual void TransformAndEvaluateMovingPoint(
                              const VirtualIndexType & index,
                              const VirtualPointType & point,
                              MovingImagePointType & mappedMovingPoint,
                              bool & pointIsValid,
                              MovingImagePixelType & movingImageValue,
                              bool computeImageGradient,
                              MovingImageDerivativesType & movingGradient,
                              ThreadIdType threadID) const;

  /** When using pre-warped images, this routine will return pixel value
   * and optionally the image gradient at a given \c index.
   * Used internally.
   * \c virtualPoint is passed in for efficiency because it's already computed
   * before calling this method.
   * \c fixedImageValue and \c movingImageValue are returned.
   * \c fixedGradient and \c movingGradient are returned if
   * \c computeImageGradient is true.
   * \c pointIsValid is returned true if the index is valid within a mask if
   * one is supplied.
   * \warning Use of masks with pre-warped images is not yet implemented. */
  virtual void EvaluatePreWarpedImagesAtIndex( const VirtualIndexType & index,
                                 const VirtualPointType & virtualPoint,
                                 const bool computeImageGradient,
                                 FixedImagePixelType & fixedImageValue,
                                 MovingImagePixelType & movingImageValue,
                                 FixedImageDerivativesType & fixedGradient,
                                 MovingImageDerivativesType & movingGradient,
                                 bool & pointIsValid,
                                 const ThreadIdType threadID ) const;

  /** Compute image derivatives at a point. */
  virtual void ComputeFixedImageDerivatives(
                                    const FixedImagePointType & mappedPoint,
                                    FixedImageDerivativesType & gradient,
                                    ThreadIdType threadID) const;
  virtual void ComputeMovingImageDerivatives(
                                    const MovingImagePointType & mappedPoint,
                                    MovingImageDerivativesType & gradient,
                                    ThreadIdType threadID) const;

  /** Compute image derivatives at an index. */
  virtual void ComputeFixedImageDerivativesAtIndex(
                                    const VirtualIndexType & index,
                                    FixedImageDerivativesType & gradient,
                                    ThreadIdType threadID) const;
  virtual void ComputeMovingImageDerivativesAtIndex(
                                    const VirtualIndexType & index,
                                    MovingImageDerivativesType & gradient,
                                    ThreadIdType threadID) const;

  /** Store derivative result from a single point calculation.
   * \warning If this method is overridden or otherwise not used
   * in a derived class, be sure to *accumulate* results in
   * \c derivative, and not assign them. */
  virtual void StoreDerivativeResult(  DerivativeType & derivative,
                                        const VirtualIndexType & virtualIndex,
                                        ThreadIdType threadID );

  /** Called from \c GetValueAndDerivativeMultiThreadedInitiate after
   * threading is complete, to count the total number of valid points
   * used during calculations, storing it in \c m_NumberOfValidPoints */
  virtual void CollectNumberOfValidPoints(void);

  FixedImageConstPointer  m_FixedImage;
  FixedTransformPointer   m_FixedTransform;
  MovingImageConstPointer m_MovingImage;
  MovingTransformPointer  m_MovingTransform;
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
  bool                               m_PrecomputeImageGradient;

  /** Gradient images to store Gaussian gradient filter output. */
  typename FixedGradientImageType::Pointer    m_FixedGaussianGradientImage;
  typename MovingGradientImageType::Pointer   m_MovingGaussianGradientImage;

  /** Image gradient calculators */
  typename FixedGradientCalculatorType::Pointer   m_FixedGradientCalculator;
  typename MovingGradientCalculatorType::Pointer  m_MovingGradientCalculator;

  /** Flag to control pre-warping of images. */
  bool                               m_PreWarpImages;

  /** Pre-warped images. */
  FixedImagePointer   m_FixedWarpedImage;
  MovingImagePointer  m_MovingWarpedImage;

  /** Resample image filters for pre-warping images */
  MovingWarpResampleImageFilterPointer    m_MovingWarpResampleImageFilter;
  FixedWarpResampleImageFilterPointer     m_FixedWarpResampleImageFilter;

  /** Derivative results holder. User a raw pointer so we can point it
   * to a user-provided object. This enables
   * safely sharing a derivative object between metrics during multi-variate
   * analsys, for memory efficiency. */
  DerivativeType*                             m_DerivativeResult;

  /** Store the number of points used during most recent value and derivative
   * calculation. */
  SizeValueType                               m_NumberOfValidPoints;

  /** Masks */
  FixedImageMaskConstPointer                  m_FixedImageMask;
  MovingImageMaskConstPointer                 m_MovingImageMask;

  /** Metric value, stored after evaluating */
  MeasureType             m_Value;

  /*
   * Multi-threading variables and methods
   */

  /** Initialize memory for threading.
   * \c derivativeReturn will be used to store the derivative results.
   * Typically this will be the user-supplied object from a call to
   * GetValueAndDerivative.
   */ //NOTE: make this private, or will a derived class want to override?
  virtual void InitializeThreadingMemory( DerivativeType & derivativeReturn );

  /** Initiates multi-threading to evaluate the current metric value
   * and derivatives.
   * Derived classes should call this in their GetValueAndDerivative method.
   * This will end up calling the derived class'
   * GetValueAndDerivativeProcessPoint for each valid point in
   * VirtualDomainRegion.
   * Pass in \c derivativeReturn from user. Results are written directly
   * into this parameter.
   * \sa GetValueAndDerivativeMultiThreadedPostProcess
   */
  virtual void GetValueAndDerivativeMultiThreadedInitiate( DerivativeType &
                                                            derivativeReturn );

  /** Default post-processing after multi-threaded calculation of
   * value and derivative. Typically called by derived classes after
   * GetValueAndDerivativeMultiThreadedInitiate. Collects the results
   * from each thread and sums them.
   * Results are stored in \c m_Value and \c m_DerivativeResult.
   * \c m_DerivativeResult is set during initialization to point to the
   * user-supplied derivative parameter in GetValueAndDerivative. Thus,
   * the derivative results are written directly to this parameter.
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
  std::vector<InternalComputationValueType>  m_MeasurePerThread;
  std::vector< DerivativeType >              m_DerivativesPerThread;
  std::vector< DerivativeType >              m_LocalDerivativesPerThread;
  std::vector< SizeValueType >               m_NumberOfValidPointsPerThread;
  /** Pre-allocated transform jacobian objects, for use as needed by dervied
   * classes for efficiency. */
  std::vector< MovingTransformJacobianType>  m_MovingTransformJacobianPerThread;
  /** FIXME: may need separate types for fixed and moving, but ok for now. */
  std::vector<
    typename MovingDeformationFieldTransformType::AffineTransformPointer >
                                                m_AffineTransformPerThread;

  ImageToImageObjectMetric();
  virtual ~ImageToImageObjectMetric() {}

  void PrintSelf(std::ostream& os, Indent indent) const
    {
    Superclass::PrintSelf(os, indent);
    os << indent << "ImageToImageObjectMetric: TODO..." << std::endl;
    }


  /* Verify that virtual domain and deformation field are the same size
   * and in the same physical space. */
  virtual void VerifyDeformationFieldSizeAndPhysicalSpace();

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

  /** Pre-warp the images for efficiency and computational stability. */
  void PreWarpImages( void );

  /** Flag to track if threading memory has been initialized since last
   * call to Initialize. */
  bool                    m_ThreadingMemoryHasBeenInitialized;

  /* The number of threads to use.
   * /warning See discussion in GetNumberOfThreads.
   */
  ThreadIdType            m_NumberOfThreads;

  //purposely not implemented
  ImageToImageObjectMetric(const Self &);
  //purposely not implemented
  void operator=(const Self &);

  //Sample point coordinates from the virtual image domain
  std::vector<VirtualPointType> m_VirtualImageCornerPoints;
};
}//namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkImageToImageObjectMetric.hxx"
#endif

#endif
