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
#include "itkSpatialObject.h"
//#include "itkImageToData.h"
#include "itkIdentityTransform.h"
#include "itkResampleImageFilter.h"

namespace itk
{

//Forward-declare this because of module dependency conflict.
//ImageToData will soon be moved to a different module, at which
// time this can be removed.
template <unsigned int VDimension, class TDataHolder>
class ImageToData;

/** \class ImageToImageObjectMetric
 *
 * Computes similarity between regions of two images, using two
 * user-supplied transforms, a 'fixed' transform and a 'moving' transform.
 *
 * \note Currently support for using only \c GetValueAndDerivative is
 * implemented. \c GetValue will follow after final implementation details
 * are worked out.
 *
 * Templated over the fixed and moving image types, as well as an optional
 * VirtualImage type to define the virtual domain. The VirtualImage type
 * defaults to TFixedImage.
 *
 * This class uses a virtual reference space. This space defines the resolution
 * at which the registration is performed, as well as the physical coordinate
 * system.  Useful for unbiased registration.  If the user does not set this
 * explicitly then it is taken from the fixed image in \c Initialize method.
 * The user can define a virtual domain by calling either
 * \c CreateVirtualDomainImage or \c SetVirtualDomainImage. The virtual
 * domain region can be changed via \c SetVirtualDomainRegion. See these
 * methods for details.
 *
 * Both transforms are initialized to an IdentityTransform, and can be
 * set by the user using \c SetFixedTranform and \c SetMovingTransform.
 *
 * At a minimum, the user must:
 *  1) Set images using \c SetFixedImage and \c SetMovingImage.
 *  3) Call \c Initialize.
 *
 * Pre-warping:
 * The \c SetPreWarpFixedImage and \c SetPreWarpMovingImage options can be set
 * for better speed. When set, these create a warped version for each image at
 * the begin of each iteration, warping each image into the virtual domain.
 * However the cost is more memory usage (VirtualDomain size for each image).
 * The fixed image is only pre-warped once, during the call to \c Initialize,
 * because it is assumed the fixed transform is not changing. The moving image
 * is pre-warped at the begin of every iteration, because it is assumed the
 * moving transform is changing (e.g. during registration).
 * By default, pre-warping is enabled for both fixed and moving images.
 *
 * Image gradient caclulations:
 * Image gradients can be calculated in one of two ways:
 * 1) Using GradientRecursiveGaussianImageFilter, by setting
 *  \c Use[Fixed|Moving]GradientRecursiveGaussianImageFilter to true. This is a
 *  smoothed gradient filter. This filter uses more memory, because it
 *  calculates all gradients at once and stores them in an image. The advantage
 *  of pre-calculation is for the fixed image gradients, since they only need be
 *  calculated once, and for metrics that need to access image gradients more
 *  than once for a particular point. The fixed image gradients are only
 *  calculated once when this option is set, during \c Initialize.
 * 2) Otherwise, the CentralDifferenceImageFunction is used. This calculation
 *  is not smoothed and gives different results than
 *  GradientRecursiveGaussianImageFilter. The advantage is that less memory is
 *  used. However for the fixed image, it means needlessly computing the image
 *  gradients at each iteration of a registration instead of just computing
 *  once at the begin.
 *
 * Both image gradient calculation methods are threaded.
 * Generally it is not recommended to use different image gradient methods for
 * the fixed and moving images because the methods return different results.
 *
 * Image masks are supported using SetMovingImageMask or SetFixedImageMask.
 *
 * Random sampling or user-supplied point lists are not yet supported, except
 * via an image mask. If the mask is at all sparse, the
 * SetPreWarp[Fixed|Moving]Image and
 * Use[Fixed|Moving]GradientRecursiveGaussianImageFilter options should be
 * disabled.
 *
 * This class is threaded.
 *
 * Derived classes:
 *
 *  Derived classes need to override at least:
 *  \c GetValueAndDerivative
 *  Pure virtual methods declared in the base class.
 *
 *  \c GetValueAndDerivativeProcessPoint must be overridden by derived
 *  classes that use \c GetValueAndDerivativeMultiThreadedInitiate.
 *
 *  See \c ImageToImageObjectMetricTest for a clear example of what a
 *  derived class must implement and do. Pre- and Post-processing are
 *  handled by the derived class in its \c GetValueAndDerivative method, as
 *  described in \c ImageToImageObjectMetricTest.
 *
 * \ingroup ITKHighDimensionalMetrics
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

  /** Graident source type */
  typedef typename Superclass::GradientSourceType GradientSourceType;

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

  typedef typename MovingTransformType::Pointer         MovingTransformPointer;
  typedef typename MovingTransformType::InputPointType  MovingInputPointType;
  typedef typename MovingTransformType::OutputPointType MovingOutputPointType;
  typedef typename MovingTransformType::ParametersType
                                                MovingTransformParametersType;

  /** Jacobian type. This is the same for all transforms */
  typedef typename FixedTransformType::JacobianType     JacobianType;

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

  /** Image derivatives types */
  typedef   CovariantVector< CoordinateRepresentationType,
                             itkGetStaticConstMacro(FixedImageDimension) >
                                                      FixedImageGradientType;
  typedef   CovariantVector< CoordinateRepresentationType,
                             itkGetStaticConstMacro(MovingImageDimension) >
                                                      MovingImageGradientType;

  typedef   CovariantVector< CoordinateRepresentationType,
                             itkGetStaticConstMacro(VirtualImageDimension) >
                                                      VirtualImageGradientType;

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
  void SetVirtualDomainRegion( VirtualRegionType & region );
  void SetVirtualDomainRegion( const VirtualRegionType & region );
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

  /** Connect the moving transform using a backwards-compatible name.
   * This assigns the input transform to the moving transform. */
  void SetTransform( MovingTransformType* transform );

  /** Get the moving transform using a backwards-compatible name */
  const MovingTransformType * GetTransform();

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

  /** Set/Get gradient computation via GradientRecursiveGaussianImageFilter
   * for fixed image. */
  itkSetMacro(UseFixedGradientRecursiveGaussianImageFilter, bool);
  itkGetConstReferenceMacro(UseFixedGradientRecursiveGaussianImageFilter, bool);
  itkBooleanMacro(UseFixedGradientRecursiveGaussianImageFilter);

  /** Set/Get gradient computation via GradientRecursiveGaussianImageFilter. */
  itkSetMacro(UseMovingGradientRecursiveGaussianImageFilter, bool);
  itkGetConstReferenceMacro(UseMovingGradientRecursiveGaussianImageFilter, bool);
  itkBooleanMacro(UseMovingGradientRecursiveGaussianImageFilter);

  /** Set/Get pre-warping of fixed image option. */
  itkSetMacro(PreWarpFixedImage, bool);
  itkGetConstReferenceMacro(PreWarpFixedImage, bool);
  itkBooleanMacro(PreWarpFixedImage);

  /** Set/Get pre-warping of Moving image option. */
  itkSetMacro(PreWarpMovingImage, bool);
  itkGetConstReferenceMacro(PreWarpMovingImage, bool);
  itkBooleanMacro(PreWarpMovingImage);

  /** Get pre-warped images */
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

  /** Computes the gradients of the fixed image, using the
   * GradientRecursiveGaussianImageFilter, assigning the output to
   * to m_FixedGaussianGradientImage. It will use either the original
   * fixed image, or the pre-warped version, depending on the setting
   * of PreWarpFixedImage. */
  virtual void ComputeFixedGaussianGradientImage(void);

  /** Computes the gradients of the moving image, using the
   * GradientRecursiveGaussianImageFilter, assigning the output to
   * to m_MovingGaussianGradientImage. It will use either the original
   * moving image, or the pre-warped version, depending on the setting
   * of PreWarpMovingImage. */
  virtual void ComputeMovingGaussianGradientImage(void);

  /** Get Fixed Gradient Image. */
  itkGetConstObjectMacro(FixedGaussianGradientImage, FixedGradientImageType);
  /** Get Moving Gradient Image. */
  itkGetConstObjectMacro(MovingGaussianGradientImage, MovingGradientImageType);

  /** Get the gradient calculators */
  itkGetConstObjectMacro( FixedGradientCalculator, FixedGradientCalculatorType);
  itkGetConstObjectMacro( MovingGradientCalculator, MovingGradientCalculatorType);

  /** Get number of valid points from most recent update */
  itkGetConstMacro( NumberOfValidPoints, SizeValueType );

  /** Get the measure value in m_Value, as computed by most recent evaluation.
   * Need to differentiate it from GetValue method in base class. */
  MeasureType GetValueResult();

  /** Return the number of parameters. Convenience methods for derived
   * classes mainly.
   * We're always optimizing the moving image transform, so return
   * its number of parameters.
   * TODO: Swithc return type to NumberOfParametersType when superclass
   * has been changed. */
  virtual unsigned int GetNumberOfParameters() const;

  /** Get a const reference to the moving transform's parameters */
  virtual const ParametersType & GetParameters() const;

  /** Update the moving transform's parameters.
   * This call is passed through directly to the transform.
   * \c factor is a scalar multiplier for each value in update, and
   * defaults to 1.0 .
   * \c derivative must be the proper size, as retrieved
   * from GetNumberOfParameters. */
  virtual void UpdateTransformParameters( DerivativeType & derivative,
                                          ParametersValueType factor = 1.0);

  /** Get the number of local parameters from the moving transform. */
  virtual unsigned int GetNumberOfLocalParameters() const;

  /** Get if the moving transform has local support. */
  virtual bool HasLocalSupport() const;

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
   * By default, this will return false and not assign any values.
   * Derived classes that use \c GetValueAndDerivativeMultiThreadedInitiate
   * to initiate process must override this method.
   * \note This method is not pure virtual because some derived classes
   * do not use \c GetValueAndDerivativeMultiThreadedInitiate, and instead
   * provide their own processing control.
   * \c mappedMovingPoint and \c mappedFixedPoint will be valid points
   * that have passed bounds checking, and lie within any mask that may
   * be assigned.
   * \c mappedFixedPixelValue, \c mappedFixedImageGradient,
   * \c mappedMovingPixelValue, and \c mappedMovingImageGradient are
   * provided for use by the derived class. Note however, that the
   * mappedFixed* values are only calculated when
   * \c m_GradientSource is set to either \c *_FIXED or \c *_BOTH, and
   * similarly respectively for mappedMoving* values. Otherwise, the
   * values are meaningless and should be ignored.
   * Results are experected to be returned in \c metricValueReturn and
   * \c localDerivativeReturn, and will be processed by this base class.
   * \c threadID may be used as needed, for example to access any per-thread
   * data cached during pre-processing by the derived class.
   * \warning The derived class should use \c m_NumberOfThreads from this base
   * class only after ImageToImageObjectMetric:: Initialize has been called, to
   * assure that the same number of threads are used.
   * \warning  This is called from the threader, and thus must be thread-safe.
   */
  virtual bool GetValueAndDerivativeProcessPoint(
        const VirtualPointType &          virtualPoint,
        const FixedImagePointType &       mappedFixedPoint,
        const FixedImagePixelType &       mappedFixedPixelValue,
        const FixedImageGradientType &    mappedFixedImageGradient,
        const MovingImagePointType &      mappedMovingPoint,
        const MovingImagePixelType &      mappedMovingPixelValue,
        const MovingImageGradientType &   mappedMovingImageGradient,
        MeasureType &                     metricValueReturn,
        DerivativeType &                  localDerivativeReturn,
        const ThreadIdType                threadID );

  /** Perform any initialization required before each evaluation of
   * value and derivative. This is distinct from Initialize, which
   * is called only once before a number of iterations, e.g. before
   * a registration loop.
   * Called from \c GetValueAndDerivativeMultiThreadedInitiate before
   * threading starts. */     //NOTE: make this private or protected?
  virtual void InitializeForIteration(void);

  /** Transform and evaluate a point into the virtual domain.
   * This function also checks if mapped point is within the mask if
   * one is set, and that is within the fixed image buffer, in which
   * case \c pointIsValid will be true on return.
   * \c mappedFixedPoint and \c mappedFixedPixelValue are  returned, and
   * \c mappedFixedImageGradient is returned if \c computeImageGradient is set.
   * All return values are in the virtual domain.
   * \note It would be better for maintainence to have a single method
   * that could work for either fixed or moving domains. However setting
   * that up is complicated because dimensionality and pixel type may
   * be different between the two. */
  virtual void TransformAndEvaluateFixedPoint(
                           const VirtualIndexType & index,
                           const VirtualPointType & point,
                           const bool computeImageGradient,
                           FixedImagePointType & mappedFixedPoint,
                           FixedImagePixelType & mappedFixedPixelValue,
                           FixedImageGradientType & mappedFixedImageGradient,
                           bool & pointIsValid ) const;

  /** See TransformAndEvaluateFixedPoint. TODO. */
  virtual void TransformAndEvaluateMovingPoint(
                           const VirtualIndexType & index,
                           const VirtualPointType & point,
                           const bool computeImageGradient,
                           MovingImagePointType & mappedMovingPoint,
                           MovingImagePixelType & mappedMovingPixelValue,
                           MovingImageGradientType & mappedMovingImageGradient,
                           bool & pointIsValid ) const;

  /** Compute image derivatives at a point. */
  virtual void ComputeFixedImageGradient(
                                    const FixedImagePointType & mappedPoint,
                                    FixedImageGradientType & gradient ) const;
  virtual void ComputeMovingImageGradient(
                                    const MovingImagePointType & mappedPoint,
                                    MovingImageGradientType & gradient ) const;

  /** Compute image derivatives at an index. */
  virtual void ComputeFixedImageGradientAtIndex(
                                    const VirtualIndexType & index,
                                    FixedImageGradientType & gradient ) const;
  virtual void ComputeMovingImageGradientAtIndex(
                                    const VirtualIndexType & index,
                                    MovingImageGradientType & gradient ) const;

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

  /** Flag to control use of precomputed Gaussian gradient filter or gradient
   * calculator for image gradient calculations. */
  bool                          m_UseFixedGradientRecursiveGaussianImageFilter;
  bool                          m_UseMovingGradientRecursiveGaussianImageFilter;

  /** Gradient images to store Gaussian gradient filter output. */
  typename FixedGradientImageType::Pointer    m_FixedGaussianGradientImage;
  typename MovingGradientImageType::Pointer   m_MovingGaussianGradientImage;

  /** Image gradient calculators */
  typename FixedGradientCalculatorType::Pointer   m_FixedGradientCalculator;
  typename MovingGradientCalculatorType::Pointer  m_MovingGradientCalculator;

  /** Flag to control pre-warping of fixed image. */
  bool                               m_PreWarpFixedImage;

  /** Flag to control pre-warping of moving image. */
  bool                               m_PreWarpMovingImage;

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
  typedef ImageToData<VirtualImageDimension, Self>
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
  std::vector<JacobianType>                  m_MovingTransformJacobianPerThread;

  ImageToImageObjectMetric();
  virtual ~ImageToImageObjectMetric();

  void PrintSelf(std::ostream& os, Indent indent) const;

  /* Verify that virtual domain and displacement field are the same size
   * and in the same physical space. */
  virtual void VerifyDisplacementFieldSizeAndPhysicalSpace();

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
  void PreWarpFixedImage( void );
  void PreWarpMovingImage( void );

  /** Flag to track if threading memory has been initialized since last
   * call to Initialize. */
  bool                    m_ThreadingMemoryHasBeenInitialized;

  /* The number of threads to use.
   * Keep private to force use of SetNumberOfThreads in derived classes.
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
