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
 * 2) otherwise, GradientRecursiveGaussianImageFilter is used by default.
 * This requires memory to hold the fully computed gradient image.
 * 3) lastly, CentralDifferenceImageFunction is used if m_UseGaussianGradient
 * option has been set to false.
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
  itkTypeMacro(DemonsImageToImageMetric, ObjectToObjectMetric);

  /** Type used internally for computations */
  typedef double  InternalComputationValueType;

  /** Type used for representing parameter values  */
  typedef typename Superclass::CoordinateRepresentationType
                                                  CoordinateRepresentationType;

  /**  Type of the parameters. */
  typedef typename Superclass::ParametersType       ParametersType;
  typedef typename Superclass::ParametersValueType  ParametersValueType;

  /** Image-accessor typedefs */
  typedef TFixedImage                             FixedImageType;
  typedef typename FixedImageType::PixelType      FixedImagePixelType;
  typedef typename FixedImageType::Pointer        FixedImagePointer;
  typedef typename FixedImageType::ConstPointer   FixedImageConstPointer;
  typedef TMovingImage                            MovingImageType;
  typedef typename MovingImageType::PixelType     MovingImagePixelType;
  typedef typename MovingImageType::Pointer       MovingImagePointer;
  typedef typename MovingImageType::ConstPointer  MovingImageConstPointer;

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
    itkGetStaticConstMacro( MovingPointSetDimension ),
    itkGetStaticConstMacro( FixedPointSetDimension )>  MovingTransformType;

  typedef Transform<CoordinateRepresentationType,
    itkGetStaticConstMacro( FixedPointSetDimension ),
    itkGetStaticConstMacro( MovingPointSetDimension )> FixedTransformType;

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

  /** Gradient calculator types. */
  typedef CentralDifferenceImageFunction<
                                FixedImageType,
                                CoordinateRepresentationType>
                                               FixedGradientCalculatorType;
  typedef typename FixedImageGradientCalculatorType::Pointer
                                               FixedGradientCalculatorPointer;
  typedef CentralDifferenceImageFunction<
                                MovingImageType,
                                CoordinateRepresentationType>
                                               MovingGradientCalculatorType;
  typedef typename MovingImageGradientCalculatorType::Pointer
                                               MovingGradientCalculatorPointer;

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
   *  explicitly then it is taken from the fixed image in Initialize method.
   *  To set all settings from a reference image, use SetVirtualDomainImage.
   */
  itkSetMacro(VirtualDomainSpacing, VirtualSpacingType);
  itkGetMacro(VirtualDomainSpacing, VirtualSpacingType);
  itkSetMacro(VirtualDomainOrigin, VirtualOriginType);
  itkGetMacro(VirtualDomainOrigin, VirtualOriginType);
  itkSetMacro(VirtualDomainDirection, VirtualDirectionType);
  itkGetMacro(VirtualDomainDirection, VirtualDirectionType);

  /** Set/Get the size of the virtual domain region over which to compute
   * the metric. Any subsequent call to SetVirtualRegion will override
   * these values. */
  void SetVirtualDomainSize( VirtualDomainIndex& index );
  const VirtualDomainSize GetVirtualDomainIndex( void );

  /** Set/Get the index of the virtual domain region over which to compute
   * the metric. Any subsequent call to SetVirtualRegion will override
   * these values. */
  void SetVirtualDomainIndex( VirtualDomainIndex& index );
  const VirtualDomainIndex GetVirtualDomainIndex( void );

  /** Set/Get the region of the virtual domain over which to compute
   * the metric. The region passed here will be modified by
   * any subsequent calls to SetVirtualDomainSize or SetVirtualDomainIndex */
  itkSetMacro(VirtualDomainRegion, VirtualRegionType);
  itkGetMacro(VirtualDomainRegion, VirtualRegionType);

  /** Set/Get all virtual domain setings at once via an image.
   * The image is expected to have already been allocated. */
  virtual void SetVirtualDomainImage(VirtualImageType* image);
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
  itkSetMacro(UseGaussianGradient, bool);
  itkGetConstReferenceMacro(UseGaussianGradient, bool);
  itkBooleanMacro(UseGaussianGradient);

  /** Computes the gradient image of the fixed and/or moving image(s),
   * depending on CoordinateSystemType setting, using the
   * GradientRecursiveGaussianImageFilter, assigning the output to
   * and assigns it to m_[Fixed|Moving]GaussianGradientImage. */
  virtual void ComputeGaussianGradient(void);

  /** Get Fixed Gradient Image. */
  itkGetConstObjectMacro(FixedGaussianGradientImage, FixedGradientImageType);
  /** Get Moving Gradient Image. */
  itkGetConstObjectMacro(MovingGaussianGradientImage, MovingGradientImageType);

  /** Return the number of parameters,
   * taking CoordinateSystemType into account */
  virtual unsigned int GetNumberOfParameters() const;

  /* Initialize the metric before calling GetValue or GetDerivative.
   * Derived classes must call this Superclass version if they override
   * this to perform their own initialization.
   * \note This is meant to be called once for a particular metric setup.
   * That is, when used in registration, this method would be called once
   * before entering the registration loop, during which GetValue or
   * GetDerivative will be called repeatedly. */
  virtual void Initialize(void) throw ( itk::ExceptionObject );

protected:

  ImageToImageObjectMetric();
  virtual ~ImageToImageObjectMetric();
  void PrintSelf(std::ostream& os, Indent indent) const;

  /** Transform a point from VirtualImage domain to FixedImage domain.
   * This function also checks if mapped point is within the mask if
   * one is set, and that is within the moving image buffer, in which
   * case \c pointIsValid will be true on return.
   * \c mappedFixedPoint and \c fixedImageValue are also returned.
   * \note It would be better for maintainence to have a single method
   * that could work for either fixed or moving domains. However setting
   * that up is complicated because dimensionality and pixel type may
   * be different between the two. */
  virtual void TransformFixedPoint(const VirtualPointType & point,
                              FixedImagePointType & mappedFixedPoint,
                              bool & pointIsValid,
                              FixedImagePixelType & FixedImageValue,
                              unsigned int threadID) const;

  /** Transform a point from VirtualImage domain to MovingImage domain.
   * This function also checks if mapped point is within the mask if
   * one is set, and that is within the moving image buffer. */
  virtual void TransformMovingPoint(const VirtualPointType & point,
                              MovingImagePointType & mappedMovingPoint,
                              bool & pointIsValid,
                              MovingImagePixelType & movingImageValue,
                              unsigned int threadID) const;

  virtual void TransformPointWithDerivatives(
                                            todo
                                              ) const;

  /** Compute image derivatives at a point. */
  virtual void ComputeFixedImageDerivatives(
                                    const FixedImagePointType & mappedPoint,
                                    FixedImageDerivativesType & gradient,
                                    unsigned int threadID) const;
  virtual void ComputeMovingImageDerivatives(
                                    const MovingImagePointType & mappedPoint,
                                    MovingImageDerivativesType & gradient,
                                    unsigned int threadID) const;

  FixedImageConstPointer  m_FixedImage;
  FixedTransformPointer   m_FixedImageTransform;
  MovingImageConstPointer m_MovingImage;
  MovingTransformPointer  m_MovingImageTransform;
  VirtualImagePointer     m_VirtualImage;

  /* needed ? I think we'll just get on the fly from GetNumberOfParameters? */
  unsigned int            m_NumberOfParameters;

  VirtualRegionType       m_VirtualDomainRegion;
  VirtualSpacingType      m_VirtualDomainSpacing;
  VirtualOriginType       m_VirtualDomainOrigin;
  VirtualDirectionType    m_VirtualDomainDirection;

  /** Pointers to interpolator bases */
  FixedInterpolatorPointer  m_FixedInterpolator;
  MovingInterpolatorPointer m_MovingInterpolator;
  /** Pointers to BSplineInterpolators */
  typename FixedBSplineInterpolatorType::Pointer  m_FixedBSplineInterpolator;
  typename MovingBSplineInterpolatorType::Pointer m_MovingBSplineInterpolator;
  /** Boolean to indicate if the fixed interpolator is BSpline. */
  bool                               m_MovingInterpolatorIsBSpline;
  /** Boolean to indicate if the fixed interpolator is BSpline. */
  bool                               m_MovingInterpolatorIsBSpline;

  /** Flag to control use of Gaussian gradient filter or gradient calculator
   * for image gradient calculations. */
  bool                               m_UseGaussianGradient;
  /** Gradient images to store Gaussian gradient filter output. */
  FixedGradientImageType::Pointer    m_FixedGaussianGradientImage;
  MovingGradientImageType::Pointer   m_MovingGaussianGradientImage;

  /* This was from DemonsImageToImage - not sure if still needed. Seems to
   * always be 1 there, but maybe b/c wasn't fully setup for vector images yet? */
  unsigned int m_InputImageVectorLength;

  /** User worker method to calculate value and derivative at a given point.
   * Must be provided by derived classes.
   * \warning  This is called from the threader, and thus must be thread-safe.
   */
  virtual inline bool GetValueAndDerivativeProcessPoint(
            const VirtualPointType &      itkNotUsed(virtualPoint),
            const FixedImagePointType &   itkNotUsed(mappedFixedPoint),
            const FixedImagePixelType &   itkNotUsed(fixedImageValue),
            const ImageDerivativesType &  itkNotUsed(fixedImageGradientValue),
            const MovingImagePointType &  itkNotUsed(mappedMovingPoint),
            const MovingImagePixelType &  itkNotUsed(movingImageValue),
            const ImageDerivativesType &  itkNotUsed(movingImageGradientValue),
            unsigned int                  itkNotUsed(threadID) )
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
  void GetValueAndDerivativeMultiThreadedInitiate();
    {
    maybe have flag in method args to zero the gradientAccumulator var,
      which defaults to true. For use from MultiVariate metric.
    init value sum to zero, and per-thread value sums
    start threader
    }

  /** Default post-processing after multi-threaded calculation of
   * value and derivative. Typically called by derived classes after
   * GetValueAndDerivativeMultiThreadedInitiate.
   * Derived classes may not call this if they require special handling.
   */
  virtual void GetValueAndDerivativeMultiThreadedPostProcess();
    {
    collect values and derivs from threads
    calc avg based on number of points calced
    }

  /** Type of the default threader used for GetValue and GetDerivative.
   * This splits an image region in per-thread sub-regions over the outermost
   * image dimension. */
  typedef ImageToData< VirtualImageDimension, Self >
                                             ValueAndDerivativeThreaderType;
  typedef typename ValueAndDerivativeThreaderType::InputObjectType
                                             ThreaderInputObjectType;
  /** Multi-threader callback used to iterate over image region by thread,
   * and call the derived class' user worker method to calculate
   * value and derivative. */
  //NOTE: this should be private?
  static void GetValueAndDerivativeMultiThreadedCallback(
                                                const ThreaderInputObjectType&,
                                                int threadId,
                                                Self * holder);

  /** Threader used to evaluate value and deriviative. */
  ValueAndDerivativeThreaderType::Pointer    m_ValueAndDerivativeThreader;

  itkSetObjectMacro(ValueAndDerivativeThreader,ValueAndDerivativeThreaderType);

  /** Intermediary threaded metric value storage. */
  std::vector<InternalComputationValueType>   m_MeasurePerThread;
  other threading stuff from GradientDescent and Obj2ObjOptimizer...
  need Set/Get routines too...
  TODO
    - in SetNumberOfThreads or better still during init for threading, must call
      BSplineInterpolator->SetNumberOfThreads so it can do what it
      needs. Old Img2Img does this in MUltiThreadingInitialize

  unsigned int            m_NumberOfThreads;

private:

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
