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
#include "itkObjectToObjectMetric.h"

namespace itk
{
/** \class ImageToImageObjectMetric
 *
 * Templated over the fixed and moving image types, as well as an optional
 * VirtualImage type to define the virtual domain. The VirtualImage type
 * defaults to TFixedImage.
 *
 */
template <class TFixedImage,class TMovingImage, class TVirtualImage = TFixedImage >
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
  typedef TMovingImage                            MovingImageType;
  typedef typename MovingImageType::PixelType     MovingImagePixelType;
  typedef typename MovingImageType::Pointer       MovingImagePointer;

  /** Types for the virtual domain */
  typedef TVirtualImage                             VirtualImageType;
  typedef typename VirtualImageType::RegionType     VirtualRegionType;
  typedef typename RegionType::SizeType             VirtualSizeType;
  typedef typename VirtualImageType::SpacingType    VirtualSpacingType;
  typedef typename VirtualImageType::PointType      VirtualOriginType;
  typedef typename VirtualImageType::PointType      VirtualPointType;
  typedef typename VirtualImageType::DirectionType  VirtualDirectionType;
  typedef typename VirtualImageType::SizeType       VirtualRadiusType;
  typedef typename VirtualImageType::IndexType      VirtualIndexType;
  typedef typename VirtualImageType::RegionType     VirtualImageRegionType;

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

  /**  Type of the Interpolator Base class */
  should this be pointer to base interpolator, and add accessors for user to set? Default
  to linear?
  typedef LinearInterpolateImageFunction< FixedImageType,
                                          CoordinateRepresentationType >
                                                      FixedInterpolatorType;
  typedef LinearInterpolateImageFunction< MovingImageType,
                                          CoordinateRepresentationType >
                                                      MovingInterpolatorType;
  typedef typename FixedInterpolatorType::Pointer     FixedInterpolatorPointer;
  typedef typename MovingInterpolatorType::Pointer    MovingInterpolatorPointer;

  /** Image derivatives types */
  typedef   CovariantVector< CoordinateRepresentationType,
                             itkGetStaticConstMacro(FixedImageDimension) >
                                                      FixedImageDerivativesType;
  typedef   CovariantVector< CoordinateRepresentationType,
                             itkGetStaticConstMacro(MovingImageDimension) >
                                                      MovingImageDerivativesType;

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

  /** Set/get images */
  should these be Object macros instead?
  itkSetMacro(FixedImage, FixedImagePointer);
  itkGetMacro(FixedImage, FixedImagePointer);
  itkSetMacro(MovingImage, MovingImagePointer);
  itkGetMacro(MovingImage, MovingImagePointer);

  /** Define the virtual reference space. This space defines the resolution,
   *  at which the registration is performed as well as the physical coordinate
   *  system.  Useful for unbiased registration.  If the user does not set this
   *  explicitly then it is taken from the fixed image in Initialize method.
   *  To set all settings from a reference image, use SetVirtualDomainImage.
   */
  itkSetMacro(VirtualDomainSpacing, VirtualSpacingType);
  itkGetMacro(VirtualDomainSpacing, VirtualSpacingType);
  itkSetMacro(VirtualDomainSize, VirtualSizeType);
  itkGetMacro(VirtualDomainSize, VirtualSizeType);
  itkSetMacro(VirtualDomainIndex, VirtualIndexType);
  itkGetMacro(VirtualDomainIndex, VirtualIndexType);
  itkSetMacro(VirtualDomainOrigin, VirtualOriginType);
  itkGetMacro(VirtualDomainOrigin, VirtualOriginType);
  itkSetMacro(VirtualDomainDirection, VirtualDirectionType);
  itkGetMacro(VirtualDomainDirection, VirtualDirectionType);

  /** Set/Get all virtual domain setings at once via an image */
  virtual void SetVirtualDomainImage(VirtualImageType* image);
  itkGetMacro(VirtualDomainImage, VirtualImageType);

  /** Connect the fixed transform. */
  itkSetObjectMacro(FixedImageTransform, FixedTransformType);
  /** Get a pointer to the fixed transform.  */
  itkGetConstObjectMacro(FixedImageTransform, FixedTransformType);
  /** Connect the moving transform. */
  itkSetObjectMacro(MovingImageTransform, MovingTransformType);
  /** Get a pointer to the moving transform.  */
  itkGetConstObjectMacro(MovingImageTransform, MovingTransformType);

  /** Return the number of parameters,
   * taking CoordinateSystemType into account */
  virtual unsigned int GetNumberOfParameters() const;

  /* Initialize the metric before calling GetValue or GetDerivative.
   * Derived classes must call this Superclass version if they override
   * this to perform their own initialization */
  virtual void Initialize(void) throw ( itk::ExceptionObject );

protected:

  ImageToImageObjectMetric();
  virtual ~ImageToImageObjectMetric();
  void PrintSelf(std::ostream& os, Indent indent) const;

private:

  //purposely not implemented
  ImageToImageObjectMetric(const Self &);
  //purposely not implemented
  void operator=(const Self &);

to sort out:
  FixedImagePointer m_FixedImage;
  TransformPointer m_FixedImageTransform;
  MovingImagePointer m_MovingImage;
  TransformPointer m_MovingImageTransform;
  FixedImagePointer m_VirtualImage;
  MeasureType*     m_ThreaderMSE;
  DerivativeType*  m_ThreaderDerivatives;
  double           m_Normalizer;

  unsigned int m_NumberOfThreads;
  unsigned int m_NumberOfParameters;

  IndexType  m_VirtualDomainIndex;
  SizeType  m_VirtualDomainSize;
  SpacingType  m_VirtualDomainSpacing;
  OriginType  m_VirtualDomainOrigin;
  DirectionType  m_VirtualDomainDirection;

  FixedInterpolatorPointer m_FixedInterpolator;
  MovingInterpolatorPointer m_MovingInterpolator;

  unsigned int m_InputImageVectorLength;
};
}//namespace itk
#endif
