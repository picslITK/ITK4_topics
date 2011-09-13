/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkDemonsImageToImageObjectMetric.h,v $
  Language:  C++
  Date:      $Date: $
  Version:   $Revision: $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkDemonsImageToImageObjectMetric_h
#define __itkDemonsImageToImageObjectMetric_h

#include "itkImageToImageObjectMetric.h"

namespace itk
{

/** \class DemonsImageToImageObjectMetric
 *
 *  \brief Class implementing rudimentary demons metric.
 *
 *  See \c GetValueAndDerivativeProcessPoint for algorithm implementation.
 *
 * \ingroup ITKHighDimensionalMetric
 */
template <class TFixedImage,
          class TMovingImage,
          class TVirtualImage = TFixedImage >
class ITK_EXPORT DemonsImageToImageObjectMetric :
public ImageToImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage>
{
public:

  /** Standard class typedefs. */
  typedef DemonsImageToImageObjectMetric                      Self;
  typedef ImageToImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage>
                                                              Superclass;
  typedef SmartPointer<Self>                                  Pointer;
  typedef SmartPointer<const Self>                            ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(DemonsImageToImageObjectMetric, ImageToImageObjectMetric);

  /** superclass types */
  typedef typename Superclass::MeasureType             MeasureType;
  typedef typename Superclass::DerivativeType          DerivativeType;
  typedef typename Superclass::VirtualPointType        VirtualPointType;
  typedef typename Superclass::FixedImagePointType     FixedImagePointType;
  typedef typename Superclass::FixedImagePixelType     FixedImagePixelType;
  typedef typename Superclass::FixedImageGradientType
                                                     FixedImageGradientType;

  typedef typename Superclass::MovingImagePointType    MovingImagePointType;
  typedef typename Superclass::MovingImagePixelType    MovingImagePixelType;
  typedef typename Superclass::MovingImageGradientType
                                                    MovingImageGradientType;

  typedef typename Superclass::MovingTransformType     MovingTransformType;
  typedef typename Superclass::JacobianType            JacobianType;

  /** Initialize. Must be called before first call to GetValue or
   *  GetValueAndDerivative, after metric settings are changed. */
  virtual void Initialize(void) throw ( itk::ExceptionObject );

  using Superclass::GetValueAndDerivative;
  void GetValueAndDerivative( MeasureType & value, DerivativeType & derivative);

  /** Evaluate and return the metric value.
   * \warning Not yet implemented. */
  using Superclass::GetValue;
  MeasureType GetValue();

protected:

  /* Worker routine to process each point */
  bool GetValueAndDerivativeProcessPoint(
                    const VirtualPointType &           virtualPoint,
                    const FixedImagePointType &        mappedFixedPoint,
                    const FixedImagePixelType &        fixedImageValue,
                    const FixedImageGradientType &  fixedImageGradient,
                    const MovingImagePointType &       mappedMovingPoint,
                    const MovingImagePixelType &       movingImageValue,
                    const MovingImageGradientType & movingImageGradient,
                    MeasureType &                      metricValueResult,
                    DerivativeType &                   localDerivativeReturn,
                    ThreadIdType                       threadID);

  DemonsImageToImageObjectMetric();
  virtual ~DemonsImageToImageObjectMetric();
  void PrintSelf(std::ostream& os, Indent indent) const;

private:

  //purposely not implemented
  DemonsImageToImageObjectMetric(const Self &);
  //purposely not implemented
  void operator=(const Self &);

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkDemonsImageToImageObjectMetric.hxx"
#endif

#endif
