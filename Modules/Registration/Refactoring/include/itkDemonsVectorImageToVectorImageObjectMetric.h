/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkDemonsVectorImageToVectorImageObjectMetric.h,v $
  Language:  C++
  Date:      $Date: $
  Version:   $Revision: $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkDemonsVectorImageToVectorImageObjectMetric_h
#define __itkDemonsVectorImageToVectorImageObjectMetric_h

#include "itkVectorImageToVectorImageObjectMetric.h"

namespace itk
{

template <class TFixedImage,
          class TMovingImage,
          class TVirtualImage = TFixedImage >
class ITK_EXPORT DemonsVectorImageToVectorImageObjectMetric :
public VectorImageToVectorImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage>
{
public:

  /** Standard class typedefs. */
  typedef DemonsVectorImageToVectorImageObjectMetric                      Self;
  typedef VectorImageToVectorImageObjectMetric<TFixedImage, TMovingImage, TVirtualImage>
                                                              Superclass;
  typedef SmartPointer<Self>                                  Pointer;
  typedef SmartPointer<const Self>                            ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(DemonsVectorImageToVectorImageObjectMetric, VectorImageToVectorImageObjectMetric);

  /** superclass types */
  typedef typename Superclass::MeasureType             MeasureType;
  typedef typename Superclass::DerivativeType          DerivativeType;
  typedef typename Superclass::VirtualPointType        VirtualPointType;
  typedef typename Superclass::FixedImagePointType     FixedImagePointType;
  typedef typename Superclass::FixedImagePixelType     FixedImagePixelType;
  typedef typename Superclass::FixedImageDerivativesType
                                                       FixedImageDerivativesType;
  typedef typename Superclass::MovingImagePointType    MovingImagePointType;
  typedef typename Superclass::MovingImagePixelType    MovingImagePixelType;
  typedef typename Superclass::MovingImageDerivativesType
                                                       MovingImageDerivativesType;
  typedef typename Superclass::MovingTransformType     MovingTransformType;
  typedef typename MovingTransformType::JacobianType
                                                       MovingTransformJacobianType;

  /** Initialize. Must be called before first call to GetValue or
   *  GetValueAndDerivative, after metric settings are changed. */
  virtual void Initialize(void) throw ( itk::ExceptionObject );

  /** Evaluate and return the value and derivative */
  void GetValueAndDerivative( MeasureType & value, DerivativeType & derivative);

  /** Evaluate and return the metric value */
  MeasureType GetValue()
  { itkExceptionMacro("GetValue not yet implemented."); }

protected:

  /* Worker routine to process each point */
  bool GetValueAndDerivativeProcessPoint(
                    const VirtualPointType &           virtualPoint,
                    const FixedImagePointType &        mappedFixedPoint,
                    const FixedImagePixelType &        fixedImageValue,
                    const FixedImageDerivativesType &  fixedImageDerivatives,
                    const MovingImagePointType &       mappedMovingPoint,
                    const MovingImagePixelType &       movingImageValue,
                    const MovingImageDerivativesType & movingImageDerivatives,
                    MeasureType &                      metricValueResult,
                    DerivativeType &                   localDerivativeReturn,
                    ThreadIdType                       threadID);

  DemonsVectorImageToVectorImageObjectMetric();
  virtual ~DemonsVectorImageToVectorImageObjectMetric();
  void PrintSelf(std::ostream& os, Indent indent) const;

private:

  //purposely not implemented
  DemonsVectorImageToVectorImageObjectMetric(const Self &);
  //purposely not implemented
  void operator=(const Self &);

  MovingTransformJacobianType m_Jacobian;

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkDemonsVectorImageToVectorImageObjectMetric.hxx"
#endif

#endif
