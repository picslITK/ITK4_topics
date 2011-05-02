/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef __itkVectorLinearInterpolateImageFunction_h
#define __itkVectorLinearInterpolateImageFunction_h

#include "itkVectorInterpolateImageFunction.h"

namespace itk
{
/**
 * \class VectorLinearInterpolateImageFunction
 * \brief Linearly interpolate a vector image at specified positions.
 *
 * VectorLinearInterpolateImageFunction linearly interpolates a vector
 * image intensity non-integer pixel position. This class is templated
 * over the input image type and the coordinate representation type.
 *
 * This function works for N-dimensional images.
 *
 * \warning This function work only for Vector images. For
 * scalar images use LinearInterpolateImageFunction.
 *
 * \ingroup ImageFunctions ImageInterpolators
 *
 * \ingroup ITK-ImageFunction
 */
template< class TInputImage, class TCoordRep = double >
class ITK_EXPORT VectorLinearInterpolateImageFunction:
  public VectorInterpolateImageFunction< TInputImage, TCoordRep >
{
public:
  /** Standard class typedefs. */
  typedef VectorLinearInterpolateImageFunction                     Self;
  typedef VectorInterpolateImageFunction< TInputImage, TCoordRep > Superclass;
  typedef SmartPointer< Self >                                     Pointer;
  typedef SmartPointer< const Self >                               ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(VectorLinearInterpolateImageFunction,
               VectorInterpolateImageFunction);

  /** InputImageType typedef support. */
  typedef typename Superclass::InputImageType InputImageType;
  typedef typename Superclass::PixelType      PixelType;
  typedef typename Superclass::ValueType      ValueType;
  typedef typename Superclass::RealType       RealType;

  /** Grab the vector dimension from the superclass. */
  itkStaticConstMacro(Dimension, unsigned int,
                      Superclass::Dimension);

  /** Dimension underlying input image. */
  itkStaticConstMacro(ImageDimension, unsigned int, Superclass::ImageDimension);

  /** Index typedef support. */
  typedef typename Superclass::IndexType      IndexType;

  /** ContinuousIndex typedef support. */
  typedef typename Superclass::ContinuousIndexType ContinuousIndexType;

  /** Output type is Vector<double,Dimension> */
  typedef typename Superclass::OutputType OutputType;

  /** Evaluate the function at a ContinuousIndex position
   *
   * Returns the linearly interpolated image intensity at a
   * specified point position. No bounds checking is done.
   * The point is assume to lie within the image buffer.
   *
   * ImageFunction::IsInsideBuffer() can be used to check bounds before
   * calling the method. */
  virtual OutputType EvaluateAtContinuousIndex(
    const ContinuousIndexType & index) const;

protected:
  VectorLinearInterpolateImageFunction();
  ~VectorLinearInterpolateImageFunction(){}
  void PrintSelf(std::ostream & os, Indent indent) const;

private:
  VectorLinearInterpolateImageFunction(const Self &); //purposely not
                                                      // implemented
  void operator=(const Self &);                       //purposely not
                                                      // implemented

  /** Number of neighbors used in the interpolation */
  static const unsigned long m_Neighbors;
};
} // end namespace itk

// Define instantiation macro for this template.
#define ITK_TEMPLATE_VectorLinearInterpolateImageFunction(_, EXPORT, TypeX, TypeY)     \
  namespace itk                                                                        \
  {                                                                                    \
  _( 2 ( class EXPORT VectorLinearInterpolateImageFunction< ITK_TEMPLATE_2 TypeX > ) ) \
  namespace Templates                                                                  \
  {                                                                                    \
  typedef VectorLinearInterpolateImageFunction< ITK_TEMPLATE_2 TypeX >                 \
  VectorLinearInterpolateImageFunction##TypeY;                                       \
  }                                                                                    \
  }

#if ITK_TEMPLATE_EXPLICIT
#include "Templates/itkVectorLinearInterpolateImageFunction+-.h"
#endif

#if ITK_TEMPLATE_TXX
#include "itkVectorLinearInterpolateImageFunction.txx"
#endif

#endif
