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
#ifndef __itkTernaryMagnitudeImageFilter_h
#define __itkTernaryMagnitudeImageFilter_h

#include "itkTernaryFunctorImageFilter.h"

namespace itk
{
/** \class TernaryMagnitudeImageFilter
 * \brief Implements pixel-wise addition of three images.
 *
 * This class is parametrized over the types of the three
 * input images and the type of the output image.
 * Numeric conversions (castings) are done by the C++ defaults.
 *
 * \ingroup IntensityImageFilters
 * \ingroup ITK-ImageIntensity
 */
namespace Function
{
template< class TInput1, class TInput2, class TInput3, class TOutput >
class Modulus3
{
public:
  Modulus3() {}
  ~Modulus3() {}
  bool operator!=(const Modulus3 &) const
  {
    return false;
  }

  bool operator==(const Modulus3 & other) const
  {
    return !( *this != other );
  }

  inline TOutput operator()(const TInput1 & A,
                            const TInput2 & B,
                            const TInput3 & C) const
  { return (TOutput)vcl_sqrt( (double)( A * A + B * B + C * C ) ); }
};
}

template< class TInputImage1, class TInputImage2,
          class TInputImage3, class TOutputImage >
class ITK_EXPORT TernaryMagnitudeImageFilter:
  public
  TernaryFunctorImageFilter< TInputImage1, TInputImage2,
                             TInputImage3, TOutputImage,
                             Function::Modulus3<
                               typename TInputImage1::PixelType,
                               typename TInputImage2::PixelType,
                               typename TInputImage3::PixelType,
                               typename TOutputImage::PixelType >   >
{
public:
  /** Standard class typedefs. */
  typedef TernaryMagnitudeImageFilter Self;
  typedef TernaryFunctorImageFilter<
    TInputImage1, TInputImage2,
    TInputImage3, TOutputImage,
    Function::Modulus3<
      typename TInputImage1::PixelType,
      typename TInputImage2::PixelType,
      typename TInputImage3::PixelType,
      typename TOutputImage::PixelType >
    >                                   Superclass;

  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Runtime information support. */
  itkTypeMacro(TernaryMagnitudeImageFilter,
               TernaryFunctorImageFilter);
protected:
  TernaryMagnitudeImageFilter() {}
  virtual ~TernaryMagnitudeImageFilter() {}
private:
  TernaryMagnitudeImageFilter(const Self &); //purposely not implemented
  void operator=(const Self &);              //purposely not implemented
};
} // end namespace itk

#endif
