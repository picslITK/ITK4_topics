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
#ifndef __itkTernaryMagnitudeSquaredImageFilter_h
#define __itkTernaryMagnitudeSquaredImageFilter_h

#include "itkTernaryFunctorImageFilter.h"

namespace itk
{
/** \class TernaryMagnitudeSquaredImageFilter
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
class ModulusSquare3
{
public:
  ModulusSquare3() {}
  ~ModulusSquare3() {}
  bool operator!=(const ModulusSquare3 &) const
  {
    return false;
  }

  bool operator==(const ModulusSquare3 & other) const
  {
    return !( *this != other );
  }

  inline TOutput operator()(const TInput1 & A,
                            const TInput2 & B,
                            const TInput3 & C) const
  { return (TOutput)( A * A + B * B + C * C ); }
};
}

template< class TInputImage1, class TInputImage2,
          class TInputImage3, class TOutputImage >
class ITK_EXPORT TernaryMagnitudeSquaredImageFilter:
  public
  TernaryFunctorImageFilter< TInputImage1, TInputImage2,
                             TInputImage3, TOutputImage,
                             Function::ModulusSquare3<
                               typename TInputImage1::PixelType,
                               typename TInputImage2::PixelType,
                               typename TInputImage3::PixelType,
                               typename TOutputImage::PixelType >   >
{
public:
  /** Standard class typedefs. */
  typedef TernaryMagnitudeSquaredImageFilter Self;
  typedef TernaryFunctorImageFilter<
    TInputImage1, TInputImage2,
    TInputImage3, TOutputImage,
    Function::ModulusSquare3<
      typename TInputImage1::PixelType,
      typename TInputImage2::PixelType,
      typename TInputImage3::PixelType,
      typename TOutputImage::PixelType >   >  Superclass;

  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Runtime information support. */
  itkTypeMacro(TernaryMagnitudeSquaredImageFilter,
               TernaryFunctorImageFilter);
protected:
  TernaryMagnitudeSquaredImageFilter() {}
  virtual ~TernaryMagnitudeSquaredImageFilter() {}
private:
  TernaryMagnitudeSquaredImageFilter(const Self &); //purposely not implemented
  void operator=(const Self &);                     //purposely not implemented
};
} // end namespace itk

#endif
