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
#ifndef __itkDivideImageFilter_h
#define __itkDivideImageFilter_h

#include "itkBinaryFunctorImageFilter.h"
#include "itkNumericTraits.h"

namespace itk
{
/** \class DivideImageFilter
 * \brief Implements an operator for pixel-wise division of two images.
 *
 * This class is parametrized over the types of the two
 * input images and the type of the output image. When the divisor is zero,
 * the division result is set to the maximum number that can be represneted  by default to
 * avoid exception. Numeric conversions (castings) are done by the C++ defaults.
 *
 * \ingroup IntensityImageFilters  Multithreaded
 * \ingroup ITK-ImageIntensity
 */

namespace Function
{
template< class TInput1, class TInput2, class TOutput >
class Div
{
public:
  Div() {}
  ~Div() {}
  bool operator!=(const Div &) const
  {
    return false;
  }

  bool operator==(const Div & other) const
  {
    return !( *this != other );
  }

  inline TOutput operator()(const TInput1 & A, const TInput2 & B) const
  {
    if ( B != (TInput2)0 )
      {
      return (TOutput)( A / B );
      }
    else
      {
      return NumericTraits< TOutput >::max(A);
      }
  }
};
}

template< class TInputImage1, class TInputImage2, class TOutputImage >
class ITK_EXPORT DivideImageFilter:
  public
  BinaryFunctorImageFilter< TInputImage1, TInputImage2, TOutputImage,
                            Function::Div<
                              typename TInputImage1::PixelType,
                              typename TInputImage2::PixelType,
                              typename TOutputImage::PixelType >   >
{
public:
  /**
   * Standard "Self" typedef.
   */
  typedef DivideImageFilter Self;

  /**
   * Standard "Superclass" typedef.
   */
  typedef BinaryFunctorImageFilter< TInputImage1, TInputImage2, TOutputImage,
                                    Function::Div<
                                      typename TInputImage1::PixelType,
                                      typename TInputImage2::PixelType,
                                      typename TOutputImage::PixelType >
                                    > Superclass;

  /**
   * Smart pointer typedef support
   */
  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer;

  /**
   * Method for creation through the object factory.
   */
  itkNewMacro(Self);

  /** Runtime information support. */
  itkTypeMacro(DivideImageFilter,
               BinaryFunctorImageFilter);

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro( IntConvertibleToInput2Check,
                   ( Concept::Convertible< int, typename TInputImage2::PixelType > ) );
  itkConceptMacro( Input1Input2OutputDivisionOperatorsCheck,
                   ( Concept::DivisionOperators< typename TInputImage1::PixelType,
                                                 typename TInputImage2::PixelType,
                                                 typename TOutputImage::PixelType > ) );
  /** End concept checking */
#endif
protected:
  DivideImageFilter() {}
  virtual ~DivideImageFilter() {}
  DivideImageFilter(const Self &) {}
  void operator=(const Self &) {}
};
} // end namespace itk

#endif
