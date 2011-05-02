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
#ifndef __itkMaximumImageFilter_h
#define __itkMaximumImageFilter_h

#include "itkBinaryFunctorImageFilter.h"

namespace itk
{
/** \class MaximumImageFilter
 * \brief Implements a pixel-wise operator Max(a,b) between two images.
 *
 * The pixel values of the output image are the maximum between the
 * corresponding pixels of the two input images.
 *
 * This class is parametrized over the types of the two
 * input images and the type of the output image.
 * Numeric conversions (castings) are done by the C++ defaults.
 *
 * \ingroup IntensityImageFilters  Multithreaded
 * \ingroup ITK-ImageIntensity
 * \wikiexample{ImageProcessing/MaximumImageFilter,Pixel wise compare two input images and set the output pixel to their max}
 */
namespace Function
{
template< class TInput1, class TInput2 = TInput1, class TOutput = TInput1 >
class Maximum
{
public:
  Maximum() {}
  ~Maximum() {}
  bool operator!=(const Maximum &) const
  {
    return false;
  }

  bool operator==(const Maximum & other) const
  {
    return !( *this != other );
  }

  inline TOutput operator()(const TInput1 & A, const TInput2 & B) const
  {
    if ( A > B )
      {
      return static_cast< TOutput >( A );
      }
    else
      {
      return static_cast< TOutput >( B );
      }
  }
};
}

template< class TInputImage1, class TInputImage2 = TInputImage1, class TOutputImage = TInputImage1 >
class ITK_EXPORT MaximumImageFilter:
  public
  BinaryFunctorImageFilter< TInputImage1, TInputImage2, TOutputImage,
                            Function::Maximum<
                              typename TInputImage1::PixelType,
                              typename TInputImage2::PixelType,
                              typename TOutputImage::PixelType >   >
{
public:
  /** Standard class typedefs. */
  typedef MaximumImageFilter Self;
  typedef BinaryFunctorImageFilter< TInputImage1, TInputImage2, TOutputImage,
                                    Function::Maximum<
                                      typename TInputImage1::PixelType,
                                      typename TInputImage2::PixelType,
                                      typename TOutputImage::PixelType >
                                    >                                 Superclass;

  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Runtime information support. */
  itkTypeMacro(MaximumImageFilter,
               BinaryFunctorImageFilter);

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro( Input1ConvertibleToOutputCheck,
                   ( Concept::Convertible< typename TInputImage1::PixelType,
                                           typename TOutputImage::PixelType > ) );
  itkConceptMacro( Input2ConvertibleToOutputCheck,
                   ( Concept::Convertible< typename TInputImage2::PixelType,
                                           typename TOutputImage::PixelType > ) );
  itkConceptMacro( Input1GreaterThanInput2Check,
                   ( Concept::GreaterThanComparable< typename TInputImage1::PixelType,
                                                     typename TInputImage2::PixelType > ) );
  /** End concept checking */
#endif
protected:
  MaximumImageFilter() {}
  virtual ~MaximumImageFilter() {}
private:
  MaximumImageFilter(const Self &); //purposely not implemented
  void operator=(const Self &);     //purposely not implemented
};
} // end namespace itk

#endif
