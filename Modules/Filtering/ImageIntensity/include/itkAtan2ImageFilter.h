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
#ifndef __itkAtan2ImageFilter_h
#define __itkAtan2ImageFilter_h

#include "itkBinaryFunctorImageFilter.h"
#include "vnl/vnl_math.h"

namespace itk
{
/** \class Atan2ImageFilter
 * \brief Computes arctangent pixel-wise from two images.
 *
 * This class is parametrized over the types of the two
 * input images and the type of the output image.
 * Numeric conversions (castings) are done by the C++ defaults.
 *
 * Both pixel input types are casted to \c double in order to be
 * used as parameters of \c vcl_atan2(). The resulting \c double value
 * is casted to the output pixel type.
 *
 * \ingroup IntensityImageFilters Multithreaded
 * \ingroup ITK-ImageIntensity
 */
namespace Functor
{
template< class TInput1, class TInput2, class TOutput >
class Atan2
{
public:
  Atan2() {}
  ~Atan2() {}
  bool operator!=(const Atan2 &) const
  {
    return false;
  }

  bool operator==(const Atan2 & other) const
  {
    return !( *this != other );
  }

  inline TOutput operator()(const TInput1 & A, const TInput2 & B) const
  {
    return static_cast< TOutput >(
             vcl_atan2(
               static_cast< double >( A ),
               static_cast< double >( B ) )
             );
  }
};
}

template< class TInputImage1, class TInputImage2, class TOutputImage >
class ITK_EXPORT Atan2ImageFilter:
  public
  BinaryFunctorImageFilter< TInputImage1, TInputImage2, TOutputImage,
                            Functor::Atan2<
                              typename TInputImage1::PixelType,
                              typename TInputImage2::PixelType,
                              typename TOutputImage::PixelType >   >
{
public:
  /** Standard class typedefs. */
  typedef Atan2ImageFilter Self;
  typedef BinaryFunctorImageFilter< TInputImage1, TInputImage2, TOutputImage,
                                    Functor::Atan2<
                                      typename TInputImage1::PixelType,
                                      typename TInputImage2::PixelType,
                                      typename TOutputImage::PixelType >
                                    >                                 Superclass;

  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Runtime information support. */
  itkTypeMacro(Atan2ImageFilter,
               BinaryFunctorImageFilter);

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro( Input1ConvertibleToDoubleCheck,
                   ( Concept::Convertible< typename TInputImage1::PixelType, double > ) );
  itkConceptMacro( Input2ConvertibleToDoubleCheck,
                   ( Concept::Convertible< typename TInputImage2::PixelType, double > ) );
  itkConceptMacro( DoubleConvertibleToOutputCheck,
                   ( Concept::Convertible< double, typename TOutputImage::PixelType > ) );
  /** End concept checking */
#endif
protected:
  Atan2ImageFilter() {}
  virtual ~Atan2ImageFilter() {}
private:
  Atan2ImageFilter(const Self &); //purposely not implemented
  void operator=(const Self &);   //purposely not implemented
};
} // end namespace itk

#endif
