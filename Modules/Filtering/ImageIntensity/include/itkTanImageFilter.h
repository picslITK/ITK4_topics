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
#ifndef __itkTanImageFilter_h
#define __itkTanImageFilter_h

#include "itkUnaryFunctorImageFilter.h"
#include "vnl/vnl_math.h"

namespace itk
{
/** \class TanImageFilter
 * \brief Computes the vcl_tan(x) pixel-wise
 *
 * \ingroup IntensityImageFilters  Multithreaded
 * \ingroup ITK-ImageIntensity
 */
namespace Function
{
template< class TInput, class TOutput >
class Tan
{
public:
  Tan() {}
  ~Tan() {}
  bool operator!=(const Tan &) const
  {
    return false;
  }

  bool operator==(const Tan & other) const
  {
    return !( *this != other );
  }

  inline TOutput operator()(const TInput & A) const
  { return (TOutput)vcl_tan( (double)A ); }
};
}
template< class TInputImage, class TOutputImage >
class ITK_EXPORT TanImageFilter:
  public
  UnaryFunctorImageFilter< TInputImage, TOutputImage,
                           Function::Tan< typename TInputImage::PixelType,
                                          typename TOutputImage::PixelType >   >
{
public:
  /** Standard class typedefs. */
  typedef TanImageFilter Self;
  typedef UnaryFunctorImageFilter<
    TInputImage, TOutputImage,
    Function::Tan< typename TInputImage::PixelType,
                   typename TOutputImage::PixelType > >  Superclass;

  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Runtime information support. */
  itkTypeMacro(TanImageFilter,
               UnaryFunctorImageFilter);

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro( InputConvertibleToDoubleCheck,
                   ( Concept::Convertible< typename TInputImage::PixelType, double > ) );
  itkConceptMacro( DoubleConvertibleToOutputCheck,
                   ( Concept::Convertible< double, typename TOutputImage::PixelType > ) );
  /** End concept checking */
#endif
protected:
  TanImageFilter() {}
  virtual ~TanImageFilter() {}
private:
  TanImageFilter(const Self &); //purposely not implemented
  void operator=(const Self &); //purposely not implemented
};
} // end namespace itk

#endif
