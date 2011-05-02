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
#ifndef __itkCompose2DCovariantVectorImageFilter_h
#define __itkCompose2DCovariantVectorImageFilter_h

#include "itkBinaryFunctorImageFilter.h"
#include "itkCovariantVector.h"

namespace itk
{
namespace Function
{
template< class TInput >
class Compose2DCovariantVector
{
public:
  typedef CovariantVector< TInput, 2 > OutputType;
  Compose2DCovariantVector() {}
  ~Compose2DCovariantVector() {}
  bool operator!=(const Compose2DCovariantVector &) const
  {
    return false;
  }

  bool operator==(const Compose2DCovariantVector & other) const
  {
    return !( *this != other );
  }

  inline OutputType operator()(const TInput & s1,
                               const TInput & s2) const
  {
    OutputType v;

    v[0] = s1;
    v[1] = s2;
    return v;
  }
};
}

/** \class Compose2DCovariantVectorImageFilter
 * \brief Implements pixel-wise composition of an 2D covariant vector pixel from two scalar images.
 *
 * This filter receives two scalar images as input. Each image containing
 * one of the 2D covariant vector components. The filter produces as output a
 * 2D covariant vector image in which the two components have been unified. The Component
 * type is preserved from the PixelType of the input images.
 *
 * \ingroup IntensityImageFilters
 * \ingroup ITK-ImageCompose
 */

template< typename TInputImage,
          typename TOutputImage =
            Image< CovariantVector< ITK_TYPENAME TInputImage::PixelType, 2 >,
                   ::itk::GetImageDimension< TInputImage >::ImageDimension > >
class ITK_EXPORT Compose2DCovariantVectorImageFilter:
  public
  BinaryFunctorImageFilter< TInputImage, TInputImage,
                            TOutputImage,
                            Function::Compose2DCovariantVector< ITK_TYPENAME TInputImage::PixelType >   >
{
public:
  /** Standard class typedefs. */
  typedef Compose2DCovariantVectorImageFilter Self;
  typedef BinaryFunctorImageFilter<
    TInputImage, TInputImage,
    TOutputImage,
    Function::Compose2DCovariantVector< ITK_TYPENAME TInputImage::PixelType > > Superclass;

  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer;

  typedef typename Superclass::OutputImageType OutputImageType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Runtime information support. */
  itkTypeMacro(Compose2DCovariantVectorImageFilter,
               BinaryFunctorImageFilter);

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro( InputHasNumericTraitsCheck,
                   ( Concept::HasNumericTraits< typename TInputImage::PixelType > ) );
  /** End concept checking */
#endif
protected:
  Compose2DCovariantVectorImageFilter() {}
  virtual ~Compose2DCovariantVectorImageFilter() {}
private:
  Compose2DCovariantVectorImageFilter(const Self &); //purposely not implemented
  void operator=(const Self &);                      //purposely not implemented
};
} // end namespace itk

#endif
