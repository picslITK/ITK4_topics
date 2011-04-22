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
#ifndef __itkValuedRegionalMinimaImageFilter_h
#define __itkValuedRegionalMinimaImageFilter_h

#include "itkValuedRegionalExtremaImageFilter.h"
#include "itkConceptChecking.h"

namespace itk
{
/** \class ValuedRegionalMinimaImageFilter
 * \brief Transforms the image so that any pixel that is not a
 * regional minima is set to the maximum value for the pixel
 * type. Pixels that are regional minima retain their value.
 *
 * Regional minima are flat zones surrounded by pixels of higher
 * value. A completely flat image will be marked as a regional minima
 * by this filter.

 * \author Richard Beare. Department of Medicine, Monash University,
 * Melbourne, Australia.
 *
 * \sa ValuedRegionalMaximaImageFilter, ValuedRegionalExtremaImageFilter,
 * \sa HMinimaImageFilter
 * \ingroup MathematicalMorphologyImageFilters
 * \ingroup ITK-Review
 */

template< class TInputImage, class TOutputImage >
class ITK_EXPORT ValuedRegionalMinimaImageFilter:
  public
  ValuedRegionalExtremaImageFilter< TInputImage, TOutputImage,
                                    std::less< typename TInputImage::PixelType >,
                                    std::less< typename TOutputImage::PixelType >
                                    >
{
public:
  typedef ValuedRegionalMinimaImageFilter Self;

  typedef ValuedRegionalExtremaImageFilter< TInputImage, TOutputImage,
                                            std::less< typename TInputImage::PixelType >,
                                            std::less< typename TOutputImage::PixelType >  > Superclass;

  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer;

  typedef TInputImage                        InputImageType;
  typedef typename InputImageType::PixelType InputImagePixelType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Runtime information support. */
  itkTypeMacro(ValuedRegionalMinimaImageFilter,
               ValuedRegionalExtremaImageFilter);

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro( InputPixelTypeComparable,
                   ( Concept::LessThanComparable< InputImagePixelType > ) );
  itkConceptMacro( InputHasPixelTraitsCheck,
                   ( Concept::HasPixelTraits< InputImagePixelType > ) );
  itkConceptMacro( InputHasNumericTraitsCheck,
                   ( Concept::HasNumericTraits< InputImagePixelType > ) );
  /** End concept checking */
#endif
protected:
  ValuedRegionalMinimaImageFilter()
  {
    this->SetMarkerValue( NumericTraits< ITK_TYPENAME TOutputImage::PixelType >::max() );
  }

  virtual ~ValuedRegionalMinimaImageFilter() {}
private:
  ValuedRegionalMinimaImageFilter(const Self &); //purposely not implemented
  void operator=(const Self &);                  //purposely not implemented
};                                               // end
                                                 // ValuedRegionalMinimaImageFilter
} //end namespace itk
#endif
