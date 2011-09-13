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
#ifndef __itkScalarToArrayCastImageFilter_h
#define __itkScalarToArrayCastImageFilter_h

#include "itkImageToImageFilter.h"

namespace itk
{
/** \class ScalarToArrayCastImageFilter
 *
 * \brief Creates the output image with vector type pixels filled with
 * the intensity values from one or more input images with scalar
 * pixels.
 *
 * This filter is templated over the input image type and
 * output image type. The each dimension of the output image pixel is
 * filled with each input image pixel's scalar pixel value. This
 * filter can be used to cast a scalar image to a vector image if
 * there is only one input image.
 *
 * \ingroup Deprecated
 * \ingroup ITKDeprecated
 */

template< class TInputImage, class TOutputImage >
class ITK_EXPORT ScalarToArrayCastImageFilter:
  public ImageToImageFilter< TInputImage, TOutputImage >
{
public:
  /** Standard class typedefs. */
  typedef ScalarToArrayCastImageFilter                    Self;
  typedef ImageToImageFilter< TInputImage, TOutputImage > Superclass;
  typedef SmartPointer< Self >                            Pointer;
  typedef SmartPointer< const Self >                      ConstPointer;

  /** Standard class macros */
  itkNewMacro(Self);
  itkTypeMacro(ScalarToArrayCastImageFilter, ImageToImageFilter);

  typedef typename Superclass::OutputImageRegionType OutputImageRegionType;
  typedef typename TOutputImage::PixelType           OutputImagePixelType;

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro( OutputHasNumericTraitsCheck,
                   ( Concept::HasNumericTraits< typename OutputImagePixelType::ValueType > ) );
  itkConceptMacro( OutputHasPixelTraitsCheck,
                   ( Concept::HasPixelTraits< OutputImagePixelType > ) );
  /** End concept checking */
#endif
protected:
  ScalarToArrayCastImageFilter();
  virtual ~ScalarToArrayCastImageFilter() {}

  void ThreadedGenerateData(const OutputImageRegionType & outputRegionForThread,
                            ThreadIdType threadId);

private:
  ScalarToArrayCastImageFilter(const Self &); //purposely not implemented
  void operator=(const Self &);               //purposely not implemented
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkScalarToArrayCastImageFilter.hxx"
#endif

#endif
