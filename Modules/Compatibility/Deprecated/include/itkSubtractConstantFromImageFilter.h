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
#ifndef __itkSubtractConstantFromImageFilter_h
#define __itkSubtractConstantFromImageFilter_h
#include "itkSubtractImageFilter.h"

namespace itk
{

/** \class SubtractConstantFromImageFilter
 *
 * \brief Add a constant to all input pixels.
 *
 * This filter is templated over the input image type
 * and the output image type.
 *
 * \author Tom Vercauteren, INRIA & Mauna Kea Technologies
 *
 * Based on filters from the Insight Journal paper:
 * http://hdl.handle.net/1926/510
 *
 * \ingroup ITKDeprecated
 * \sa SubtractImageFilter
 */
template <class TInputImage, class TConstant, class TOutputImage>
class ITK_EXPORT SubtractConstantFromImageFilter :
      public
SubtractImageFilter<TInputImage, Image<TConstant, TInputImage::ImageDimension>, TOutputImage>
{
public:
  typedef SubtractConstantFromImageFilter                           Self;
  typedef SubtractImageFilter<TInputImage, Image<TConstant, TInputImage::ImageDimension>, TOutputImage>
                                                                    Superclass;
  typedef SmartPointer<Self>                                        Pointer;
  typedef SmartPointer<const Self>                                  ConstPointer;

  /** method for creation through object factory */
  itkNewMacro(Self);
  /** Run-time type information (and related methods). */
  itkTypeMacro(SubtractConstantFromImageFilter, SubtractImageFilter);

protected:
  SubtractConstantFromImageFilter() {}
  virtual ~SubtractConstantFromImageFilter() {}
};

}
#endif
