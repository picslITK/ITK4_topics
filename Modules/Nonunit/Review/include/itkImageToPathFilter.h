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
#ifndef __itkImageToPathFilter_h
#define __itkImageToPathFilter_h

#include "itkImage.h"
#include "itkPathSource.h"

namespace itk
{
/** \class ImageToPathFilter
 * \brief Base class for filters that take an image as input and produce an path as output.
 *
 * ImageToPathFilter is the base class for all process objects that output
 * path data and require image data as input. Specifically, this class
 * defines the SetInput() method for defining the input to a filter.
 *
 * \ingroup ImageFilters
 * \ingroup ITK-Review
 */
template< class TInputImage, class TOutputPath >
class ITK_EXPORT ImageToPathFilter:public PathSource< TOutputPath >
{
public:
  /** Standard class typedefs. */
  typedef ImageToPathFilter          Self;
  typedef PathSource< TOutputPath >  Superclass;
  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro(ImageToPathFilter, PathSource);

  /** Some convenient typedefs. */
  typedef TInputImage                           InputImageType;
  typedef typename InputImageType::Pointer      InputImagePointer;
  typedef typename InputImageType::ConstPointer InputImageConstPointer;
  typedef typename InputImageType::RegionType   InputImageRegionType;
  typedef typename InputImageType::PixelType    InputImagePixelType;

  /** ImageDimension constants */
  itkStaticConstMacro(InputImageDimension, unsigned int,
                      TInputImage::ImageDimension);

  /** Set/Get the image input of this process object.  */
  virtual void SetInput(const InputImageType *image);

  virtual void SetInput(unsigned int, const TInputImage *image);

  const InputImageType * GetInput(void);

  const InputImageType * GetInput(unsigned int idx);

protected:
  ImageToPathFilter();
  ~ImageToPathFilter();

  virtual void PrintSelf(std::ostream & os, Indent indent) const;

private:
  ImageToPathFilter(const Self &); //purposely not implemented
  void operator=(const Self &);    //purposely not implemented
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkImageToPathFilter.txx"
#endif

#endif
